#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region over 3 ecoregions
# 10-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    -habitat suitability, habitat area, current density layers
#	   - generic resistance layer
#    - protected areas layer
#	   - natural areas layer 
#	   - ecoregions layer
#  
#   Outputs: prioritization solution rasters for scenarios, intermediate files
#    
#   * Notes - the code is pretty long. 
#           1) Load workspace
#           2) Landscape data is loaded & prepped
#           3) Load and process species feature data
#           4) Set prioritization budgets and costs
#           5) Run prioritization scenarios and merge ecoregion-level solutions 
#             5.1) Generic resistance
#             5.2) Combine species Res, then calc Density
#             5.3) Combine species densities
#             5.4) Minimize shortfall
#             5.5) Maximize utility 
#           - ignore message "Discarded datum Unknown based on GRS80 ellipsoid in CRS 
#               definition, but +towgs84= values preserved"" 
#               This is due to sp/sf compatibility bug but doesn't affect output
#           - Protected areas are omitted from input maps and ignored in solutions    
#                                                                   
# Script by C Tucker for ApexRMS 									
#############################################################################################

## 1) Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(raster)
library(sp)
library(prioritizr)
library(vegan)
library(Rsymphony)
library(zonator)
library(corrplot)

dev.new()

setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

## Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Results/PrioritizationMaps"

## Functions 
rescaleR <- function(x, new.min = 0, new.max = 1) {
  x.min = suppressWarnings(min(x, na.rm=TRUE))
  x.max = suppressWarnings(max(x, na.rm=TRUE))
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

evaltext <- function(x, y){
  sapply(x, 
         FUN = function(X){
           eval(parse(text = paste0(X, y)))}
  )}

evaltext2 <- function(x, y){
  sapply(y, 
         FUN = function(Y){
           eval(parse(text = paste0(x, Y)))}
  )}


## 2) Load files ---------------------------------------------------------

# Load species list
speciesID <- read.csv(
  file.path(
    paste0(rawDataDir, "/Focal Species"), 
    "Species.csv"), 
  stringsAsFactors = FALSE)
specieslist <- speciesID$Code

# Focal area
LULC <- raster(file.path(procDataDir, "LULC_FocalArea.tif")) # 1465929 cells
naturalAreasFocal <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(procDataDir, "LULCbinary_FocalArea.tif"))

# Protected areas
protectedAreas <- raster(file.path(procDataDir, "protectedAreasNatural_FocalArea.tif"))
#values so 0 = protected, 1 = unprotected, reflecting cost of inclusion in network
protectedAreasNA <-  reclassify(protectedAreas, rcl=matrix(c(0, 1, 1, 0), ncol=2, byrow=T)) %>%
                       reclassify(., rcl=matrix(c(0, NA, 1, 1), ncol=2, byrow=T))

# Ecoregions 
ecoregions <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
  calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
  crop(., naturalAreasFocal) %>%
  mask(., naturalAreasFocal)
# Ecoregions - zone1 Adirondacks, zone 3 = StL lowlands, Zone 4 = appalachians

# Set up 3 ecoregion zones (1, 3 and 4)
zone1 <- calc(ecoregions, fun=function(x){ifelse(x==1, 1, NA)})
zone3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
zone4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})

naturalAreasBinaryFocal1 <-  mask(naturalAreasBinaryFocal, zone1)
naturalAreasBinaryFocal3 <-  mask(naturalAreasBinaryFocal, zone3)
naturalAreasBinaryFocal4 <-  mask(naturalAreasBinaryFocal, zone4)

protectedAreasNA1 <-  mask(protectedAreasNA, zone1)
protectedAreasNA3 <-  mask(protectedAreasNA, zone3)
protectedAreasNA4 <-  mask(protectedAreasNA, zone4)


## 3) Load and generate feature files ---------------------------------------------------------

## Using for-loop over zones, for i species with 3 input layers each
# all layers normalized
# log scaling for current density values only (many orders of magnitude var)
# Range scaling reqd for input into prioritizr

ecozones <- c(1, 3, 4)

for(j in ecozones){ # run for both zones, all species
  
  # ID  inputs for the correct zone
  zone <- eval(parse(text=paste0("zone", j))) 
  naturalAreasZ <- eval(parse(text=paste0("naturalAreasBinaryFocal", j))) 
  protectedAreasZ <- eval(parse(text=paste0("protectedAreasNA", j))) 
  
  ## For species-level inputs
  for(i in specieslist){ # run for all species
    species <- i
    
density <- raster(file.path(procDataDir, paste0(species, "_Resistance_FocalAreaBuffer_out_cum_curmap.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
      calc(., fun = log) %>% #density has log normal distbn
      scale(.) %>%
      calc(., fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam2 <- paste0(species, "_density")
assign(nam2, density)		
    
habitatSuitability <- raster(file.path(procDataDir, paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
      scale(.) %>%
      calc(., fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam3 <- paste0(species, "_habitatSuitability")
assign(nam3, habitatSuitability)
    
habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
      scale(.) %>%
      calc(., fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam4 <- paste0(species, "_habitatArea")
assign(nam4, habitatArea)
    
rm(density, habitatSuitability, habitatArea)
  
   } # Finish species data loop
  
# Combine species level inputs into stacks, by zone and feature type 
Suitability <- stack(evaltext(specieslist, "_habitatSuitability"))  
names(Suitability) <- specieslist
namm1 <- paste0("Suitability", j)
assign(namm1, Suitability)
  
Area <- stack(evaltext(specieslist, "_habitatArea"))  
names(Area) <- specieslist
namm2 <- paste0("Area", j)
assign(namm2, Area)
  
Density <- stack(evaltext(specieslist, "_density"))  
names(Density) <- specieslist
namm3 <- paste0("Density", j)
assign(namm3, Density)
  
# Combine for all features into complete stack
All <- c(evaltext(specieslist, "_habitatSuitability"), 
        evaltext(specieslist, "_habitatArea"),
        evaltext(specieslist, "_density"))
All <- stack(All)
names(All) <- c(paste0(specieslist, "_habitatSuitability"), 
                  paste0(specieslist, "_habitatArea"),
                  paste0(specieslist, "_density"))
namm4 <- paste0("All", j)
assign(namm4, All)
  
## For generic resistance layer
genericRes <- raster(file.path(procDataDir, paste0("Generic_Resistance_FocalAreaBuffer_out_cum_curmap.tif"))) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
namn3 <- paste0("genericRes", j)
assign(namn3, genericRes)
  
  
## For processed current density inputs
SumResDensity <- raster(file.path(procDataDir, "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("SumResDensity", j)
assign(nammm1, SumResDensity)
  
MeanResDensity <- raster(file.path(procDataDir, "combined_Resistance_Mean_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm2 <- paste0("MeanResDensity", j)
assign(nammm2, MeanResDensity)
  
MaxResDensity <- raster(file.path(procDataDir, "combined_Resistance_Max_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm3 <- paste0("MaxResDensity", j)
assign(nammm3, MaxResDensity)
  
rm(SumResDensity, MeanResDensity, MaxResDensity)
  
} # Finish ecoregion loop

# Output is: All3 and All4 (rasterStacks with all features), 
# Area3/Area4, Suitability3/Suitability 4. and Density3/4 (rasterStacks for types of features)
# Sum, Mean, max ResDensity layers, GenericRes layer


## 3.1) Combine transformed rasters into single map per feature type -------------------------------

# Density input from generic resistance
genericResAll <- mosaic(genericRes3, genericRes4, fun="max", na.rm=TRUE) %>%
                  mosaic(., genericRes1, fun="max", na.rm=TRUE)

# Density inputs from current density summarized as sum, mean, max resistances 
SumResDensity <- mosaic(SumResDensity3, SumResDensity4, fun="max", na.rm=TRUE) %>%
                    mosaic(., SumResDensity1, fun="max", na.rm=TRUE)
MeanResDensity <- mosaic(MeanResDensity3, MeanResDensity4, fun="max", na.rm=TRUE) %>%
                    mosaic(., MeanResDensity1, fun="max", na.rm=TRUE)
MaxResDensity <- mosaic(MaxResDensity3, MaxResDensity4, fun="max", na.rm=TRUE) %>%
                    mosaic(., MaxResDensity1, fun="max", na.rm=TRUE)
ResDensity <- stack(SumResDensity, MeanResDensity, MaxResDensity)
names(ResDensity) <- c("SumResDensity", "MeanResDensity", "MaxResDensity")

# Species level features summarized as sum, mean, or max
Suitability <- mosaic(Suitability3, Suitability4, fun="max", na.rm=TRUE) %>%
                  mosaic(., Suitability1, fun="max", na.rm=TRUE)
names(Suitability) <- specieslist

Area <- mosaic(Area3, Area4, fun="max", na.rm=TRUE)%>%
            mosaic(., Area1, fun="max", na.rm=TRUE)
names(Area) <- specieslist

Density <- mosaic(Density3, Density4, fun="max", na.rm=TRUE)%>%
              mosaic(., Density1, fun="max", na.rm=TRUE)
names(Density) <- specieslist

All <- mosaic(All3, All4, fun="max", na.rm=TRUE) %>%
          mosaic(., All1, fun="max", na.rm=TRUE)
names(All) <- c(paste0(specieslist, "_habitatSuitability"), 
                paste0(specieslist, "_habitatArea"),
                paste0(specieslist, "_density"))


## 3.2 Output intermediate rasters----------------------------------------------------------

writeRaster(Suitability, 
            file.path("Results",  "Suitability.tif"), 
            overwrite=TRUE)  
writeRaster(Density, 
            file.path("Results",  "Density.tif"), 
            overwrite=TRUE)  
writeRaster(Area, 
            file.path("Results",  "Area.tif"), 
            overwrite=TRUE)  
writeRaster(Area, 
            file.path("Results",  "ResDensity.tif"), 
            overwrite=TRUE)  
writeRaster(Area, 
            file.path("Results",  "genericResAll.tif"), 
            overwrite=TRUE)  

## 4) Set scenario target values ----------------------------------------------------------------

# Budget
Budget <- 0.20  ## proportion of the unprotected natural areas

# Identify all natural areas available for prioritization and omit PA areas
costLayer1 <- naturalAreasBinaryFocal1 %>%
  mask(., protectedAreasNA1, inv=TRUE) #omit protected area cells
costLayer3 <- naturalAreasBinaryFocal3 %>%
  mask(., protectedAreasNA3, inv=TRUE) #omit protected area cells
costLayer4 <- naturalAreasBinaryFocal4 %>%
  mask(., protectedAreasNA4, inv=TRUE) #omit protected area cells

# Ecoregion budgets, set by by ecoregion size (number of natural area pixels in zone)
NumSitesGoal1 <- round(Budget * cellStats(costLayer1, sum), 0)
NumSitesGoal3 <- round(Budget * cellStats(costLayer3, sum), 0)
NumSitesGoal4 <- round(Budget * cellStats(costLayer4, sum), 0)

#Set up solver for single layer objective problems

prob <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(z) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}


## 5) Prioritization scenarios, evaluating ecoregions separately------------------------------------

## 5.1) Scenario 1-Generic Resistance layer--------------------------------------------------------------------------------
# Front-end


genericResSol1 <- prob(costLayer1, genericRes1, NumSitesGoal1) %>% #input is the cost surface + features 
                  solve(.)

#constraints for following 2 are relaxed, hard to solve
genericResSol3 <- prob(costLayer3, genericRes3, NumSitesGoal3, first=T, gap=1.5) %>% #input is the cost surface + features 
                  solve(.)
genericResSol4 <- prob(costLayer4, genericRes4, NumSitesGoal4, first=T, gap=2) %>% #input is the cost surface + features 
                   solve(.)
  
# Final maps  per layer input type, ecozone, and statistic type
genericResSol <- mosaic(genericResSol3, genericResSol4, fun="max", na.rm=TRUE) %>%
                  mosaic(., genericResSol1, fun="max", na.rm=TRUE)

    
## 5.2) Scenario 2--Summarizes Resistance into Current Density ------------------------------------------------------------
# Current density layers calculated as the sum, mean, or max of species resistance layers
# Front-end
# Input = 1 layer, which is a summary of resistances to produce density

# Sum of resistances
fctSumResDensitySol1 <- prob(costLayer1, SumResDensity1, NumSitesGoal1) %>% #input is the cost surface + features 
                        solve(.)
#zone3
fctSumResDensitySol3 <- prob(costLayer3, SumResDensity3, NumSitesGoal3) %>% #input is the cost surface + features 
                        solve(.)
# zone 4
fctSumResDensitySol4 <- prob(costLayer4, SumResDensity4, NumSitesGoal3, gap=1.5) %>% #input is the cost surface + features 
                        solve(.)
# Final map, all ecoregions + protected areas
FinalSumResDensity <- mosaic(fctSumResDensitySol3, fctSumResDensitySol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., fctSumResDensitySol1, fun="max", na.rm=TRUE)

# Mean of resistances
fctMeanResDensitySol1 <- prob(costLayer1, MeanResDensity1, NumSitesGoal1) %>% #input is the cost surface + features 
                        solve(.)
#zone3
fctMeanResDensitySol3 <- prob(costLayer3, MeanResDensity3, NumSitesGoal3, gap=1.5)%>% #input is the cost surface + features 
                        solve(.)
# zone 4
fctMeanResDensitySol4 <- prob(costLayer4, MeanResDensity4, NumSitesGoal3, gap=1.5)%>% #input is the cost surface + features 
                        solve(.)
# Final map, all ecoregions + protected areas
FinalMeanResDensity <- mosaic(fctMeanResDensitySol3, fctMeanResDensitySol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., fctMeanResDensitySol1, fun="max", na.rm=TRUE)

# Max of resistances
fctMaxResDensitySol1 <- prob(costLayer1, MaxResDensity1, NumSitesGoal1) %>% #input is the cost surface + features 
                         solve(.)
#zone3
fctMaxResDensitySol3 <- prob(costLayer3, MaxResDensity3, NumSitesGoal3, gap=1.5)%>%  #input is the cost surface + features 
                        solve(.)
# zone 4
fctMaxResDensitySol4 <- prob(costLayer4, MaxResDensity4, NumSitesGoal3, gap=1.5)%>%  #input is the cost surface + features 
                          solve(.)
# Final map, all ecoregions + protected areas
FinalMaxResDensity <- mosaic(fctMaxResDensitySol3, fctMaxResDensitySol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., fctMaxResDensitySol1, fun="max", na.rm=TRUE)

      
## Scenario 5.3 - Summarize species densities --------------------------------------------------
# Take features from multiple species and feature types and summarize into a single layer. 
# Calculate  a summary value per pixel across all species and feature layers 
# Back end approach
      
statChoice <- c("sum", "max", "mean")
      
for(k in statChoice){ # loop over stat choices
        
    for(l in ecozones){   
          
zone <- eval(parse(text=paste0("zone", l))) 
costlayerZ <- eval(parse(text=paste0("costLayer", l))) 
NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
          
## Habitat Suitability only
inputs <- eval(parse(text=paste0("Suitability", l)))
fctSuitability <- stackApply(inputs, nlayers(inputs), k)
fctSuitabilitySol <- prob(costlayerZ, fctSuitability, NumSitesGoalZ) %>% #input is the cost surface + features 
                      solve(.)
assign(paste(k, "Suitability", l, sep="_"), fctSuitabilitySol)
          
## Density  only
inputs <- eval(parse(text=paste0("Density", l)))
fctDensity <- stackApply(inputs, nlayers(inputs), k)
fctDensitySol <- prob(costlayerZ, fctDensity, NumSitesGoalZ)%>%  #input is the cost surface + features 
                  solve(.)
assign(paste(k, "Density", l, sep="_"), fctDensitySol)
          
## Habitat area only
inputs <- eval(parse(text=paste0("Area", l)))
fctArea <- stackApply(inputs, nlayers(inputs), k)
fctAreaSol <- prob(costlayerZ, fctArea, NumSitesGoalZ) %>% #input is the cost surface + features 
                solve(.)
assign(paste(k, "Area", l, sep="_"), fctAreaSol)
          
## All
inputs <- eval(parse(text=paste0("All", l)))
fctAll <- stackApply(inputs, nlayers(inputs), k)
fctAllSol <- prob(costlayerZ, fctAll, NumSitesGoalZ) %>% #input is the cost surface + features 
                solve(.)
assign(paste(k, "All", l, sep="_"), fctAllSol)
          
rm(zone, costlayerZ, NumSitesGoalZ, inputs, fctSuitability, fctSuitabilitySol, 
      fctDensitySol, fctAreaSol, fctAllSol)  
          
      } # End ecozone loop 
        
# Final maps  per layer input type, ecozone, and statistic type
SuitabilitySol <- mosaic(eval(parse(text=paste0(k, "_Suitability_3"))), 
                          eval(parse(text=paste0(k, "_Suitability_4"))), 
                          fun="max", na.rm=TRUE) %>%
                    mosaic(., eval(parse(text=paste0(k, "_Suitability_1"))), 
                          fun="max", na.rm=TRUE)
assign(paste0(k, "_FinalSuitability"), SuitabilitySol)
        
DensitySol <- mosaic(eval(parse(text=paste0(k, "_Density_3"))), 
                             eval(parse(text=paste0(k, "_Density_4"))), 
                             fun="max", na.rm=TRUE) %>%
              mosaic(., eval(parse(text=paste0(k, "_Density_1"))), 
                              fun="max", na.rm=TRUE)
assign(paste0(k, "_FinalDensity"), DensitySol)
        
AreaSol <- mosaic(eval(parse(text=paste0(k, "_Area_3"))), 
                  eval(parse(text=paste0(k, "_Area_4"))), 
                  fun="max", na.rm=TRUE) %>%
            mosaic(., eval(parse(text=paste0(k, "_Area_1"))), 
                         fun="max", na.rm=TRUE)
assign(paste0(k, "_FinalArea"), AreaSol)
        
AllSol <- mosaic(eval(parse(text=paste0(k, "_All_3"))), 
                eval(parse(text=paste0(k, "_All_4"))), 
                fun="max", na.rm=TRUE) %>%
          mosaic(., eval(parse(text=paste0(k, "_All_1"))),
                         fun="max", na.rm=TRUE)
assign(paste0(k, "_FinalAll"), AllSol)
        
  } #end stat choice loop     
      
      
## Scenario 4 -  using minimize_shortfall_objective in prioritizer----------------------------------------------
#Back end
# Using all input layer, apply a multiobjective algorithm - minimize shortfall
    
    ## Suitability
    # zone 3
    minShortSuitProb3 <- problem(costLayer3, Suitability3) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
      add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(minShortSuitProb3) # < 30s
    minShortSuitSol3 <- solve(minShortSuitProb3)
    # zone 4
    minShortSuitProb4 <- problem(costLayer4, Suitability4) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
      add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(minShortSuitProb4) # < 30s
    minShortSuitSol4 <- solve(minShortSuitProb4)
    # Final map,  all ecoregions + protected areas
    minShortSuitSol <- mosaic(minShortSuitSol3, minShortSuitSol4, fun="max", na.rm=TRUE)
    minShortSuitSolPA <- mosaic(minShortSuitSol, protectedAreasNA, fun="max", na.rm=TRUE)
    minShort_FinalSuitability <- stack(minShortSuitSol, minShortSuitSolPA)
    
    ## Density  
    # zone 3
    minShortdensityProb3 <- problem(costLayer3, Density3) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
      add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(minShortdensityProb3) # < 30s
    minShortdensitySol3 <- solve(minShortdensityProb3)
    # zone 4
    minShortdensityProb4 <- problem(costLayer4, Density4) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
      add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(minShortdensityProb4) # < 30s
    minShortdensitySol4 <- solve(minShortdensityProb4)
    # Final map, all ecoregions + protected areas
    minShortDensitySol <- mosaic(minShortdensitySol3, minShortdensitySol4, fun="max", na.rm=TRUE)
    minShortDensitySolPA <- mosaic(minShortDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
    minShort_FinalDensity <- stack(minShortDensitySol, minShortDensitySolPA)
    
    ## Area
    # zone 3
    minShortAreaProb3 <- problem(costLayer3, Area3) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
      add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(minShortAreaProb3) # < 30s
    minShortAreaSol3 <- solve(minShortAreaProb3)
    # zone 4
    minShortAreaProb4 <- problem(costLayer4, Area4) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
      add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(minShortAreaProb4) # < 30s
    minShortAreaSol4 <- solve(minShortAreaProb4)
    # Final map, all ecoregions + protected areas
    minShortAreaSol <- mosaic(minShortAreaSol3, minShortAreaSol4, fun="max", na.rm=TRUE)
    minShortAreaSolPA <- mosaic(minShortAreaSol, protectedAreasNA, fun="max", na.rm=TRUE)
    minShort_FinalArea <- stack(minShortAreaSol, minShortAreaSolPA)
    
    ## All
    # zone 3
    minShortProb3 <- problem(costLayer3, All3) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
      add_relative_targets(0.50) %>%  # set high to avoid being a constraint??
      add_binary_decisions() %>% #inclusion vs no-inclusion	
      add_default_solver()
    presolve_check(minShortProb3)
    minShortSol3 <- solve(minShortProb3)
    # zone 4
    minShortProb4 <- problem(costLayer4, All4) %>% #input is the cost surface + features 
      add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
      add_relative_targets(0.50) %>%  # set high to avoid being a constraint??
      add_binary_decisions() %>% #inclusion vs no-inclusion	
      add_default_solver()
    presolve_check(minShortProb4)
    minShortSol4 <- solve(minShortProb4)
    # Final map,  all ecoregions + protected areas
    minShortAllSol <- mosaic(minShortSol3, minShortSol4, fun="max", na.rm=TRUE)
    minShortAllSolPA <- mosaic(minShortAllSol, protectedAreasNA, fun="max", na.rm=TRUE)
    minShort_FinalAll <- stack(minShortAllSol, minShortAllSolPA)
    
    
    ## Scenario 4 -  using add_max_utility objective in prioritizer----------------------------------------------
    
    ## Suitability
    # zone 3
    maxUtilitySuitProb3 <- problem(costLayer3, Suitability3) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(maxUtilitySuitProb3) # < 30s
    maxUtilitySuitSol3 <- solve(maxUtilitySuitProb3)
    # zone 4
    maxUtilitySuitProb4 <- problem(costLayer4, Suitability4) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(maxUtilitySuitProb4) # < 30s
    maxUtilitySuitSol4 <- solve(maxUtilitySuitProb4)
    # Final map,  all ecoregions + protected areas
    maxUtilitySuitSol <- mosaic(maxUtilitySuitSol3, maxUtilitySuitSol4, fun="max", na.rm=TRUE)
    maxUtilitySuitSolPA <- mosaic(maxUtilitySuitSol, protectedAreasNA, fun="max", na.rm=TRUE)
    maxUtility_FinalSuitability <- stack(maxUtilitySuitSol, maxUtilitySuitSolPA)
    
    ## Density  
    # zone 3
    maxUtilitydensityProb3 <- problem(costLayer3, Density3) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(maxUtilitydensityProb3) # < 30s
    maxUtilitydensitySol3 <- solve(maxUtilitydensityProb3)
    # zone 4
    maxUtilitydensityProb4 <- problem(costLayer4, Density4) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(maxUtilitydensityProb4) # < 30s
    maxUtilitydensitySol4 <- solve(maxUtilitydensityProb4)
    # Final map,  all ecoregions + protected areas
    maxUtilityDensitySol <- mosaic(maxUtilitydensitySol3, maxUtilitydensitySol4, fun="max", na.rm=TRUE)
    maxUtilityDensitySolPA <- mosaic(maxUtilityDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
    maxUtility_FinalDensity <- stack(maxUtilityDensitySol, maxUtilityDensitySolPA)
    
    ## Area
    # zone 3
    maxUtilityAreaProb3 <- problem(costLayer3, Area3) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(maxUtilityAreaProb3) # < 30s
    maxUtilityAreaSol3 <- solve(maxUtilityAreaProb3)
    # zone 4
    maxUtilityAreaProb4 <- problem(costLayer4, Area4) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
      add_binary_decisions() #inclusion vs no-inclusion	
    presolve_check(maxUtilityAreaProb4) # < 30s
    maxUtilityAreaSol4 <- solve(maxUtilityAreaProb4)
    # Final map, all ecoregions + protected areas
    maxUtilityAreaSol <- mosaic(maxUtilityAreaSol3, maxUtilityAreaSol4, fun="max", na.rm=TRUE)
    maxUtilityAreaSolPA <- mosaic(maxUtilityAreaSol, protectedAreasNA, fun="max", na.rm=TRUE)
    maxUtility_FinalArea <- stack(maxUtilityAreaSol, maxUtilityAreaSolPA)
    
    ## All
    # zone 3
    maxUtilityProb3 <- problem(costLayer3, All3) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
      add_binary_decisions() %>% #inclusion vs no-inclusion	
      add_default_solver()
    presolve_check(maxUtilityProb3)
    maxUtilitySol3 <- solve(maxUtilityProb3)
    # zone 4
    maxUtilityProb4 <- problem(costLayer4, All4) %>% #input is the cost surface + features 
      add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
      add_binary_decisions() %>% #inclusion vs no-inclusion	
      add_default_solver()
    presolve_check(maxUtilityProb4)
    maxUtilitySol4 <- solve(maxUtilityProb4)
    # Final map,  all ecoregions + protected areas
    maxUtilityAllSol <- mosaic(maxUtilitySol3, maxUtilitySol4, fun="max", na.rm=TRUE)
    maxUtilityAllSolPA <- mosaic(maxUtilityAllSol, protectedAreasNA, fun="max", na.rm=TRUE)
    maxUtility_FinalAll <- stack(maxUtilityAllSol, maxUtilityAllSolPA)
    
    
    ## Combine and summarize output files --------------------------------------------------------------------
    
    numModels <- 23
    
    # Raster stack all solns incl PA
    outputAll <- stack(sum_FinalAll[[1]],
                       sum_FinalSuitability[[1]], 
                       sum_FinalDensity[[1]], 
                       sum_FinalArea[[1]],
                       max_FinalAll[[1]],
                       max_FinalSuitability[[1]],
                       max_FinalDensity[[1]],
                       max_FinalArea[[1]],
                       mean_FinalAll[[1]],
                       mean_FinalSuitability[[1]],
                       mean_FinalDensity[[1]],
                       mean_FinalArea[[1]],
                       FinalSumResDensity[[1]],
                       FinalMeanResDensity[[1]],
                       FinalMaxResDensity[[1]],
                       minShort_FinalAll[[1]],
                       minShort_FinalSuitability[[1]],
                       minShort_FinalDensity[[1]],  
                       minShort_FinalArea[[1]],
                       maxUtility_FinalAll[[1]],
                       maxUtility_FinalSuitability[[1]],
                       maxUtility_FinalDensity[[1]],  
                       maxUtility_FinalArea[[1]])
    
    mapNames <- c("SumAll",
                  "SumSuit",
                  "SumDensity",
                  "SumArea",
                  "MaxAll",
                  "MaxSuit",
                  "MaxDensity",
                  "MaxArea",
                  "MeanAll",
                  "MeanSuit",
                  "MeanDensity",
                  "MeanArea",
                  "SumResistDensity",
                  "MeanResistDensity",
                  "MaxResistDensity",
                  "MSAll",
                  "MSSuit",
                  "MSDensity",
                  "MSArea",
                  "MUAll",
                  "MUSuit",
                  "MUDensity",
                  "MUArea")
    names(outputAll) <- mapNames
    outputAll <- mask(outputAll, Suitability[[1]])
    
    # Collect layers with PA added
    outputAllPA <-  stack(sum_FinalAll[[2]],
                          sum_FinalSuitability[[2]], 
                          sum_FinalDensity[[2]], 
                          sum_FinalArea[[2]],
                          max_FinalAll[[2]],
                          max_FinalSuitability[[2]],
                          max_FinalDensity[[2]],
                          max_FinalArea[[2]],
                          mean_FinalAll[[2]],
                          mean_FinalSuitability[[2]],
                          mean_FinalDensity[[2]],
                          mean_FinalArea[[2]],
                          FinalSumResDensity[[2]],
                          FinalMeanResDensity[[2]],
                          FinalMaxResDensity[[2]],
                          minShort_FinalAll[[2]],
                          minShort_FinalSuitability[[2]],
                          minShort_FinalDensity[[2]],  
                          minShort_FinalArea[[2]],
                          maxUtility_FinalAll[[2]],
                          maxUtility_FinalSuitability[[2]],
                          maxUtility_FinalDensity[[2]],  
                          maxUtility_FinalArea[[2]])
    
    mapNamesPA <- paste0(c("SumAll",
                           "SumSuit",
                           "SumDensity",
                           "SumArea",
                           "MaxAll",
                           "MaxSuit",
                           "MaxDensity",
                           "MaxArea",
                           "MeanAll",
                           "MeanSuit",
                           "MeanDensity",
                           "MeanArea",
                           "SumResistDensity",
                           "MeanResistDensity",
                           "MaxResistDensity",
                           "MSAll",
                           "MSSuit",
                           "MSDensity",
                           "MSArea",
                           "MUAll",
                           "MUSuit",
                           "MUDensity",
                           "MUArea"), "_PA")
    names(outputAllPA) <- mapNamesPA
    outputAllPA <- mask(outputAllPA, Suitability[[1]])
    
    
    
    ## Export solution rasters ----------------------------------------------------------------------
    
    
 
    
    # output results maps
    writeRaster(outputAll, 
                filename=file.path(outDir, 
                                   paste0(names(outputAll), "_prioritizationMap.tif")), 
                bylayer=TRUE,
                overwrite=TRUE)
    
    writeRaster(outputAllPA, 
                filename=file.path(outDir, 
                                   paste0(names(outputAll), "_prioritizationMap.tif")), 
                bylayer=TRUE,
                overwrite=TRUE)
    
    ## End script ---------------------------------------
    