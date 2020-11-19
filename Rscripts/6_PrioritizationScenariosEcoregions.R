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
#           2) Landscape data is loaded & processed
#           3) Load and process species feature data
#           4) Process and output intermediate data
#           5) Set prioritization budgets and costs
#           6) Run prioritization scenarios and merge ecoregion-level solutions 
#             6.1) Generic resistance
#             6.2) Combine species Res, then calc Density
#             6.3) Combine species densities
#             6.4) Minimize shortfall
#             6.5) Maximize utility 
#           7) Combine and outputsolution files
#     * ignore message "Discarded datum Unknown based on GRS80 ellipsoid in CRS 
#       definition, but +towgs84= values preserved".
#       This is due to sp/sf compatibility bug but doesn't affect output
#     * Protected areas are omitted from input maps and ignored in solutions    
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

  # Using for-loop over ecoregion zones, and over i species
  # all layers normalized
  # log scaling for current density values only (many orders of magnitude var)
  # Range scaling reqd for input into prioritizr

ecozones <- c(1, 3, 4)

for(j in ecozones){ # multispecies inputs
  
  # ID  inputs for the correct zone
  zone <- eval(parse(text=paste0("zone", j))) 
  naturalAreasZ <- eval(parse(text=paste0("naturalAreasBinaryFocal", j))) 
  protectedAreasZ <- eval(parse(text=paste0("protectedAreasNA", j))) 
  
## 3a) For species-level inputs-----------------------------------------------------------------
  for(i in specieslist){ 
    species <- i

density <- raster(file.path(procDataDir, "CircuitscapeOutputs", paste0(species, "_Resistance_FocalAreaBuffer_out_cum_curmap.tif"))) %>%
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

    if(i=="MAAM"&&j=="1"){
     
       # Replace MAAM Area raster due to single appropriate patch in zone1
      MAAM_habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
        crop(., naturalAreasZ) %>%
        mask(., naturalAreasZ) %>%
        mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
        calc(., fun=function(x){ifelse(x==271, 1, NA)})
      
      }else{
    
habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
      scale(.) %>%
      calc(., fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam4 <- paste0(species, "_habitatArea")
assign(nam4, habitatArea)
      }

rm(density, habitatSuitability, habitatArea)
  
   } # Finish species data loop

  
## 3b) For generic resistance layer-----------------------------------------------------------------
genericRes <- raster(file.path(procDataDir, "CircuitscapeOutputs", paste0("Generic_Resistance_FocalAreaBuffer_out_cum_curmap.tif"))) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
namn3 <- paste0("genericRes", j)
assign(namn3, genericRes)
  
  
## 3c) For processed current density inputs----------------------------------------------------------
SumResDensity <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("SumResDensity", j)
assign(nammm1, SumResDensity)
  
MeanResDensity <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Mean_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm2 <- paste0("MeanResDensity", j)
assign(nammm2, MeanResDensity)
  
MaxResDensity <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Max_FocalAreaBuffer_out_cum_curmap.tif")) %>%
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
  

## 3d) For ecoprofile current density inputs-----------------------------------------------------
  # Scenario 1
ecoprofileBird <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofile1Bird", j)
assign(nammm1, ecoprofileBird)

ecoprofileMammal <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Mean_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm2 <- paste0("ecoprofile1Mammal", j)
assign(nammm2, ecoprofileBird)

ecoprofileAmphibian <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Max_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm3 <- paste0("ecoprofile1Amphibian", j)
assign(nammm3, ecoprofileAmphibian)

  # Scenario2
ecoprofileOmnivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileOmnivore", j)
assign(nammm1, ecoprofile2Omnivore)

ecoprofileInsectivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileInsectivore", j)
assign(nammm1, ecoprofile2Insectivore)

ecoprofileCarnivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileCarnivore", j)
assign(nammm1, ecoprofile2Carnivore)

ecoprofileHerbivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileHerbivore", j)
assign(nammm1, ecoprofile2Herbivore)


## Collect intermediate data

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

    } # Finish ecoregion loop

  # Output is: All3 and All4 (rasterStacks with all features), 
  # Area3/Area4, Suitability3/Suitability 4. and Density3/4 (rasterStacks for types of features)
  # Sum, Mean, max ResDensity layers, GenericRes layer, ecoprofile layers


## 4.1) Combine transformed rasters into single map per feature type -------------------------------

  # Density input from generic resistance
genericResAll <- mosaic(genericRes3, genericRes4, fun="max", na.rm=TRUE) %>%
                  mosaic(., genericRes1, fun="max", na.rm=TRUE)

  # Density inputs from current density (already summarized as sum, mean, max resistances) 
SumResDensity <- mosaic(SumResDensity3, SumResDensity4, fun="max", na.rm=TRUE) %>%
                    mosaic(., SumResDensity1, fun="max", na.rm=TRUE)
MeanResDensity <- mosaic(MeanResDensity3, MeanResDensity4, fun="max", na.rm=TRUE) %>%
                    mosaic(., MeanResDensity1, fun="max", na.rm=TRUE)
MaxResDensity <- mosaic(MaxResDensity3, MaxResDensity4, fun="max", na.rm=TRUE) %>%
                    mosaic(., MaxResDensity1, fun="max", na.rm=TRUE)
ResDensity <- stack(SumResDensity, MeanResDensity, MaxResDensity)
names(ResDensity) <- c("SumResDensity", "MeanResDensity", "MaxResDensity")

  # Species level features (to be summarized as sum, mean, or max)
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


## 4.2) Output intermediate rasters----------------------------------------------------------

writeRaster(Suitability, 
            file.path(outDir,  "Suitability.tif"), 
            overwrite=TRUE)  
writeRaster(Density, 
            file.path(outDir,  "Density.tif"), 
            overwrite=TRUE)  
writeRaster(Area, 
            file.path(outDir,  "Area.tif"), 
            overwrite=TRUE)  
writeRaster(ResDensity, 
            file.path(outDir,  "ResDensity.tif"), 
            overwrite=TRUE)  
writeRaster(genericResAll, 
            file.path(outDir,  "genericResAll.tif"), 
            overwrite=TRUE)  

## 5) Set scenario target values ----------------------------------------------------------------

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


## 6) Prioritization scenarios, evaluating ecoregions separately------------------------------------

  # Set up solver for single layer objective problems
prob <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(z) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}

## 6.1) Generic Resistance layer--------------------------------------------------------------------------------

  # Solver struggles with these so I loosen the gap size and return first identified soln
genericResSol1 <- prob(costLayer1, genericRes1, NumSitesGoal1) %>% #input is the cost surface + features 
                      solve(.)
  # constraints for following 2 are relaxed, hard to solve
genericResSol3 <- prob(costLayer3, genericRes3, NumSitesGoal3, first=T, gap=10) %>% #input is the cost surface + features 
                      solve(.)
genericResSol4 <- prob(costLayer4, genericRes4, NumSitesGoal4, first=T, gap=5) %>% #input is the cost surface + features 
                    solve(.)
  # Final maps  per layer input type, ecozone, and statistic type
genericResSol <- mosaic(genericResSol3, genericResSol4, fun="max", na.rm=TRUE) %>%
                  mosaic(., genericResSol1, fun="max", na.rm=TRUE)

    
## 6.2) Summarizes Resistance into Current Density ------------------------------------------------------------

  # Current density layers calculated as the sum, mean, or max of species resistance layers
  # Input = 1 layer, which is a summary of resistances to produce density

  # Sum of resistances
fctSumResDensitySol1 <- prob(costLayer1, SumResDensity1, NumSitesGoal1) %>% #input is the cost surface + features 
                        solve(.)
  # zone3
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
fctMeanResDensitySol3 <- prob(costLayer3, MeanResDensity3, NumSitesGoal3, gap=1.5)%>% #input is the cost surface + features 
                        solve(.)
fctMeanResDensitySol4 <- prob(costLayer4, MeanResDensity4, NumSitesGoal3, gap=1.5)%>% #input is the cost surface + features 
                        solve(.)
FinalMeanResDensity <- mosaic(fctMeanResDensitySol3, fctMeanResDensitySol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., fctMeanResDensitySol1, fun="max", na.rm=TRUE)

  # Max of resistances
fctMaxResDensitySol1 <- prob(costLayer1, MaxResDensity1, NumSitesGoal1) %>% #input is the cost surface + features 
                         solve(.)
fctMaxResDensitySol3 <- prob(costLayer3, MaxResDensity3, NumSitesGoal3, gap=1.5)%>%  #input is the cost surface + features 
                        solve(.)
fctMaxResDensitySol4 <- prob(costLayer4, MaxResDensity4, NumSitesGoal3, gap=3) %>%  #input is the cost surface + features 
                          solve(.)
FinalMaxResDensity <- mosaic(fctMaxResDensitySol3, fctMaxResDensitySol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., fctMaxResDensitySol1, fun="max", na.rm=TRUE)

      
## 6.3) Summarize species densities --------------------------------------------------

  # Take features from multiple species and feature types and summarize into a single layer. 
  # Calculate  a summary value per pixel across all species and feature layers 

statChoice <- c("sum", "max", "mean")
      
for(k in statChoice){ # loop over stat choices
        
    for(l in ecozones){   
          
zone <- eval(parse(text=paste0("zone", l))) 
costlayerZ <- eval(parse(text=paste0("costLayer", l))) 
NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
          
  # Habitat Suitability only
inputs <- eval(parse(text=paste0("Suitability", l)))
fctSuitability <- stackApply(inputs, nlayers(inputs), k)
fctSuitabilitySol <- prob(costlayerZ, fctSuitability, NumSitesGoalZ) %>% #input is the cost surface + features 
                      solve(.)
assign(paste(k, "Suitability", l, sep="_"), fctSuitabilitySol)
          
  # Density  only
inputs <- eval(parse(text=paste0("Density", l)))
fctDensity <- stackApply(inputs, nlayers(inputs), k)
fctDensitySol <- prob(costlayerZ, fctDensity, NumSitesGoalZ)%>%  #input is the cost surface + features 
                  solve(.)
assign(paste(k, "Density", l, sep="_"), fctDensitySol)
          
  # Habitat area only
inputs <- eval(parse(text=paste0("Area", l)))
fctArea <- stackApply(inputs, nlayers(inputs), k)
fctAreaSol <- prob(costlayerZ, fctArea, NumSitesGoalZ) %>% #input is the cost surface + features 
                solve(.)
assign(paste(k, "Area", l, sep="_"), fctAreaSol)
          
  # All
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
      
      
## 6.4) Using minimize_shortfall_objective in prioritizer--------------------------------

  # Using all input layer, apply a multiobjective algorithm - minimize shortfall
    
  # Set up minimum shortfall solver for single layer objective problems
probMS <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_min_shortfall_objective(z) %>% #minimize cost surface
    add_relative_targets(0.5) %>%
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}


  # Suitability
minShortSuitSol1 <- probMS(costLayer1, Suitability1, NumSitesGoal1) %>%
                      solve(.)
minShortSuitSol3 <- probMS(costLayer3, Suitability3, NumSitesGoal3) %>%
                      solve(.)
minShortSuitSol4 <- probMS(costLayer4, Suitability4, NumSitesGoal4) %>%
                     solve(.)
minShortSuitSol <- mosaic(minShortSuitSol3, minShortSuitSol4, fun="max", na.rm=TRUE) %>%
                      mosaic(., minShortSuitSol1, fun="max", na.rm=TRUE)

  # Density  
minShortDensitySol1 <- probMS(costLayer1, Density1, NumSitesGoal1) %>%
                        solve(.)
minShortDensitySol3 <- probMS(costLayer3, Density3, NumSitesGoal3) %>%
                        solve(.)
minShortDensitySol4 <- probMS(costLayer4, Density4, NumSitesGoal4) %>%
                        solve(.)
minShortDensitySol <- mosaic(minShortDensitySol3, minShortDensitySol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., minShortDensitySol1, fun="max", na.rm=TRUE)

  # Area
minShortAreaSol1 <- probMS(costLayer1, Area1, NumSitesGoal1) %>%
                      solve(.)
minShortAreaSol3 <- probMS(costLayer3, Area3, NumSitesGoal3, gap=2) %>%
                      solve(.)
minShortAreaSol4 <- probMS(costLayer4, Area4, NumSitesGoal4) %>%
                      solve(.)
minShortAreaSol <- mosaic(minShortAreaSol3, minShortAreaSol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., minShortAreaSol1, fun="max", na.rm=TRUE)

  # All
minShortAllSol1 <- probMS(costLayer1, All1, NumSitesGoal1, gap=2) %>%
                    solve(.)
minShortAllSol3 <- probMS(costLayer3, All3, NumSitesGoal3, gap=2) %>%
                    solve(.)
minShortAllSol4 <- probMS(costLayer4, All4, NumSitesGoal4, gap=2) %>%
                    solve(.)
minShortAllSol <- mosaic(minShortAllSol3, minShortAllSol4, fun="max", na.rm=TRUE) %>%
                    mosaic(., minShortAllSol1, fun="max", na.rm=TRUE)    
    

## 6.5) Scenario using max_utility_objective in prioritizer--------------------------------

# Using all input layer, apply a multiobjective algorithm - max utility

probMU <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(z) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}
   
  # Suitability
maxUtilitySuitSol1 <- probMU(costLayer1, Suitability1, NumSitesGoal1) %>%
                      solve(.)
maxUtilitySuitSol3 <- probMU(costLayer3, Suitability3, NumSitesGoal3) %>%
                      solve(.)
maxUtilitySuitSol4 <- probMU(costLayer4, Suitability4, NumSitesGoal4) %>%
                      solve(.)
maxUtilitySuitSol <- mosaic(maxUtilitySuitSol3, maxUtilitySuitSol4, fun="max", na.rm=TRUE) %>%
                      mosaic(., maxUtilitySuitSol1, fun="max", na.rm=TRUE)

  # Density  
maxUtilityDensitySol1 <- probMU(costLayer1, Density1, NumSitesGoal1) %>%
                          solve(.)
maxUtilityDensitySol3 <- probMU(costLayer3, Density3, NumSitesGoal3) %>%
                          solve(.)
maxUtilityDensitySol4 <- probMU(costLayer4, Density4, NumSitesGoal4) %>%
                          solve(.)
maxUtilityDensitySol <- mosaic(maxUtilityDensitySol3, maxUtilityDensitySol4, fun="max", na.rm=TRUE) %>%
                          mosaic(., maxUtilityDensitySol1, fun="max", na.rm=TRUE)

  # Area
maxUtilityAreaSol1 <- probMU(costLayer1, Area1, NumSitesGoal1, gap=2) %>%
                        solve(.)
maxUtilityAreaSol3 <- probMU(costLayer3, Area3, NumSitesGoal3) %>%
                        solve(.)
maxUtilityAreaSol4 <- probMU(costLayer4, Area4, NumSitesGoal4) %>%
                        solve(.)
maxUtilityAreaSol <- mosaic(maxUtilityAreaSol3, maxUtilityAreaSol4, fun="max", na.rm=TRUE) %>%
                        mosaic(., maxUtilityAreaSol1, fun="max", na.rm=TRUE)

  # All
maxUtilityAllSol1 <- probMU(costLayer1, All1, NumSitesGoal1, gap=2) %>%
                        solve(.)
maxUtilityAllSol3 <- probMU(costLayer3, All3, NumSitesGoal3) %>%
                         solve(.)
maxUtilityAllSol4 <- probMU(costLayer4, All4, NumSitesGoal4) %>%
                          solve(.)
maxUtilityAllSol <- mosaic(maxUtilityAllSol3, maxUtilityAllSol4, fun="max", na.rm=TRUE) %>%
                          mosaic(., maxUtilityAllSol1, fun="max", na.rm=TRUE)    

    
## 7 Combine and summarize output files --------------------------------------------------------------------

  # Raster stack all solns 
outputAll <- stack(genericResSol,
                   FinalSumResDensity,
                   FinalMeanResDensity,
                   FinalMaxResDensity,
                   sum_FinalSuitability,
                   sum_FinalDensity,
                   sum_FinalArea,
                   sum_FinalAll,
                   mean_FinalSuitability,
                   mean_FinalDensity,
                   mean_FinalArea,
                   mean_FinalAll,
                   max_FinalSuitability,
                   max_FinalDensity,
                   max_FinalArea,
                   max_FinalAll,
                   minShortSuitSol,
                   minShortDensitySol,
                   minShortAreaSol,
                   minShortAllSol,
                   maxUtilitySuitSol,
                   maxUtilityDensitySol,
                   maxUtilityAreaSol,
                   maxUtilityAllSol)
    
    mapNames <- c("genericResSol",
                  "FinalSumResDensity",
                  "FinalMeanResDensity",
                  "FinalMaxResDensity",
                  "sum_FinalSuitability",
                  "sum_FinalDensity",
                  "sum_FinalArea",
                  "sum_FinalAll",
                  "mean_FinalSuitability",
                  "mean_FinalDensity",
                  "mean_FinalArea",
                  "mean_FinalAll",
                  "max_FinalSuitability",
                  "max_FinalDensity",
                  "max_FinalArea",
                  "max_FinalAll",
                  "minShortSuitSol",
                  "minShortDensitySol",
                  "minShortAreaSol",
                  "minShortAllSol",
                  "maxUtilitySuitSol",
                  "maxUtilityDensitySol",
                  "maxUtilityAreaSol",
                  "maxUtilityAllSol")
    names(outputAll) <- mapNames
    outputAll <- mask(outputAll, Suitability[[1]])
    
    
  # Export solution rasters 
writeRaster(outputAll, 
              filename=file.path(outDir, 
                                   paste0(names(outputAll), "_prioritizationMap.tif")), 
                bylayer=TRUE,
                overwrite=TRUE)
   

## End script ---------------------------------------
    