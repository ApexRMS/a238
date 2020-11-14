#####################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region
# 10-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    -habitat suitability, habitat area, current density layers
#	- protected areas layer
#	- natural areas layer 
#	- ecoregions layer
#  
#   Outputs: prioritization solution rasters for scenarios
#    
#   * Notes - the code is really long.
#           - ignore s w "Discarded datum Unknown based on GRS80 ellipsoid in CRS definition,
#               but +towgs84= values preserved"" Due to sp/sf incompatibility but doesn't affect
#               data.
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

## Workspace ---------------------------------------------------------

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

#setwd("~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity/")
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


## Load files ---------------------------------------------------------

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
protectedAreasBinary <-  reclassify(protectedAreas, rcl=matrix(c(0, 1, 1, 0), ncol=2, byrow=T))
#values so 1 = protected, 0 = unprotected, reflecting contribution to network
protectedAreasNA <-  reclassify(protectedAreasBinary, rcl=matrix(c(0, NA, 1, 1), ncol=2, byrow=T))

  # Ecoregions
ecoregions <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
              calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
              crop(., naturalAreasFocal) %>%
              mask(., naturalAreasFocal) %>%
              calc(., fun=function(x){ifelse(x==1, 3, x)}) #(zone 1 is combined with zone 3)

  # Set up 2 ecoregion zones (3 and 4)
zone3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
zone4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})

naturalAreasBinaryFocal3 <-  mask(naturalAreasBinaryFocal, zone3)
naturalAreasBinaryFocal4 <-  mask(naturalAreasBinaryFocal, zone4)

protectedAreasNA3 <-  mask(protectedAreasNA, zone3)
protectedAreasNA4 <-  mask(protectedAreasNA, zone4)
protectedAreas3 <- mask(protectedAreasBinary, zone3)
protectedAreas4 <- mask(protectedAreasBinary, zone4)

## Generate feature files ---------------------------------------------------------

## Using for loop over j zones, for i species with 3 input layers each
  # all layers normalized
  # log scaling for current density values (many orders of magnitude var)
  # Range scaling for input into prioritizr


for(j in c(3, 4)){ # run for both zones, all species

zone <- eval(parse(text=paste0("zone", j))) 
naturalAreasZ <- eval(parse(text=paste0("naturalAreasBinaryFocal", j))) 

  for(i in specieslist){ # run for all species
  
      species <- i
        
density <- raster(file.path(procDataDir, paste0(species, "_Resistance_FocalAreaBuffer_out_cum_curmap.tif"))) %>%
        crop(., naturalAreasZ) %>%
        mask(., naturalAreasZ) %>%
        calc(., fun = log) %>% #density has log normal distbn
				scale(.) %>%
				calc(., fun = rescaleR) %>%
        calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam2 <- paste0(species, "_density")
assign(nam2, density)		

habitatSuitability <- raster(file.path(procDataDir, paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      scale(.) %>%
      calc(., fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam3 <- paste0(species, "_habitatSuitability")
assign(nam3, habitatSuitability)

habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam4 <- paste0(species, "_habitatArea")
assign(nam4, habitatArea)

rm(density, habitatSuitability, habitatArea)
	} # Finish species data loop
  #Ignore warnings

## Combine into stacks, by zone and feature type 
Suitability <- stack(evaltext(specieslist, "_habitatSuitability"))  
names(Suitability) <- specieslist
namm <- paste0("Suitability", j)
assign(namm, Suitability)

Area <- stack(evaltext(specieslist, "_habitatArea"))  
names(Area) <- specieslist
namm <- paste0("Area", j)
assign(namm, Area)

Density <- stack(evaltext(specieslist, "_density"))  
names(Density) <- specieslist
namm <- paste0("Density", j)
assign(namm, Density)

## Combine for all features into complete stack
All <- c(evaltext(specieslist, "_habitatSuitability"), 
                      evaltext(specieslist, "_habitatArea"),
                      evaltext(specieslist, "_density"))
All <- stack(All)

names(All) <- c(paste0(specieslist, "_habitatSuitability"), 
                             paste0(specieslist, "_habitatArea"),
                             paste0(specieslist, "_density"))
namm <- paste0("All", j)
assign(namm, All)

  # Ignore warnings (sp/sf proj issue)
} # Finish ecoregion loop

  # Output is: All3 and All4 (rasterStacks with all features), 
  # Area3/Area4, Suitability3/Suitability 4. and Density3/4 (rasterStacks for types of features)


## Load and process summary current density files for scenario 2------------------------------------------------
  # Current density calculated on summarized (sum, mean, max) resistance layer

for(p in c(3, 4)){ # run for both zones, all species
  
  zone <- eval(parse(text=paste0("zone", p))) 
  naturalAreasZ <- eval(parse(text=paste0("naturalAreasBinaryFocal", p))) 

SumResDensity <- raster(file.path(procDataDir, "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
namm <- paste0("SumResDensity", p)
assign(namm, SumResDensity)

MeanResDensity <- raster(file.path(procDataDir, "combined_Resistance_Mean_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
namm <- paste0("MeanResDensity", p)
assign(namm, MeanResDensity)

MaxResDensity <- raster(file.path(procDataDir, "combined_Resistance_Max_FocalAreaBuffer_out_cum_curmap.tif")) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  calc(., fun = log) %>% #density has log normal distbn
  scale(.) %>%
  calc(., fun = rescaleR) %>%
  calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
namm <- paste0("MaxResDensity", p)
assign(namm, MaxResDensity)

} 


## Merge ecoregion rasters to create complete rasters covering entire Monteregie ---------------

  # For calculating representation later
Suitability <- merge(Suitability3, Suitability4)
names(Suitability) <- specieslist

Area <- merge(Area3, Area4)
names(Area) <- specieslist

Density <- merge(Density3, Density4)
names(Density) <- specieslist

All <- merge(All3, All4)
names(All) <- c(paste0(specieslist, "_habitatSuitability"), 
                paste0(specieslist, "_habitatArea"),
                paste0(specieslist, "_density"))


## Set target values ----------------------------------------------------------------

# Budget
Budget <- 0.20  
costLayer3 <- naturalAreasBinaryFocal3
costLayer4 <- naturalAreasBinaryFocal4

# Ecoregion budgets, set by by ecoregion size (number of natural area pixels in zone)
NumSitesGoal3 <- mask(costLayer3, protectedAreasNA3, inv=TRUE)
NumSitesGoal3 <- round(Budget * cellStats(NumSitesGoal3, sum), 0)
NumSitesGoal4 <- mask(costLayer4, protectedAreasNA4, inv=TRUE)
NumSitesGoal4 <- round(Budget * cellStats(NumSitesGoal4, sum), 0)

###################################################################################################
## Prioritization scenarios evaluating ecoregions separately--------------------------------------

#Set up solvers for zones

prob3 <- function(x, y){
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
    add_default_solver()
}

prob4 <- function(x, y){
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
    add_default_solver()
}


## Scenario 1 - take features from multiple species and feature types and summarize into a single layer. 
  # Back end approach  
  # Calculate  a summary value per pixel across all species and feature layers 
  # Choice of potential summary metrics - sum, median, max

statChoice <- c("sum", "max", "mean")


for(k in statChoice){ # loop over stat choices
    
## Habitat Suitability only
fctSuitability3 <- stackApply(Suitability3, nlayers(Suitability3), k)
fctSuitabilityProb3 <- prob3(costLayer3, fctSuitability3) #input is the cost surface + features 
fctSuitabilitySol3 <- solve(fctSuitabilityProb3)
# zone 4
fctSuitability4 <- stackApply(Suitability4, nlayers(Suitability4), k)
fctSuitabilityProb4 <- prob4(costLayer4, fctSuitability4) #input is the cost surface + features 
fctSuitabilitySol4 <- solve(fctSuitabilityProb4)
# Final map, sum suit, all ecoregions + protected areas
fctSuitabilitySol <- mosaic(fctSuitabilitySol3, fctSuitabilitySol4, fun="max", na.rm=TRUE) 
fctSuitabilitySolPA <- mosaic(fctSuitabilitySol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctSuitabilitySol, fctSuitabilitySolPA)
assign(paste0(k, "_FinalSuitability"), final)
rm(final, fctSuitability3, fctSuitability4, fctSuitabilitySol, fctSuitabilityProb3, fctSuitabilityProb4, fctSuitabilitySolPA)

## Density only 
fctDensity3 <- stackApply(Density3, nlayers(Density3), k)
fctDensityProb3 <- prob3(costLayer3, fctDensity3) #input is the cost surface + features 
fctDensitySol3 <- solve(fctDensityProb3)
# zone 4
fctDensity4 <- stackApply(Density4, nlayers(Density4), k)
fctDensityProb4 <- prob4(costLayer4, fctDensity4) #input is the cost surface + features 
fctDensitySol4 <- solve(fctDensityProb4)
# Final map, sum suit, all ecoregions + protected areas
fctDensitySol <- mosaic(fctDensitySol3, fctDensitySol4, fun="max", na.rm=TRUE) 
fctDensitySolPA <- mosaic(fctDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctDensitySol, fctDensitySolPA)
assign(paste0(k, "_FinalDensity"), final)
rm(final, fctDensity3, fctDensity4, fctDensitySol, fctDensityProb3, fctDensityProb4, fctDensitySolPA)

## Habitat area only
fctArea3 <- stackApply(Area3, nlayers(Area3), k)
fctAreaProb3 <- prob3(costLayer3, fctArea3) #input is the cost surface + features 
fctAreaSol3 <- solve(fctAreaProb3)
# zone 4
fctArea4 <- stackApply(Area4, nlayers(Area4), k)
fctAreaProb4 <- prob4(costLayer4, fctArea4) #input is the cost surface + features 
fctAreaSol4 <- solve(fctAreaProb4)
# Final map, sum suit, all ecoregions + protected areas
fctAreaSol <- mosaic(fctAreaSol3, fctAreaSol4, fun="max", na.rm=TRUE) 
fctAreaSolPA <- mosaic(fctAreaSol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctAreaSol, fctAreaSolPA)
assign(paste0(k, "_FinalArea"), final)
rm(final, fctArea3, fctArea4, fctAreaSol, fctAreaProb3, fctAreaProb4, fctAreaSolPA)

## All
fctAll3 <- stackApply(All3, nlayers(All3), k)
fctAllProb3 <- prob3(costLayer3, fctAll3) #input is the cost surface + features 
fctAllSol3 <- solve(fctAllProb3)
# zone 4
fctAll4 <- stackApply(All4, nlayers(All4), k)
fctAllProb4 <- prob4(costLayer4, fctAll4) #input is the cost surface + features 
fctAllSol4 <- solve(fctAllProb4)
# Final map, sum suit, all ecoregions + protected areas
fctAllSol <- mosaic(fctAllSol3, fctAllSol4, fun="max", na.rm=TRUE) 
fctAllSolPA <- mosaic(fctAllSol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctAllSol, fctAllSolPA)
assign(paste0(k, "_FinalAll"), final)
rm(final, fctAll3, fctAll4, fctAllSol, fctAllProb3, fctAllProb4, fctAllSolPA)

} # End statistic choice loop  

  # Outputs - k_FinalAll stacks with a solution, and solution w PA added in

## Scenario 2--------------------------------------------------------------------------------
  # Front-end
  # Input = 1 layer, which is a summary of resistances to produce density
  # Summarization via 3 methods, mean, max, sum

  # Sum of resistances
fctSumResDensityProb3 <- prob3(costLayer3, SumResDensity3) #input is the cost surface + features 
fctSumResDensitySol3 <- solve(fctSumResDensityProb3)
# zone 4
fctSumResDensityProb4 <- prob4(costLayer4, SumResDensity4) #input is the cost surface + features 
fctSumResDensitySol4 <- solve(fctSumResDensityProb4)
# Final map, sum suit, all ecoregions + protected areas
fctSumResDensitySol <- mosaic(fctSumResDensitySol3, fctSumResDensitySol4, fun="max", na.rm=TRUE) 
fctSumResDensitySolPA <- mosaic(fctSumResDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctSumResDensitySol, fctSumResDensitySolPA)
assign("FinalSumResDensity", final)
rm(final)

  # Mean of resistances
fctMeanResDensityProb3 <- prob3(costLayer3, MeanResDensity3) #input is the cost surface + features 
fctMeanResDensitySol3 <- solve(fctMeanResDensityProb3)
# zone 4
fctMeanResDensityProb4 <- prob4(costLayer4, MeanResDensity4) #input is the cost surface + features 
fctMeanResDensitySol4 <- solve(fctMeanResDensityProb4)
# Final map, sum suit, all ecoregions + protected areas
fctMeanResDensitySol <- mosaic(fctMeanResDensitySol3, fctMeanResDensitySol4, fun="max", na.rm=TRUE) 
fctMeanResDensitySolPA <- mosaic(fctMeanResDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctMeanResDensitySol, fctMeanResDensitySolPA)
assign("FinalMeanResDensity", final)
rm(final)

  # Max of resistances
fctMaxResDensityProb3 <- prob3(costLayer3, MaxResDensity3) #input is the cost surface + features 
fctMaxResDensitySol3 <- solve(fctMaxResDensityProb3)
# zone 4
fctMaxResDensityProb4 <- prob4(costLayer4, MaxResDensity4) #input is the cost surface + features 
fctMaxResDensitySol4 <- solve(fctMaxResDensityProb4)
# Final map, sum suit, all ecoregions + protected areas
fctMaxResDensitySol <- mosaic(fctMaxResDensitySol3, fctMaxResDensitySol4, fun="max", na.rm=TRUE) 
fctMaxResDensitySolPA <- mosaic(fctMaxResDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctMaxResDensitySol, fctMaxResDensitySolPA)
assign("FinalMaxResDensity", final)
rm(final)
  # Outputs - statResDensitySolAll containing solution and solution + PA


## Scenario 3 -  using minimize_shortfall_objective in prioritizer----------------------------------------------
  #Back end
  # Using all input layer, apply a multiobjective algorithm
  # minimize shortfall

## Suitability
  # zone 3
minShortSuitProb3 <- problem(costLayer3, Suitability3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(minShortSuitProb3) # < 30s
minShortSuitSol3 <- solve(minShortSuitProb3)
  # zone 4
minShortSuitProb4 <- problem(costLayer4, Suitability4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(minShortSuitProb4) # < 30s
minShortSuitSol4 <- solve(minShortSuitProb4)
  # Final map, sum suit, all ecoregions + protected areas
minShortSuitSol <- mosaic(minShortSuitSol3, minShortSuitSol4, fun="max", na.rm=TRUE)
minShortSuitSolPA <- mosaic(minShortSuitSol, protectedAreasNA, fun="max", na.rm=TRUE)
minShort_FinalSuitability <- stack(minShortSuitSol, minShortSuitSolPA)

## Density  
  # zone 3
minShortdensityProb3 <- problem(costLayer3, Density3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(minShortdensityProb3) # < 30s
minShortdensitySol3 <- solve(minShortdensityProb3)
  # zone 4
minShortdensityProb4 <- problem(costLayer4, Density4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(minShortdensityProb4) # < 30s
minShortdensitySol4 <- solve(minShortdensityProb4)
  # Final map, sum suit, all ecoregions + protected areas
minShortDensitySol <- mosaic(minShortdensitySol3, minShortdensitySol4, fun="max", na.rm=TRUE)
minShortDensitySolPA <- mosaic(minShortDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
minShort_FinalDensity <- stack(minShortDensitySol, minShortDensitySolPA)

## Area
  # zone 3
minShortAreaProb3 <- problem(costLayer3, Area3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(minShortAreaProb3) # < 30s
minShortAreaSol3 <- solve(minShortAreaProb3)
  # zone 4
minShortAreaProb4 <- problem(costLayer4, Area4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.5) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(minShortAreaProb4) # < 30s
minShortAreaSol4 <- solve(minShortAreaProb4)
# Final map, sum suit, all ecoregions + protected areas
minShortAreaSol <- mosaic(minShortAreaSol3, minShortAreaSol4, fun="max", na.rm=TRUE)
minShortAreaSolPA <- mosaic(minShortAreaSol, protectedAreasNA, fun="max", na.rm=TRUE)
minShort_FinalArea <- stack(minShortAreaSol, minShortAreaSolPA)

## All
 # zone 3
minShortProb3 <- problem(costLayer3, All3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.50) %>%  # set high to avoid being a constraint??
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(minShortProb3)
minShortSol3 <- solve(minShortProb3)
 # zone 4
minShortProb4 <- problem(costLayer4, All4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.50) %>%  # set high to avoid being a constraint??
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(minShortProb4)
minShortSol4 <- solve(minShortProb4)
  # Final map, sum suit, all ecoregions + protected areas
minShortAllSol <- mosaic(minShortSol3, minShortSol4, fun="max", na.rm=TRUE)
minShortAllSolPA <- mosaic(minShortAllSol, protectedAreasNA, fun="max", na.rm=TRUE)
minShort_FinalAll <- stack(minShortAllSol, minShortAllSolPA)


## Scenario 4 -  using add_max_utility objective in prioritizer----------------------------------------------

## Suitability
# zone 3
maxUtilitySuitProb3 <- problem(costLayer3, Suitability3) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(maxUtilitySuitProb3) # < 30s
maxUtilitySuitSol3 <- solve(maxUtilitySuitProb3)
# zone 4
maxUtilitySuitProb4 <- problem(costLayer4, Suitability4) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(maxUtilitySuitProb4) # < 30s
maxUtilitySuitSol4 <- solve(maxUtilitySuitProb4)
# Final map, sum suit, all ecoregions + protected areas
maxUtilitySuitSol <- mosaic(maxUtilitySuitSol3, maxUtilitySuitSol4, fun="max", na.rm=TRUE)
maxUtilitySuitSolPA <- mosaic(maxUtilitySuitSol, protectedAreasNA, fun="max", na.rm=TRUE)
maxUtility_FinalSuitability <- stack(maxUtilitySuitSol, maxUtilitySuitSolPA)

## Density  
# zone 3
maxUtilitydensityProb3 <- problem(costLayer3, Density3) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(maxUtilitydensityProb3) # < 30s
maxUtilitydensitySol3 <- solve(maxUtilitydensityProb3)
# zone 4
maxUtilitydensityProb4 <- problem(costLayer4, Density4) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(maxUtilitydensityProb4) # < 30s
maxUtilitydensitySol4 <- solve(maxUtilitydensityProb4)
# Final map, sum suit, all ecoregions + protected areas
maxUtilityDensitySol <- mosaic(maxUtilitydensitySol3, maxUtilitydensitySol4, fun="max", na.rm=TRUE)
maxUtilityDensitySolPA <- mosaic(maxUtilityDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
maxUtility_FinalDensity <- stack(maxUtilityDensitySol, maxUtilityDensitySolPA)

## Area
# zone 3
maxUtilityAreaProb3 <- problem(costLayer3, Area3) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(maxUtilityAreaProb3) # < 30s
maxUtilityAreaSol3 <- solve(maxUtilityAreaProb3)
# zone 4
maxUtilityAreaProb4 <- problem(costLayer4, Area4) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(maxUtilityAreaProb4) # < 30s
maxUtilityAreaSol4 <- solve(maxUtilityAreaProb4)
# Final map, sum suit, all ecoregions + protected areas
maxUtilityAreaSol <- mosaic(maxUtilityAreaSol3, maxUtilityAreaSol4, fun="max", na.rm=TRUE)
maxUtilityAreaSolPA <- mosaic(maxUtilityAreaSol, protectedAreasNA, fun="max", na.rm=TRUE)
maxUtility_FinalArea <- stack(maxUtilityAreaSol, maxUtilityAreaSolPA)

## All
# zone 3
maxUtilityProb3 <- problem(costLayer3, All3) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal3) %>% #minimize cost surface
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(maxUtilityProb3)
maxUtilitySol3 <- solve(maxUtilityProb3)
# zone 4
maxUtilityProb4 <- problem(costLayer4, All4) %>% #input is the cost surface + features 
  add_max_utility_objective(NumSitesGoal4) %>% #minimize cost surface
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(maxUtilityProb4)
maxUtilitySol4 <- solve(maxUtilityProb4)
# Final map, sum suit, all ecoregions + protected areas
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


  # Output intermediate rasters
writeRaster(Suitability, filename=file.path("Results", "Suitability.tif"), 
            format="raster", 
            overwrite=TRUE)
writeRaster(Density, filename=file.path("Results", "Density.tif"), 
            format="raster", 
            overwrite=TRUE) 
writeRaster(Area, filename=file.path("Results", "Area.tif"), 
            format="raster", 
            overwrite=TRUE)
writeRaster(All, filename=file.path("Results", "All.tif"), 
            format="raster", 
            overwrite=TRUE)

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
