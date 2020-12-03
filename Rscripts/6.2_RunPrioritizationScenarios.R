#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Run prioritization scenarioswith variable targets
# 10-2020                                       					
#
#	  Inputs:                                
#    - Load input rasters from 6.1_PrepareInputsPrioritizationScenario
#    - protected areas layer
#	   - natural areas layer 
#	   - ecoregions layer
# 
#   Outputs:
#     - solution maps for 21 scenarios
#   
#   Structure:  
#               1) Workspace
#               2) Load input files
#               3) Set prioritization costs and targets
#     -----------Set budget here (Line 197)------------
#               4) Run prioritization scenarios and merge ecoregion-level solutions 
#                     4.1) Generic resistance
#                     4.2) Ecoprofiles
#                     4.3) Combine species Res, then calc Density
#                     4.4) Combine species densities
#                     4.5) Multi objective, Minimize shortfall
#                     4.6) Multi objective, Maximize utility 
#               5) Combine and output solution files
#     
#    * ignore message "Discarded datum Unknown based on GRS80 ellipsoid in CRS 
#       definition, but +towgs84= values preserved".
#       This is due to sp/sf compatibility issue but doesn't affect output
#     * Protected areas are omitted from input maps and ignored in solutions    
#     * Periodically a scenario will give you a message of "infeasible solution returned, try relaxing 
#       parameters." Generally you can just re-run the specific code and it will converge
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

setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

options(warn = -1)

## Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
inputDir <- "Data/Processed/PrioritizationInputs"
outDir <- "Results/PrioritizationSolutions"

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


## Load prioritizR solver function
source("a238/Rscripts/prioritzRSolverFunctions.R")


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
LULCterrestrial <- calc(LULC, fun=function(x){ifelse(x==700, NA, x)})
#sizeLULC <- 1465929
naturalAreasFocal <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(procDataDir, "LULCbinary_FocalArea.tif"))
#sizenatural <- 450287

# Protected areas
protectedAreas <- raster(file.path(procDataDir, "protectedAreasTerrestrial_FocalArea.tif"))%>%
  crop(., naturalAreasFocal)

# Ecoregions - zone1 Adirondacks, zone 4 = StL lowlands, Zone 3 = appalachians
ecoregions <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
  calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
  crop(., naturalAreasFocal) %>%
  mask(., naturalAreasFocal)
#Eco regions all terrestrial areas
ecoregionsLULC <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
  calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
  crop(., LULCterrestrial) %>%
  mask(., LULCterrestrial)

# Divide landscapes into 3 ecoregion zones (1, 3 and 4)
zone1 <- calc(ecoregions, fun=function(x){ifelse(x==1, 1, NA)})
zone3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
zone4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})
#
naturalAreasBinaryFocal1 <-  mask(naturalAreasBinaryFocal, zone1)
naturalAreasBinaryFocal3 <-  mask(naturalAreasBinaryFocal, zone3)
naturalAreasBinaryFocal4 <-  mask(naturalAreasBinaryFocal, zone4)
#
protectedAreasNA1 <-  mask(protectedAreas, zone1)
protectedAreasNA3 <-  mask(protectedAreas, zone3)
protectedAreasNA4 <-  mask(protectedAreas, zone4)
#

# Calculate zone-specific % protected areas for deliverables (as # of cells)
zone1LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==1, 1, NA)})
zone3LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==3, 1, NA)})
zone4LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==4, 1, NA)})
#
zonePA1 <- cellStats(protectedAreasNA1, sum, na.rm=T)#/cellStats(zone1LULC, sum, na.rm=T)
zonePA3 <- cellStats(protectedAreasNA3, sum, na.rm=T)#/cellStats(zone3LULC, sum, na.rm=T)
zonePA4 <- cellStats(protectedAreasNA4, sum, na.rm=T)#/cellStats(zone4LULC, sum, na.rm=T)

# Load outputs from Script 6.1-------------------------------------------------
genericRes <- stack(file.path(inputDir,  "genericResAll.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
ResDensity <- stack(file.path(inputDir,  "ResDensity.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(ResDensity) <- c("SumResDensity", 
                       "MeanResDensity")

ResTaxon <- stack(file.path(inputDir,  "ResTaxon.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(ResTaxon) <- c("ecoprofileBird", 
                     "ecoprofileMammal", 
                     "ecoprofileAmphibian")

ResTrophic <- stack(file.path(inputDir,  "ResTrophic.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(ResTrophic) <- c("ecoprofileOmnivore", 
                       "ecoprofileInsectivore", 
                       "ecoprofileCarnivore", 
                       "ecoprofileHerbivore")

Suitability <- stack(file.path(inputDir,  "SuitabilityR.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
Density <- stack(file.path(inputDir,  "DensityR.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
Area <- stack(file.path(inputDir,  "AreaR.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(Suitability) <- names(Density) <- names(Area) <- specieslist

All <- stack(file.path(inputDir,  "AllR.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 

## Process input layers for specific scenarios------------------------------------------
# Sum  of species Layers
SuitabilityS <- stackApply(Suitability, nlayers(Suitability), "sum", na.rm=TRUE) 
DensityS <- stackApply(Density, nlayers(Density), "sum", na.rm=TRUE) 
AreaS <- stackApply(Area, nlayers(Area), "sum", na.rm=TRUE) 
AllS <-  stackApply(All, nlayers(All), "sum", na.rm=TRUE) 


# Mean  of species Layers
SuitabilityM <- stackApply(Suitability, nlayers(Suitability), "mean", na.rm=TRUE) 
DensityM <- stackApply(Density, nlayers(Density), "mean", na.rm=TRUE) 
AreaM <- stackApply(Area, nlayers(Area), "mean", na.rm=TRUE) 
AllM <-  stackApply(All, nlayers(All), "mean", na.rm=TRUE) 


## 3) Set prioritization targets  ----------------------------------------------------------------


## Budget - Change budget as desired ##
Budget <- 0.1 

# Identify all natural areas available for prioritization and omit PA areas
costLayer1 <- naturalAreasBinaryFocal1 %>%
  mask(., protectedAreasNA1, inv=TRUE) #omit protected area cells
costLayer3 <- naturalAreasBinaryFocal3 %>%
  mask(., protectedAreasNA3, inv=TRUE) #omit protected area cells
costLayer4 <- naturalAreasBinaryFocal4 %>%
  mask(., protectedAreasNA4, inv=TRUE) #omit protected area cells

# Calculate ecoregion budgets, set by by ecoregion size (number of natural area pixels in zone)
NumSitesGoal1 <- round((Budget * cellStats(zone1LULC, sum, na.rm=T)), 0) - zonePA1
NumSitesGoal1 <- ifelse(NumSitesGoal1 <= 0, 1, NumSitesGoal1)
NumSitesGoal3 <- round((Budget * cellStats(zone3LULC, sum, na.rm=T)), 0) - zonePA3
NumSitesGoal3 <- ifelse(NumSitesGoal3 <= 0, 1, NumSitesGoal3)
NumSitesGoal4 <- round((Budget * cellStats(zone4LULC, sum, na.rm=T)), 0) - zonePA4
NumSitesGoal4 <- ifelse(NumSitesGoal4 <= 0, 1, NumSitesGoal4)


## 4) Prioritization scenarios, evaluating ecoregions separately------------------------------------

## 4.1) Generic Resistance layer--------------------------------------------------------------------------------

  genericResSol1 <- genericSolver(genericRes, costLayer1, NumSitesGoal1)
  genericResSol3 <- genericSolver(genericRes, costLayer3, NumSitesGoal3)
  genericResSol4 <- genericSolver(genericRes, costLayer4, NumSitesGoal4)
  # Final map generic resistance
  genericResSol <- mosaic(genericResSol3, genericResSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., genericResSol1, fun="max", na.rm=TRUE)
 

## 4.2) Ecoprofile layers--------------------------------------------------------------------------------

# Ecoprofile Scenario 1
  resTaxonCombSol1 <-  stackApply(ResTaxon, nlayers(ResTaxon), "mean", na.rm=TRUE) %>%
    calc(., fun = rescaleR) %>%
    genericSolver(., costLayer1, NumSitesGoal1)
  resTaxonCombSol3 <-   stackApply(ResTaxon, nlayers(ResTaxon), "mean", na.rm=TRUE) %>%
    calc(., fun = rescaleR) %>%
    genericSolver(., costLayer3, NumSitesGoal3)
  resTaxonCombSol4 <- stackApply(ResTaxon, nlayers(ResTaxon), "mean", na.rm=TRUE) %>%
    calc(., fun = rescaleR) %>%
    genericSolver(., costLayer4, NumSitesGoal4)
  FinalEcoprofileTaxon <- mosaic(resTaxonCombSol3, resTaxonCombSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., resTaxonCombSol1, fun="max", na.rm=TRUE)


# Ecoprofile Scenario 2
  resTrophicCombSol1 <- stackApply(ResTrophic, nlayers(ResTaxon), "mean", na.rm=TRUE) %>%
    calc(., fun = rescaleR) %>%
    genericSolver(., costLayer1, NumSitesGoal1)
  resTrophicCombSol3 <-   stackApply(ResTrophic, nlayers(ResTaxon), "mean", na.rm=TRUE) %>%
    calc(., fun = rescaleR) %>%
    genericSolver(., costLayer3, NumSitesGoal3)
  resTrophicCombSol4 <- stackApply(ResTrophic, nlayers(ResTaxon), "mean", na.rm=TRUE) %>%
    calc(., fun = rescaleR) %>%
    genericSolver(., costLayer4, NumSitesGoal4)
  FinalEcoprofileTrophic <- mosaic(resTrophicCombSol3, resTrophicCombSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., resTrophicCombSol1, fun="max", na.rm=TRUE)


## 4.3) Summarizes Resistance into Current Density ------------------------------------------------------------

# Current density layers calculated as the sum, mean of species resistance layers
# Input = 1 layer, which is a summary of resistances to produce density
  SumResDensity1 <- genericSolver(ResDensity[["SumResDensity"]], costLayer1, NumSitesGoal1)
  SumResDensity3 <- genericSolver(ResDensity[["SumResDensity"]], costLayer3, NumSitesGoal3)
  SumResDensity4 <- genericSolver(ResDensity[["SumResDensity"]], costLayer4, NumSitesGoal4)
  # Final map generic resistance
  FinalSumResDensity <- mosaic(SumResDensity3, SumResDensity4, fun="max", na.rm=TRUE) %>%
    mosaic(., SumResDensity1, fun="max", na.rm=TRUE)

  MeanResDensity1 <- genericSolver(ResDensity[["MeanResDensity"]], costLayer1, NumSitesGoal1)
  MeanResDensity3 <- genericSolver(ResDensity[["MeanResDensity"]], costLayer3, NumSitesGoal3)
  MeanResDensity4 <- genericSolver(ResDensity[["MeanResDensity"]], costLayer4, NumSitesGoal4)
  # Final map generic resistance
  FinalMeanResDensity <- mosaic(MeanResDensity3, MeanResDensity4, fun="max", na.rm=TRUE) %>%
    mosaic(., MeanResDensity1, fun="max", na.rm=TRUE)


## 4.4) Summarize species densities --------------------------------------------------

# Take current density features from multiple species and feature types and summarize into a single layer. 
# Calculate  a summary value per pixel across all species and feature layers 

# Habitat Suitability only
  SuitabilitySol1 <- genericSolver2(SuitabilityS, costLayer1, NumSitesGoal1)
  SuitabilitySol3 <- genericSolver2(SuitabilityS, costLayer3, NumSitesGoal3)
  SuitabilitySol4 <- genericSolver2(SuitabilityS, costLayer4, NumSitesGoal4)
  sum_FinalSuitability <- mosaic(SuitabilitySol3, SuitabilitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., SuitabilitySol1,  fun="max", na.rm=TRUE)

# Density  only
  DensitySol1 <- genericSolver2(DensityS, costLayer1, NumSitesGoal1)
  DensitySol3 <- genericSolver2(DensityS, costLayer3, NumSitesGoal3)
  DensitySol4 <- genericSolver2(DensityS, costLayer4, NumSitesGoal4)
  sum_FinalDensity <- mosaic(DensitySol3, DensitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., DensitySol1,  fun="max", na.rm=TRUE)

# Habitat area only
  AreaSol1 <- genericSolver2(AreaS, costLayer1, NumSitesGoal1)
  AreaSol3 <- genericSolver2(AreaS, costLayer3, NumSitesGoal3)
  AreaSol4 <- genericSolver2(AreaS, costLayer4, NumSitesGoal4)
  sum_FinalArea <- mosaic(AreaSol3, AreaSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., AreaSol1,  fun="max", na.rm=TRUE)

# All
  AllSol1 <- genericSolver2(AllS, costLayer1, NumSitesGoal1)
  AllSol3 <- genericSolver2(AllS, costLayer3, NumSitesGoal3)
  AllSol4 <- genericSolver2(AllS, costLayer4, NumSitesGoal4)
  sum_FinalAll <- mosaic(AllSol3, AllSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., AllSol1,  fun="max", na.rm=TRUE)

## Mean
# Habitat suitability
  SuitabilitySol1 <- genericSolver2(SuitabilityM, costLayer1, NumSitesGoal1)
  SuitabilitySol3 <- genericSolver2(SuitabilityM, costLayer3, NumSitesGoal3)
  SuitabilitySol4 <- genericSolver2(SuitabilityM, costLayer4, NumSitesGoal4)
  mean_FinalSuitability <- mosaic(SuitabilitySol3, SuitabilitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., SuitabilitySol1,  fun="max", na.rm=TRUE)

# Density  only
  DensitySol1 <- genericSolver2(DensityM, costLayer1, NumSitesGoal1)
  DensitySol3 <- genericSolver2(DensityM, costLayer3, NumSitesGoal3)
  DensitySol4 <- genericSolver2(DensityM, costLayer4, NumSitesGoal4)
  mean_FinalDensity <- mosaic(DensitySol3, DensitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., DensitySol1,  fun="max", na.rm=TRUE)

# Habitat area only
  AreaSol1 <- genericSolver2(AreaM, costLayer1, NumSitesGoal1)
  AreaSol3 <- genericSolver2(AreaM, costLayer3, NumSitesGoal3)
  AreaSol4 <- genericSolver2(AreaM, costLayer4, NumSitesGoal4)
  mean_FinalArea <- mosaic(AreaSol3, AreaSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., AreaSol1,  fun="max", na.rm=TRUE)

# All
  AllSol1 <- genericSolver2(AllM, costLayer1, NumSitesGoal1)
  AllSol3 <- genericSolver2(AllM, costLayer3, NumSitesGoal3)
  AllSol4 <- genericSolver2(AllM, costLayer4, NumSitesGoal4)
  mean_FinalAll <- mosaic(AllSol3, AllSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., AllSol1,  fun="max", na.rm=TRUE)


## 4.5) Using minimize_shortfall_objective in prioritizer--------------------------------
# Uses MSSolver
# Suitability
  MSSuitabilitySol1 <- MSSolver(Suitability, costLayer1, NumSitesGoal1)
  MSSuitabilitySol3 <- MSSolver(Suitability, costLayer3, NumSitesGoal3)
  MSSuitabilitySol4 <- MSSolver(Suitability, costLayer4, NumSitesGoal4)
  minShortSuitSol <- mosaic(MSSuitabilitySol3, MSSuitabilitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MSSuitabilitySol1,  fun="max", na.rm=TRUE)

# Density  only
  MSDensitySol1 <- MSSolver(Density, costLayer1, NumSitesGoal1)
  MSDensitySol3 <- MSSolver(Density, costLayer3, NumSitesGoal3)
  MSDensitySol4 <- MSSolver(Density, costLayer4, NumSitesGoal4)
  minShortDensitySol <- mosaic(MSDensitySol3, MSDensitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MSDensitySol1,  fun="max", na.rm=TRUE)  

# Habitat area only
  MSAreaSol1 <- MSSolver(Area, costLayer1, NumSitesGoal1)
  MSAreaSol3 <- MSSolver(Area, costLayer3, NumSitesGoal3)
  MSAreaSol4 <- MSSolver(Area, costLayer4, NumSitesGoal4)
  minShortAreaSol <- mosaic(MSAreaSol3, MSAreaSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MSAreaSol1,  fun="max", na.rm=TRUE)

#All
  MSAllSol1 <- MSSolver(All, costLayer1, NumSitesGoal1)
  MSAllSol3 <- MSSolver(All, costLayer3, NumSitesGoal3)
  MSAllSol4 <- MSSolver(All, costLayer4, NumSitesGoal4)
  minShortAllSol <- mosaic(MSAllSol3, MSAllSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MSAllSol1,  fun="max", na.rm=TRUE)    


## 4.6) Scenario using max_utility_objective in prioritizer--------------------------------

# Uses MUSolver
# Suitability
  MUSuitabilitySol1 <- MUSolver(Suitability, costLayer1, NumSitesGoal1)
  MUSuitabilitySol3 <- MUSolver(Suitability, costLayer3, NumSitesGoal3)
  MUSuitabilitySol4 <- MUSolver(Suitability, costLayer4, NumSitesGoal4)
  maxUtilitySuitSol <- mosaic(MUSuitabilitySol3, MUSuitabilitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MUSuitabilitySol1,  fun="max", na.rm=TRUE)

  # Density  only
  MUDensitySol1 <- MUSolver(Density, costLayer1, NumSitesGoal1)
  MUDensitySol3 <- MUSolver(Density, costLayer3, NumSitesGoal3)
  MUDensitySol4 <- MUSolver(Density, costLayer4, NumSitesGoal4)
  maxUtilityDensitySol <- mosaic(MUDensitySol3, MUDensitySol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MUDensitySol1,  fun="max", na.rm=TRUE)  

# Habitat area only
  MUAreaSol1 <- MUSolver(Area, costLayer1, NumSitesGoal1)
  MUAreaSol3 <- MUSolver(Area, costLayer3, NumSitesGoal3)
  MUAreaSol4 <- MUSolver(Area, costLayer4, NumSitesGoal4)
  maxUtilityAreaSol <- mosaic(MUAreaSol3, MUAreaSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MUAreaSol1,  fun="max", na.rm=TRUE)

#All
  MUAllSol1 <- MUSolver(All, costLayer1, NumSitesGoal1)
  MUAllSol3 <- MUSolver(All, costLayer3, NumSitesGoal3)
  MUAllSol4 <- MUSolver(All, costLayer4, NumSitesGoal4)
  maxUtilityAllSol <- mosaic(MUAllSol3, MUAllSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., MUAllSol1,  fun="max", na.rm=TRUE) 


## 5 Combine and summarize output files --------------------------------------------------------------------

# collect subfolder name
subfolder <- paste0("MapsBudget", Budget)

#"Nice names for output files"
properNames <- c("Generic-Resistance", #1
                                 "Ecoprofile-Taxon", #2
                                 "Ecoprofile-Trophic", #3
                                 "Sum-Species-Resistance", #4
                                 "Mean-Species-Resistance", #5
                                 "Sum-Species-Suitability", #6
                                 "Sum-Species-Density", #7
                                 "Sum-Species-Area", #8
                                 "Sum-Species-All", #9
                                 "Mean-Species-Suitability", #10
                                 "Mean-Species-Density", #11
                                 "Mean-Species-Area", #12
                                 "Mean-Species-All", #13
                                 "Minimize-Shortfall-Species-Suitability", #14
                                 "Minimize-Shortfall-Species-Density", #15
                                 "Minimize-Shortfall-Species-Area", #16
                                 "Minimize-Shortfall-Species-All", #17
                                 "Maximum-Utility-Species-Suitability", #18
                                 "Maximum-Utility-Species-Density", #19
                                 "Maximum-Utility-Species-Area", #20
                                 "Maximum-Utility-Species-All") #21
# Collect outputs - messy names I gave things...
outputAllNames <- c("genericResSol",
                    "FinalEcoprofileTaxon",
                    "FinalEcoprofileTrophic",
                    "FinalSumResDensity",
                    "FinalMeanResDensity",
                    "sum_FinalSuitability",
                    "sum_FinalDensity",
                    "sum_FinalArea",
                    "sum_FinalAll",
                    "mean_FinalSuitability",
                    "mean_FinalDensity",
                    "mean_FinalArea",
                    "mean_FinalAll",
                    "minShortSuitSol",
                    "minShortDensitySol",
                    "minShortAreaSol",
                    "minShortAllSol",
                    "maxUtilitySuitSol",
                    "maxUtilityDensitySol",
                    "maxUtilityAreaSol",
                    "maxUtilityAllSol")

# Match with models selected to get correct names and outputs
outputAll <- eval(parse(text=outputAllNames[[1]]))
for(c in 2:length(outputAllNames)){
  outputAll <- addLayer(outputAll, eval(parse(text=outputAllNames[[c]])))
}

# Export individual maps
writeRaster(outputAll, 
            filename=file.path(outDir, subfolder, paste0(properNames, Budget, ".tif")), 
            bylayer=TRUE,
            overwrite=TRUE)

# Export solution rasters as single rasterStack
writeRaster(outputAll, 
            filename=file.path(outDir, paste0("Allsolutions_", Budget, ".tif")), 
            overwrite=TRUE)


## End script ---------------------------------------------------------------------------