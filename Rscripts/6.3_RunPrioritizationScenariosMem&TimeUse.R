#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region over 3 ecoregions
# 10-2020                                       					
# 
# Part C: Calculate time and memory usage per scenario for moderate target. 
#
#	  Inputs:                                
#    - Load input rasters from 6.1_PrepareInputsPrioritizationScenario
#    - protected areas layer
#	   - natural areas layer 
#	   - ecoregions layer
# 
#   Outputs:
#     - solution maps
#   
#   Structure:  1) Workspace
#               2) Load input files
#               3) Set prioritization budgets and costs
#                     3.1) Loop over targets
#               4) Run prioritization scenarios and merge ecoregion-level solutions 
#                     4.1) Generic resistance
#                     4.2) Ecoprofiles
#                     4.3) Combine species Res, then calc Density
#                     4.4) Combine species densities
#                     4.5) Multi objective, Minimize shortfall
#                     4.6) Multi objective, Maximize utility 
#               5) Combine and output solution files
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
library(purrr)
library(pryr)

setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

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

options(warn = -1)

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
#sizeLULC <- 1465929
naturalAreasFocal <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(procDataDir, "LULCbinary_FocalArea.tif"))
#sizenatural <- 450287

# Protected areas
protectedAreas <- raster(file.path(procDataDir, "protectedAreasTerrestrial_FocalArea.tif"))%>%
                  crop(., naturalAreasFocal)

# Ecoregions 
ecoregions <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
  calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
  crop(., naturalAreasFocal) %>%
  mask(., naturalAreasFocal)
# Ecoregions - zone1 Adirondacks, zone 4 = StL lowlands, Zone 3 = appalachians

ecoregionsLULC <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
  calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
  crop(., LULC) %>%
  mask(., LULC)

#
ecozones <- c(1, 3, 4)

# Set up 3 ecoregion zones (1, 3 and 4)
zone1 <- calc(ecoregions, fun=function(x){ifelse(x==1, 1, NA)})
zone3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
zone4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})

naturalAreasBinaryFocal1 <-  mask(naturalAreasBinaryFocal, zone1)
naturalAreasBinaryFocal3 <-  mask(naturalAreasBinaryFocal, zone3)
naturalAreasBinaryFocal4 <-  mask(naturalAreasBinaryFocal, zone4)

protectedAreasNA1 <-  mask(protectedAreas, zone1)
protectedAreasNA3 <-  mask(protectedAreas, zone3)
protectedAreasNA4 <-  mask(protectedAreas, zone4)

  # Calculate zone-specific % protected areas

zone1LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==1, 1, NA)})
zone3LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==3, 1, NA)})
zone4LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==4, 1, NA)})

zonePA1 <- cellStats(protectedAreasNA1, sum, na.rm=T)/cellStats(zone1LULC, sum, na.rm=T)
zonePA3 <- cellStats(protectedAreasNA3, sum, na.rm=T)/cellStats(zone3LULC, sum, na.rm=T)
zonePA4 <- cellStats(protectedAreasNA4, sum, na.rm=T)/cellStats(zone4LULC, sum, na.rm=T)

# Load outputs from Script 6.1

genericRes <- stack(file.path(inputDir,  "genericResAll.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
ResDensity <- stack(file.path(inputDir,  "ResDensity.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(ResDensity) <- c("SumResDensity", "MeanResDensity")

ResTaxon <- stack(file.path(inputDir,  "ResTaxon.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(ResTaxon) <- c("ecoprofileBird", "ecoprofileMammal", "ecoprofileAmphibian")

ResTrophic <- stack(file.path(inputDir,  "ResTrophic.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(ResTrophic) <- c("ecoprofileOmnivore", "ecoprofileInsectivore", "ecoprofileCarnivore", "ecoprofileHerbivore")

Suitability <- stack(file.path(inputDir,  "Suitability.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
Density <- stack(file.path(inputDir,  "Density.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
Area <- stack(file.path(inputDir,  "Area.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(Suitability) <- names(Density) <- names(Area) <- specieslist

All <- stack(file.path(inputDir,  "All.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 



## 3) Set prioritization targets  ----------------------------------------------------------------

# Identify all natural areas available for prioritization and omit PA areas
costLayer1 <- naturalAreasBinaryFocal1 %>%
  mask(., protectedAreasNA1, inv=TRUE) #omit protected area cells
costLayer3 <- naturalAreasBinaryFocal3 %>%
  mask(., protectedAreasNA3, inv=TRUE) #omit protected area cells
costLayer4 <- naturalAreasBinaryFocal4 %>%
  mask(., protectedAreasNA4, inv=TRUE) #omit protected area cells

## Budget


Budget <- 0.17  

# Ecoregion budgets, set by by ecoregion size (number of natural area pixels in zone)
NumSitesGoal1 <- round((Budget-zonePA1) * cellStats(costLayer1, sum), 0) %>%
                ifelse(. <= 0, 1, .)
NumSitesGoal3 <- round((Budget-zonePA3) * cellStats(costLayer3, sum), 0)%>%
  ifelse(. <= 0, 1, .)
NumSitesGoal4 <- round((Budget-zonePA3) * cellStats(costLayer4, sum), 0)%>%
  ifelse(. <= 0, 1, .)
  

## 4) Prioritization scenarios, evaluating ecoregions separately------------------------------------
  
## 4.1) Generic Resistance layer--------------------------------------------------------------------------------


genericTime <- system.time({
  genericResSol1 <- genericSolver(genericRes, costLayer1, NumSitesGoal1)
  genericResSol3 <- genericSolver(genericRes, costLayer3, NumSitesGoal3)
  genericResSol4 <- genericSolver(genericRes, costLayer4, NumSitesGoal4)
  # Final map generic resistance
  genericResSol <- mosaic(genericResSol3, genericResSol4, fun="max", na.rm=TRUE) %>%
    mosaic(., genericResSol1, fun="max", na.rm=TRUE)
  })[3]
)
print(genericTime)
object_size(genericResSol1, genericResSol3, genericResSol4, genericRes, genericResSol)  


## 4.2) Ecoprofile layers--------------------------------------------------------------------------------
  
  # Ecoprofile Scenario 1

genericTime <- system.time({
    
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
  })[3]
print(genericTime)
object_size(resTaxonCombSol1, resTaxonCombSol3, resTaxonCombSol4, ResTaxon, FinalEcoprofileTaxon)  


  # Ecoprofile Scenario 2

genericTime <- system.time({ 
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
  })[3]

print(genericTime)
object_size(resTrophicCombSol1, resTrophicCombSol3, resTrophicCombSol4, ResTrophic, FinalEcoprofileTrophic)  

## 4.3) Summarizes Resistance into Current Density ------------------------------------------------------------
  
  # Current density layers calculated as the sum, mean, or max of species resistance layers
  # Input = 1 layer, which is a summary of resistances to produce density
genericTime <- system.time({
  SumResDensity1 <- genericSolver(ResDensity[["SumResDensity"]], costLayer1, NumSitesGoal1)
  SumResDensity3 <- genericSolver(ResDensity[["SumResDensity"]], costLayer3, NumSitesGoal3)
  SumResDensity4 <- genericSolver(ResDensity[["SumResDensity"]], costLayer4, NumSitesGoal4)
  # Final map generic resistance
  FinalSumResDensity <- mosaic(SumResDensity3, SumResDensity4, fun="max", na.rm=TRUE) %>%
    mosaic(., SumResDensity1, fun="max", na.rm=TRUE)
  })[3]
print(genericTime)
object_size(SumResDensity1, SumResDensity3, SumResDensity4, ResDensity[["SumResDensity"]], FinalSumResDensity)  


genericTime <- system.time({
  MeanResDensity1 <- genericSolver(ResDensity[["MeanResDensity"]], costLayer1, NumSitesGoal1)
  MeanResDensity3 <- genericSolver(ResDensity[["MeanResDensity"]], costLayer3, NumSitesGoal3)
  MeanResDensity4 <- genericSolver(ResDensity[["MeanResDensity"]], costLayer4, NumSitesGoal4)
  # Final map generic resistance
  FinalMeanResDensity <- mosaic(MeanResDensity3, MeanResDensity4, fun="max", na.rm=TRUE) %>%
    mosaic(., MeanResDensity1, fun="max", na.rm=TRUE)
  })[3]

print(genericTime)
object_size(MeanResDensity1, MeanResDensity3, MeanResDensity4, ResDensity[["MeanResDensity"]], FinalMeanResDensity)  
## 4.4) Summarize species densities --------------------------------------------------
  
  # Take current density features from multiple species and feature types and summarize into a single layer. 
  # Calculate  a summary value per pixel across all species and feature layers 
  
statChoice <- c("sum", "mean")
  
for(k in statChoice){ # loop over stat choices
 
 
    # Habitat Suitability only
genericTime <- system.time({    
    SuitabilityK <- stackApply(Suitability, nlayers(Suitability), k, na.rm=TRUE) %>%
                    calc(., fun = rescaleR) 
    SuitabilitySol1 <- genericSolver(SuitabilityK, costLayer1, NumSitesGoal1)
    SuitabilitySol3 <- genericSolver(SuitabilityK, costLayer3, NumSitesGoal3)
    SuitabilitySol4 <- genericSolver(SuitabilityK, costLayer4, NumSitesGoal4)
    SuitabilitySol <- mosaic(SuitabilitySol3, SuitabilitySol4, fun="max", na.rm=TRUE) %>%
      mosaic(., SuitabilitySol1,  fun="max", na.rm=TRUE)
    assign(paste0(k, "_FinalSuitability"), SuitabilitySol)
})[3]
print(genericTime)
object_size(Suitability, SuitabilityK, SuitabilitySol1, SuitabilitySol3, SuitabilitySol4, SuitabilitySol)  
  
    # Density  only
genericTime <- system.time({
    DensityK <- stackApply(Density, nlayers(Density), k, na.rm=TRUE) %>%
                  calc(., fun = rescaleR) 
    DensitySol1 <- genericSolver(DensityK, costLayer1, NumSitesGoal1)
    DensitySol3 <- genericSolver(DensityK, costLayer3, NumSitesGoal3)
    DensitySol4 <- genericSolver(DensityK, costLayer4, NumSitesGoal4)
    DensitySol <- mosaic(DensitySol3, DensitySol4, fun="max", na.rm=TRUE) %>%
      mosaic(., DensitySol1,  fun="max", na.rm=TRUE)
    assign(paste0(k, "_FinalDensity"), DensitySol)
})[3]
object_size(Density, DensityK, DensitySol1, DensitySol3, DensitySol4, DensitySol)  
print(genericTime)
    
    # Habitat area only
genericTime <- system.time({
    AreaK <- stackApply(Area, nlayers(Area), k, na.rm=TRUE) %>%
                calc(., fun = rescaleR) 
    AreaSol1 <- genericSolver(AreaK, costLayer1, NumSitesGoal1)
    AreaSol3 <- genericSolver(AreaK, costLayer3, NumSitesGoal3)
    AreaSol4 <- genericSolver(AreaK, costLayer4, NumSitesGoal4)
    AreaSol <- mosaic(AreaSol3, AreaSol4, fun="max", na.rm=TRUE) %>%
      mosaic(., AreaSol1,  fun="max", na.rm=TRUE)
    assign(paste0(k, "_FinalArea"), AreaSol)
  })[3]
object_size(Area, AreaK, AreaSol1, AreaSol3, AreaSol4, AreaSol) 
print(genericTime)
    
genericTime <- system.time({
    AllK <- stackApply(All, nlayers(All), k, na.rm=TRUE) %>%
               calc(., fun = rescaleR) 
    AllSol1 <- genericSolver(AllK, costLayer1, NumSitesGoal1)
    AllSol3 <- genericSolver(AllK, costLayer3, NumSitesGoal3)
    AllSol4 <- genericSolver(AllK, costLayer4, NumSitesGoal4)
    AllSol <- mosaic(AllSol3, AllSol4, fun="max", na.rm=TRUE) %>%
      mosaic(., AllSol1,  fun="max", na.rm=TRUE)
    assign(paste0(k, "_FinalAll"), AllSol)
  })[3]
print(genericTime)
 object_size(All, AllK, AllSol1, AllSol3, AllSol4, AllSol)  
    
  } #end stat choice loop     
  

## 4.5) Using minimize_shortfall_objective in prioritizer--------------------------------

genericTime <- system.time({
  # Suitability
  MSSuitabilitySol1 <- MSSolver(Suitability, costLayer1, NumSitesGoal1)
  MSSuitabilitySol3 <- MSSolver(Suitability, costLayer3, NumSitesGoal3)
  MSSuitabilitySol4 <- MSSolver(Suitability, costLayer4, NumSitesGoal4)
  minShortSuitSol <- mosaic(MSSuitabilitySol3, MSSuitabilitySol4, fun="max", na.rm=TRUE) %>%
                      mosaic(., MSSuitabilitySol1,  fun="max", na.rm=TRUE)
  })[3]
object_size(MSSuitabilitySol1, MSSuitabilitySol1, MSSuitabilitySol3, MSSuitabilitySol4, minShortSuitSol, Suitability)  
print(genericTime)

# Density  only
genericTime <- system.time({
  MSDensitySol1 <- MSSolver(Density, costLayer1, NumSitesGoal1)
  MSDensitySol3 <- MSSolver(Density, costLayer3, NumSitesGoal3)
  MSDensitySol4 <- MSSolver(Density, costLayer4, NumSitesGoal4)
  minShortDensitySol <- mosaic(MSDensitySol3, MSDensitySol4, fun="max", na.rm=TRUE) %>%
                  mosaic(., MSDensitySol1,  fun="max", na.rm=TRUE)  
    })[3]
print(genericTime)
object_size(MSDensitySol1, MSDensitySol1, MSDensitySol3, MSDensitySol4, minShortDensitySol, Density)  
  
  # Habitat area only
genericTime <- system.time({
  MSAreaSol1 <- MSSolver(Area, costLayer1, NumSitesGoal1)
  MSAreaSol3 <- MSSolver(Area, costLayer3, NumSitesGoal3)
  MSAreaSol4 <- MSSolver(Area, costLayer4, NumSitesGoal4)
  minShortAreaSol <- mosaic(MSAreaSol3, MSAreaSol4, fun="max", na.rm=TRUE) %>%
                mosaic(., MSAreaSol1,  fun="max", na.rm=TRUE)
    })[3]
object_size(MSAreaSol1, MSAreaSol1, MSAreaSol3, MSAreaSol4, minShortAreaSol, Area)  
print(genericTime)
  
  #All
genericTime <- system.time({
  MSAllSol1 <- MSSolver(All, costLayer1, NumSitesGoal1)
  MSAllSol3 <- MSSolver(All, costLayer3, NumSitesGoal3)
  MSAllSol4 <- MSSolver(All, costLayer4, NumSitesGoal4)
  minShortAllSol <- mosaic(MSAllSol3, MSAllSol4, fun="max", na.rm=TRUE) %>%
            mosaic(., MSAllSol1,  fun="max", na.rm=TRUE)    
    })[3]
object_size(MSAllSol1, MSAllSol1, MSAllSol3, MSAllSol4, minShortAllSol, All)  
  print(genericTime)
 
## 4.6) Scenario using max_utility_objective in prioritizer--------------------------------
  
  genericSize <- object_size(
    genericTime <- system.time({
  MUSuitabilitySol1 <- MUSolver(Suitability, costLayer1, NumSitesGoal1)
  MUSuitabilitySol3 <- genericSolver(Suitability, costLayer3, NumSitesGoal3)
  MUSuitabilitySol4 <- genericSolver(Suitability, costLayer4, NumSitesGoal4)
  maxUtilitySuitSol <- mosaic(MUSuitabilitySol3, MUSuitabilitySol4, fun="max", na.rm=TRUE) %>%
                      mosaic(., MUSuitabilitySol1,  fun="max", na.rm=TRUE)
    })[3]
  )
  object_size(MUSuitabilitySol1, MUSuitabilitySol1, MUSuitabilitySol3, MUSuitabilitySol4, maxUtilitySuitSol, Suitability)  
  print(genericTime)
  
  # Density  only
genericTime <- system.time({
  MUDensitySol1 <- MUSolver(Density, costLayer1, NumSitesGoal1)
  MUDensitySol3 <- MUSolver(Density, costLayer3, NumSitesGoal3)
  MUDensitySol4 <- MUSolver(Density, costLayer4, NumSitesGoal4)
  maxUtilityDensitySol <- mosaic(MUDensitySol3, MUDensitySol4, fun="max", na.rm=TRUE) %>%
                 mosaic(., MUDensitySol1,  fun="max", na.rm=TRUE)  
    })[3]
  object_size(MUDensitySol1, MUDensitySol1, MUDensitySol3, MUDensitySol4, maxUtilityDensitySol, Density)  
  print(genericTime)
  
  # Habitat area only
genericTime <- system.time({
  MUAreaSol1 <- MUSolver(Area, costLayer1, NumSitesGoal1)
  MUAreaSol3 <- MUSolver(Area, costLayer3, NumSitesGoal3)
  MUAreaSol4 <- MUSolver(Area, costLayer4, NumSitesGoal4)
  maxUtilityAreaSol <- mosaic(MUAreaSol3, MUAreaSol4, fun="max", na.rm=TRUE) %>%
                mosaic(., MUAreaSol1,  fun="max", na.rm=TRUE)
    })[3]
  object_size(MUAreaSol1, MUAreaSol1, MUAreaSol3, MUAreaSol4, maxUtilityAreaSol, Area)  
  print(genericTime)
  

  #All
genericTime <- system.time({ 
  MUAllSol1 <- MUSolver(All, costLayer1, NumSitesGoal1)
  MUAllSol3 <- MUSolver(All, costLayer3, NumSitesGoal3)
  MUAllSol4 <- MUSolver(All, costLayer4, NumSitesGoal4)
  maxUtilityAllSol <- mosaic(MUAllSol3, MUAllSol4, fun="max", na.rm=TRUE) %>%
                mosaic(., MUAllSol1,  fun="max", na.rm=TRUE) 
    })[3]
object_size(MUAllSol1, MUAllSol1, MUAllSol3, MUAllSol4, maxUtilityAllSol, All)  
  print(genericTime)
 
  
  
## End script ---------------------------------------------------------------------------






