#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region over 3 ecoregions
# 10-2020                                       					
# 
# Part B: Run prioritization case  with variable targets 
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

defaultW <- getOption("warn")
options(warn = -1)

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
# Ecoregions - zone1 Adirondacks, zone 4 = StL lowlands, Zone 3 = appalachians

ecozones <- c(1, 3, 4)

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


# Load outputs from Script 6.1

genericRes <- stack(file.path(inputDir,  "genericResAll.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
ResDensity <- stack(file.path(inputDir,  "ResDensity.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(ResDensity) <- c("SumResDensity", "MeanResDensity", "MaxResDensity")

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

# Budget
#Budget <- 0.20  ## proportion of the unprotected natural areas
Budgets <- c(0.05, 0.17, 0.20, 0.25, 0.30, 0.5, 0.7, 0.9) #0.1

# Identify all natural areas available for prioritization and omit PA areas
costLayer1 <- naturalAreasBinaryFocal1 %>%
  mask(., protectedAreasNA1, inv=TRUE) #omit protected area cells
costLayer3 <- naturalAreasBinaryFocal3 %>%
  mask(., protectedAreasNA3, inv=TRUE) #omit protected area cells
costLayer4 <- naturalAreasBinaryFocal4 %>%
  mask(., protectedAreasNA4, inv=TRUE) #omit protected area cells


## 3.1) Begin variable target loop  
for(v in 1:length(Budgets)){ # Per budget level

  Budget <- Budgets[v]

# Ecoregion budgets, set by by ecoregion size (number of natural area pixels in zone)
NumSitesGoal1 <- round(Budget * cellStats(costLayer1, sum), 0)
NumSitesGoal3 <- round(Budget * cellStats(costLayer3, sum), 0)
NumSitesGoal4 <- round(Budget * cellStats(costLayer4, sum), 0)


# Solver for single layer objective problems
prob <- function(x, y, z, first=TRUE, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(z) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}


## 4) Prioritization scenarios, evaluating ecoregions separately------------------------------------

## 4.1) Generic Resistance layer--------------------------------------------------------------------------------
#code will be updated once real data available

for(l in ecozones){   
  
  zone <- eval(parse(text=paste0("zone", l))) 
  costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
  NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
  
  # Solver struggles with these so I loosen the gap size and return first identified soln
  genericResSol <- mask(genericRes, costLayerZ) %>%
    prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=7.5) %>% #input is the cost surface + features 
    solve(.)
  assign(paste0("genericResSol", l), genericResSol)
}
rm(l, costLayerZ, NumSitesGoalZ)

# Final map generic resistance
genericResSol <- mosaic(genericResSol3, genericResSol4, fun="max", na.rm=TRUE) %>%
  mosaic(., genericResSol1, fun="max", na.rm=TRUE)


## 4.3) Summarizes Resistance into Current Density ------------------------------------------------------------

# Current density layers calculated as the sum, mean, or max of species resistance layers
# Input = 1 layer, which is a summary of resistances to produce density

for(l in ecozones){   
  
  zone <- eval(parse(text=paste0("zone", l))) 
  costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
  NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
  # Sum of resistances
  fctSumResDensitySol <-  mask(ResDensity[["SumResDensity"]], costLayerZ) %>%
    prob(costLayerZ, ., NumSitesGoalZ, gap=6) %>% #input is the cost surface + features 
    solve(.)
  assign(paste0("SumResDensity", l), fctSumResDensitySol)
} # end loop
rm(l, costLayerZ, NumSitesGoalZ)

# Final maps, sum, mean, max density
FinalSumResDensity <- mosaic(SumResDensity3, SumResDensity4, fun="max", na.rm=TRUE) %>%
  mosaic(., SumResDensity1, fun="max", na.rm=TRUE)

## 4.4) Summarize species densities --------------------------------------------------

# Take current density features from multiple species and feature types and summarize into a single layer. 
# Calculate  a summary value per pixel across all species and feature layers 

statChoice <- c("sum")

for(k in statChoice){ # loop over stat choices
  
  for(l in ecozones){   
    
    zone <- eval(parse(text=paste0("zone", l))) 
    costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
    NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
    
     # Density  only
    fctDensitySol <- mask(Density, costLayerZ) %>%
      stackApply(., nlayers(.), k, na.rm=TRUE) %>%
      prob(costLayerZ, ., NumSitesGoalZ) %>% #input is the cost surface + features 
      solve(.)
    assign(paste(k, "Density", l, sep="_"), fctDensitySol)
 
        # All
    fctAllSol <-  mask(All, costLayerZ) %>%
      stackApply(., nlayers(.), k, na.rm=TRUE) %>%
      prob(costLayerZ, ., NumSitesGoalZ) %>% #input is the cost surface + features 
      solve(.)
    assign(paste(k, "All", l, sep="_"), fctAllSol)
    
    
    rm(zone, costLayerZ, NumSitesGoalZ, fctSuitabilitySol, fctDensitySol, fctAreaSol, fctAllSol)  
  } # End ecozone loop 
  
  # Final solution maps per k stat
 
  DensitySol <- mosaic(eval(parse(text=paste0(k, "_Density_3"))), 
                       eval(parse(text=paste0(k, "_Density_4"))), 
                       fun="max", na.rm=TRUE) %>%
    mosaic(., eval(parse(text=paste0(k, "_Density_1"))), 
           fun="max", na.rm=TRUE)
  assign(paste0(k, "_FinalDensity"), DensitySol)
  
   AllSol <- mosaic(eval(parse(text=paste0(k, "_All_3"))), 
                   eval(parse(text=paste0(k, "_All_4"))), 
                   fun="max", na.rm=TRUE) %>%
    mosaic(., eval(parse(text=paste0(k, "_All_1"))),
           fun="max", na.rm=TRUE)
  assign(paste0(k, "_FinalAll"), AllSol)
  
} #end stat choice loop     
rm(l)


## 4.5) Using minimize_shortfall_objective in prioritizer--------------------------------

# Set up minimum shortfall solver for single layer objective problems
probMS <- function(x, y, z, first=T, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_min_shortfall_objective(z) %>% #minimize cost surface
    add_relative_targets(0.5) %>%
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}

# Using all input layers, apply a multiobjective algorithm - minimize shortfall
for(l in ecozones){   
  
  zone <- eval(parse(text=paste0("zone", l))) 
  costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
  NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
  
   # Density  
  minShortDensitySol <- mask(Density, costLayerZ) %>%
    probMS(costLayerZ, ., NumSitesGoalZ) %>%
    solve(.)
  assign(paste0("minShortDensitySol", l), minShortDensitySol)
  # All
  minShortAllSol <- mask(All, costLayerZ) %>%
    probMS(costLayerZ, ., NumSitesGoalZ) %>%
    solve(.)
  assign(paste0("minShortAllSol", l), minShortAllSol)
}
rm(l)

# Final maps 
minShortDensitySol <- mosaic(minShortDensitySol3, minShortDensitySol4, fun="max", na.rm=TRUE) %>%
  mosaic(., minShortDensitySol1, fun="max", na.rm=TRUE)
minShortAllSol <- mosaic(minShortAllSol3, minShortAllSol4, fun="max", na.rm=TRUE) %>%
  mosaic(., minShortAllSol1, fun="max", na.rm=TRUE)   

## 4.6) Scenario using max_utility_objective in prioritizer--------------------------------

# Define multiobjective solver, max utility

probMU <- function(x, y, z, first=T, gap=0.9){ #z is # is number of sites to protect
  problem(x, y) %>% #input is the cost surface + features 
    add_max_utility_objective(z) %>% #minimize cost surface
    add_binary_decisions() %>% #inclusion vs no-inclusion	
    add_rsymphony_solver(., gap=gap, first_feasible=first) 
}

# Using all input layer, apply multiobjective algorithm - max utility
for(l in ecozones){   
  
  zone <- eval(parse(text=paste0("zone", l))) 
  costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
  NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
  # Suitability
  # Density  
  maxUtilityDensitySol <- mask(Density, costLayerZ) %>%
    probMU(costLayerZ, ., NumSitesGoalZ) %>%
    solve(.)
  assign(paste0("maxUtilityDensitySol", l), maxUtilityDensitySol)
  # Area
   # All
  maxUtilityAllSol <- mask(All, costLayerZ) %>%
    probMU(costLayerZ, ., NumSitesGoalZ, gap=2) %>%
    solve(.)
  assign(paste0("maxUtilityAllSol", l), maxUtilityAllSol)
}


# Final maps 
maxUtilityDensitySol <- mosaic(maxUtilityDensitySol3, maxUtilityDensitySol4, fun="max", na.rm=TRUE) %>%
  mosaic(., maxUtilityDensitySol1, fun="max", na.rm=TRUE)
maxUtilityAllSol <- mosaic(maxUtilityAllSol3, maxUtilityAllSol4, fun="max", na.rm=TRUE) %>%
  mosaic(., maxUtilityAllSol1, fun="max", na.rm=TRUE)   



## 7 Combine and summarize output files --------------------------------------------------------------------

# Raster stack all solns 
outputAll <- stack(genericResSol,
                   FinalSumResDensity,
                   sum_FinalDensity,
                   sum_FinalAll,
                   minShortDensitySol,
                   minShortAllSol,
                   maxUtilityDensitySol,
                   maxUtilityAllSol)

mapNames <- c("genericResSol",
              "FinalSumResDensity",
              "sum_FinalDensity",
              "sum_FinalAll",
              "minShortDensitySol",
              "minShortAllSol",
              "maxUtilityDensitySol",
              "maxUtilityAllSol")
names(outputAll) <- mapNames
outputAll <- mask(outputAll, Suitability[[1]])


# Export solution rasters 
writeRaster(outputAll, 
            filename=file.path(outDir,  paste0("Allsolutions_", Budget, ".tif")), 
            overwrite=TRUE)


rm(genericResSol,
   FinalSumResDensity,
   sum_FinalDensity,
   sum_FinalAll,
   minShortDensitySol,
   minShortAllSol,
   maxUtilityDensitySol,
   maxUtilityAllSol)
} # End budget loop

## End script ---------------------------------------
