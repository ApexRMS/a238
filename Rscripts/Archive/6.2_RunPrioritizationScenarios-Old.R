#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region over 3 ecoregions
# 10-2020                                       					
# 
# Part B: Run prioritization case studies
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
inputDir <- "Data/Processed/InputsPrioritizR"
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
  
ecozones <- c(1, 3, 4)
  
  # Set up 3 ecoregion zones (1, 3 and 4)
zone1 <- calc(ecoregions, fun=function(x){ifelse(x==1, 1, NA)})
zone3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
zone4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})

naturalAreasBinaryFocal1 <-  mask(naturalAreasBinaryFocal, zone1)
naturalAreasBinaryFocal3 <-  mask(naturalAreasBinaryFocal, zone3)
naturalAreasBinaryFocal4 <-  mask(naturalAreasBinaryFocal, zone4)


  # Load outputs from Script 6.1

genericRes <- stack(file.path(inputDir,  "genericResAll.tif"))

ResDensity <- stack(file.path(inputDir,  "ResDensity.tif"))
names(ResDensity) <- c("SumResDensity", "MeanResDensity", "MaxResDensity")

ResTaxon <- stack(file.path(inputDir,  "ResTaxon.tif"))
names(ResTaxon) <- c("ecoprofileBird", "ecoprofileMammal", "ecoprofileAmphibian")

ResTrophic <- stack(file.path(inputDir,  "ResTrophic.tif"))
names(ResTrophic) <- c("ecoprofileOmnivore", "ecoprofileInsectivore", "ecoprofileCarnivore", "ecoprofileHerbivore")

Suitability <- stack(file.path(inputDir,  "Suitability.tif"))
Density <- stack(file.path(inputDir,  "Density.tif"))
Area <- stack(file.path(inputDir,  "Area.tif"))
names(Suitability) <- names(Density) <- names(Area) <- specieslist

All <- stack(file.path(inputDir,  "All.tif"))


## 3) Set prioritization targets  ----------------------------------------------------------------

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


  # Solver for single layer objective problems
prob <- function(x, y, z, first=F, gap=0.9){ #z is # is number of sites to protect
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

## 4.2) Ecoprofile layers--------------------------------------------------------------------------------

for(l in ecozones){   
  
  zone <- eval(parse(text=paste0("zone", l))) 
  costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
  NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
  
  # Scenario 1
fctecoprofileBird <-  mask(ResTaxon[["ecoprofileBird"]], costLayerZ) %>%
                      prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=10) %>% #input is the cost surface + features 
                      solve(.)
assign(paste0("ecoprofileBird", l), fctecoprofileBird)

fctecoprofileMammal <-  mask(ResTaxon[["ecoprofileMammal"]], costLayerZ) %>%
                        prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=10) %>% #input is the cost surface + features 
                        solve(.)
assign(paste0("ecoprofileMammal", l), fctecoprofileMammal)

fctecoprofileAmphibian <-  mask(ResTaxon[["ecoprofileAmphibian"]], costLayerZ) %>%
                            prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=10) %>% #input is the cost surface + features 
                            solve(.)
assign(paste0("ecoprofileAmphibian", l), fctecoprofileAmphibian)
#

fctecoprofileTaxonComb <- mask(ResTaxon, costLayerZ) %>%
                          stackApply(., nlayers(.), "mean", na.rm=TRUE) %>%
                          prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=10) %>% #input is the cost surface + features 
                          solve(.)
assign(paste0("resTaxonCombSol", l), fctecoprofileTaxonComb)

rm(fctecoprofileTaxonComb, fctecoprofileAmphibian, fctecoprofileMammal, fctecoprofileBird)
}
rm(l, costLayerZ, NumSitesGoalZ)

# Final map, ecoprofile taxons 
FinalEcoprofileBird <- mosaic(fctecoprofileBird3, fctecoprofileBird4, fun="max", na.rm=TRUE) %>%
  mosaic(., fctecoprofileBird1, fun="max", na.rm=TRUE)
FinalEcoprofileMammal <- mosaic(fctecoprofileMammal3, fctecoprofileMammal4, fun="max", na.rm=TRUE) %>%
  mosaic(., fctecoprofileMammal1, fun="max", na.rm=TRUE)
FinalEcoprofileAmphibian <- mosaic(fctecoprofileAmphibian3, fctecoprofileAmphibian4, fun="max", na.rm=TRUE) %>%
  mosaic(., fctecoprofileAmphibian1, fun="max", na.rm=TRUE)
FinalEcoprofileTaxon <- mosaic(resTaxonCombSol3, resTaxonCombSol4, fun="max", na.rm=TRUE) %>%
  mosaic(., resTaxonCombSol1, fun="max", na.rm=TRUE)


# Scenario 2
for(l in ecozones){   
  
  zone <- eval(parse(text=paste0("zone", l))) 
  costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
  NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
  
fctecoprofileHerbivore <-  mask(ResTrophic[["ecoprofileHerbivore"]], costLayerZ) %>%
                           prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=20) %>% #input is the cost surface + features 
                           solve(.)
assign(paste0("ecoprofileHerbivore", l), fctecoprofileHerbivore)

fctecoprofileOmnivore <-  mask(ResTrophic[["ecoprofileOmnivore"]], costLayerZ) %>%
                          prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=10) %>% #input is the cost surface + features 
                          solve(.)
assign(paste0("ecoprofileOmnivore", l), fctecoprofileOmnivore)

fctecoprofileInsectivore <-  mask(ResTrophic[["ecoprofileInsectivore"]], costLayerZ) %>%
                            prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=10) %>% #input is the cost surface + features 
                            solve(.)
assign(paste0("ecoprofileInsectivore", l), fctecoprofileInsectivore)

fctecoprofileCarnivore <-  mask(ResTrophic[["ecoprofileCarnivore"]], costLayerZ) %>%
                            prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=10) %>% #input is the cost surface + features 
                            solve(.)
assign(paste0("ecoprofileCarnivore", l), fctecoprofileCarnivore)
#
fctecoprofileTrophicComb <- mask(ResTrophic, costLayerZ) %>%
                            stackApply(., nlayers(.), "mean", na.rm=TRUE) %>%
                            prob(costLayerZ, ., NumSitesGoalZ, first=T, gap=20) %>% #input is the cost surface + features 
                            solve(.)
assign(paste0("resTrophicComb", l), fctecoprofileTrophicComb)

rm(fctecoprofileTrophicComb, fctecoprofileCarnivore, fctecoprofileInsectivore, 
   fctecoprofileOmnivore, fctecoprofileHerbivore) 
}
rm(l, costLayerZ, NumSitesGoalZ)


# Final map, ecoprofile trophic 
FinalEcoprofileHerbivore <- mosaic(fctecoprofileHerbivore3, fctecoprofileHerbivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., fctecoprofileHerbivore1, fun="max", na.rm=TRUE)
FinalEcoprofileInsectivore <- mosaic(fctecoprofileInsectivore3, fctecoprofileInsectivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., fctecoprofileInsectivore1, fun="max", na.rm=TRUE)
FinalEcoprofileOmnivore <- mosaic(fctecoprofileOmnivore3, fctecoprofileOmnivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., fctecoprofileOmnivore1, fun="max", na.rm=TRUE)
FinalEcoprofileCarnivore <- mosaic(fctecoprofileCarnivore3, fctecoprofileCarnivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., fctecoprofileCarnivore1, fun="max", na.rm=TRUE)
FinalEcoprofileTrophic <- mosaic(resTrophicComb3, resTrophicComb4, fun="max", na.rm=TRUE) %>%
  mosaic(., resTrophicComb1, fun="max", na.rm=TRUE)

## 4.3) Summarizes Resistance into Current Density ------------------------------------------------------------

# Current density layers calculated as the sum, mean, or max of species resistance layers
# Input = 1 layer, which is a summary of resistances to produce density

for(l in ecozones){   
  
  zone <- eval(parse(text=paste0("zone", l))) 
  costLayerZ <- eval(parse(text=paste0("costLayer", l))) 
  NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
# Sum of resistances
fctSumResDensitySol <-  mask(ResDensity[["SumResDensity"]], costLayerZ) %>%
                        prob(costLayerZ, ., NumSitesGoalZ, gap=3) %>% #input is the cost surface + features 
                        solve(.)
assign(paste0("SumResDensity", l), fctSumResDensitySol)

# Mean of resistances
fctMeanResDensitySol <-   mask(ResDensity[["MeanResDensity"]], costLayerZ) %>%
                          prob(costLayerZ, ., NumSitesGoalZ, gap=3) %>% #input is the cost surface + features 
                          solve(.)
assign(paste0("MeanResDensity", l), fctMeanResDensitySol)

# Max of resistances
fctMaxResDensitySol <- mask(ResDensity[["MaxResDensity"]], costLayerZ) %>%
                        prob(costLayerZ, ., NumSitesGoalZ, gap=3) %>% #input is the cost surface + features 
                        solve(.)
assign(paste0("MaxResDensity", l), fctMaxResDensitySol)

rm(fctSumResDensitySol, fctMeanResDensitySol, fctMaxResDensitySol) 

} # end loop
rm(l, costLayerZ, NumSitesGoalZ)

  # Final maps, sum, mean, max density
FinalSumResDensity <- mosaic(fctSumResDensitySol3, fctSumResDensitySol4, fun="max", na.rm=TRUE) %>%
                      mosaic(., fctSumResDensitySol1, fun="max", na.rm=TRUE)
FinalMeanResDensity <- mosaic(fctMeanResDensitySol3, fctMeanResDensitySol4, fun="max", na.rm=TRUE) %>%
                      mosaic(., fctMeanResDensitySol1, fun="max", na.rm=TRUE)
FinalMaxResDensity <- mosaic(fctMaxResDensitySol3, fctMaxResDensitySol4, fun="max", na.rm=TRUE) %>%
                      mosaic(., fctMaxResDensitySol1, fun="max", na.rm=TRUE)

     
## 4.4) Summarize species densities --------------------------------------------------

  # Take features from multiple species and feature types and summarize into a single layer. 
  # Calculate  a summary value per pixel across all species and feature layers 

statChoice <- c("sum", "max", "mean")
      
for(k in statChoice){ # loop over stat choices
        
    for(l in ecozones){   
          
zone <- eval(parse(text=paste0("zone", l))) 
costlayerZ <- eval(parse(text=paste0("costLayer", l))) 
NumSitesGoalZ <- eval(parse(text=paste0("NumSitesGoal", l))) 
          
  # Habitat Suitability only
fctSuitabilitySol <- mask(Suitability, costLayerZ) %>%
                      stackApply(., nlayers(.), k, na.rm=TRUE) %>%
                      prob(costlayerZ, ., NumSitesGoalZ) %>% #input is the cost surface + features 
                      solve(.)
assign(paste(k, "Suitability", l, sep="_"), fctSuitabilitySol)

  # Density  only
fctDensitySol <- mask(Density, costLayerZ) %>%
                    stackApply(., nlayers(.), k, na.rm=TRUE) %>%
                    prob(costlayerZ, ., NumSitesGoalZ) %>% #input is the cost surface + features 
                    solve(.)
assign(paste(k, "Density", l, sep="_"), fctDensitySol)

  # Habitat area only
fctAreaSol <- mask(Area, costLayerZ) %>%
              stackApply(., nlayers(.), k, na.rm=TRUE) %>%
              prob(costlayerZ, ., NumSitesGoalZ) %>% #input is the cost surface + features 
              solve(.)
assign(paste(k, "Area", l, sep="_"), fctAreaSol)
          
  # All
fctAllSol <-  mask(All, costLayerZ) %>%
              stackApply(., nlayers(.), k, na.rm=TRUE) %>%
              prob(costlayerZ, ., NumSitesGoalZ) %>% #input is the cost surface + features 
              solve(.)
assign(paste(k, "All", l, sep="_"), fctAllSol)

          
rm(zone, costlayerZ, NumSitesGoalZ, inputs, fctSuitability, fctSuitabilitySol, 
      fctDensitySol, fctAreaSol, fctAllSol)  
      } # End ecozone loop 
        
# Final  solution maps per k stat
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
      
      
## 4.5) Using minimize_shortfall_objective in prioritizer--------------------------------

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
    

## 4.6) Scenario using max_utility_objective in prioritizer--------------------------------

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
    