#####################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region for 5 focal species per 2 ecoregions
# 10-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    -habitat suitability, habitat area, current density layers
#	- protected areas layer
#	- natural areas layer 
#	- ecoregions layer
#  
#   Outputs:
#    
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
outDir <- "Results"

## Function definitions

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

## Load species list
#speciesID <- read.csv(
#					file.path(
#					paste0(rawDataDir, "/Focal Species"), 
#					"Species.csv"), 
#					stringsAsFactors = FALSE)
#specieslist <- speciesID$Code

specieslist <- c("PLCI", "RASY", "URAM", "MAAM", "BLBR")


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

## Using for loop over all species 5 for 3 layers each
  # all layers normalized
  # log scaling for density values (many orders of magnitude var)
  # Range scaling for input into prioritizr


for(j in c(3, 4)){ # run for both zones, all species

zone <- eval(parse(text=paste0("zone", j))) 
naturalAreasZ <- eval(parse(text=paste0("naturalAreasBinaryFocal", j))) 
protectedAreasZ <- eval(parse(text=paste0("protectedAreasNA", j))) 

  for(i in specieslist){ # run for all species
  
      species <- i
        
density <- raster(file.path(procDataDir, paste0(species, "_curmap_FocalArea.tif"))) %>%
        crop(., naturalAreasZ) %>%
        mask(., naturalAreasZ) %>%
        mask(., protectedAreasZ, inv=TRUE) %>%
        calc(., fun = log) %>% #density has log normal distbn
				scale(.) %>%
				calc(., fun = rescaleR) %>%
        calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam2 <- paste0(species, "_density")
assign(nam2, density)		

habitatSuitability <- raster(file.path(procDataDir, paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      mask(., protectedAreasZ, inv=TRUE) %>%
      calc(., fun = log) %>% #density has log normal distbn
      scale(.) %>%
      calc(., fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nam3 <- paste0(species, "_habitatSuitability")
assign(nam3, habitatSuitability)

habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
  crop(., naturalAreasZ) %>%
  mask(., naturalAreasZ) %>%
  mask(., protectedAreasZ, inv=TRUE) %>%
  calc(., fun = log) %>% #density has log normal distbn
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

} # Finish ecoregion loop

  # Output is All3 and All 4 (rasterStacks with all features), and rasterStacks for Area3/4, Suitability3/4. 
  # and Density3/4

## Build complete rasters covering entire Monteregie 
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


## Set Prioritizr target values ----------------------------------------------------------------

# Budget
Budget <- 0.25  
costLayer3 <- naturalAreasBinaryFocal3
costLayer4 <- naturalAreasBinaryFocal4

# New budgets, set by by ecoregion size (numbe of natural area pixels in zone)
NumSitesGoal3 <- round(Budget * cellStats(naturalAreasBinaryFocal3, sum), 0)
NumSitesGoal4 <- round(Budget * cellStats(naturalAreasBinaryFocal4, sum), 0)


## Prioritization scenarios evaluating ecoregions separately----------------------------------------

## Scenario 1 - take features from multiple species and kinds of data and summarize into a single layer. 
    # Calculate  a summary value per pixel across all species and feature layers 
    # Rescale for prioritizer

 # Choice of potential summary metrics - sum, median, max

statChoice <- c("sum", "max", "mean")

for(k in statChoice){
    
## Habitat Suitability
fctSuitability3 <- stackApply(Suitability3, nlayers(Suitability3), k)
minValues <- sort(values(fctSuitability3), decreasing=TRUE)[NumSitesGoal3]
fctSuitabilitySolutionBudget3 <- fctSuitability3
values(fctSuitabilitySolutionBudget3) <- 0
fctSuitabilitySolutionBudget3[Which(fctSuitability3 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal3 - cellStats(fctSuitabilitySolutionBudget3, sum)
fctSuitabilitySolutionBudget3[sample(Which(fctSuitability3 == minValues, cells=TRUE), overage)] <- 1

fctSuitability4 <- stackApply(Suitability4, nlayers(Suitability4), k)
minValues <- sort(values(fctSuitability4), decreasing=TRUE)[NumSitesGoal4]
fctSuitabilitySolutionBudget4 <- fctSuitability4
values(fctSuitabilitySolutionBudget4) <- 0
fctSuitabilitySolutionBudget4[Which(fctSuitability4 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal4 - cellStats(fctSuitabilitySolutionBudget4, sum)
fctSuitabilitySolutionBudget4[sample(Which(fctSuitability4 == minValues, cells=TRUE), overage)] <- 1

# Final map, sum suit, Suitability ecoregions + protected Suitabilitys
fctSuitabilitySol <- mosaic(fctSuitabilitySolutionBudget3, fctSuitabilitySolutionBudget4, fun="max", na.rm=TRUE) %>%
                     mask(., naturalAreasBinaryFocal)
fctSuitabilitySolPA <- mosaic(fctSuitabilitySol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctSuitabilitySol, fctSuitabilitySolPA)
assign(paste0(k, "_FinalSuitability"), final)
rm(fctSuitability3, fctSuitability4, fctSuitabilitySol, fctSuitabilitySolPA, final)

## Density
fctDensity3 <- stackApply(Density3, nlayers(Density3), k)
minValues <- sort(values(fctDensity3), decreasing=TRUE)[NumSitesGoal3]
fctDensitySolutionBudget3 <- fctDensity3
values(fctDensitySolutionBudget3) <- 0
fctDensitySolutionBudget3[Which(fctDensity3 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal3 - cellStats(fctDensitySolutionBudget3, sum)
fctDensitySolutionBudget3[sample(Which(fctDensity3 == minValues, cells=TRUE), overage)] <- 1

fctDensity4 <- stackApply(Density4, nlayers(Density4), k)
minValues <- sort(values(fctDensity4), decreasing=TRUE)[NumSitesGoal4]
fctDensitySolutionBudget4 <- fctDensity4
values(fctDensitySolutionBudget4) <- 0
fctDensitySolutionBudget4[Which(fctDensity4 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal4 - cellStats(fctDensitySolutionBudget4, sum)
fctDensitySolutionBudget4[sample(Which(fctDensity4 == minValues, cells=TRUE), overage)] <- 1

# Final map, sum suit, Density ecoregions + protected Densitys
fctDensitySol <- mosaic(fctDensitySolutionBudget3, fctDensitySolutionBudget4, fun="max", na.rm=TRUE) %>%
                mask(., naturalAreasBinaryFocal)
fctDensitySolPA <- mosaic(fctDensitySol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctDensitySol, fctDensitySolPA)
assign(paste0(k, "_FinalDensity"), final)
rm(fctDensity3, fctDensity4, fctDensitySol, fctDensitySolPA, final)

## Habitat area
fctArea3 <- stackApply(Area3, nlayers(Area3), k)
minValues <- sort(values(fctArea3), decreasing=TRUE)[NumSitesGoal3]
fctAreaSolutionBudget3 <- fctArea3
values(fctAreaSolutionBudget3) <- 0
fctAreaSolutionBudget3[Which(fctArea3 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal3 - cellStats(fctAreaSolutionBudget3, sum)
fctAreaSolutionBudget3[sample(Which(fctArea3 == minValues, cells=TRUE), overage)] <- 1

fctArea4 <- stackApply(Area4, nlayers(Area4), k)
minValues <- sort(values(fctArea4), decreasing=TRUE)[NumSitesGoal4]
fctAreaSolutionBudget4 <- fctArea4
values(fctAreaSolutionBudget4) <- 0
fctAreaSolutionBudget4[Which(fctArea4 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal4 - cellStats(fctAreaSolutionBudget4, sum)
fctAreaSolutionBudget4[sample(Which(fctArea4 == minValues, cells=TRUE), overage)] <- 1

# Final map, sum suit, Area ecoregions + protected Areas
fctAreaSol <- mosaic(fctAreaSolutionBudget3, fctAreaSolutionBudget4, fun="max", na.rm=TRUE) %>%
              mask(., naturalAreasBinaryFocal)
fctAreaSolPA <- mosaic(fctAreaSol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctAreaSol, fctAreaSolPA)
assign(paste0(k, "_FinalArea"), final)
rm(fctArea3, fctArea4, fctAreaSol, fctAreaSolPA, final)


## All
fctAll3 <- stackApply(All3, nlayers(All3), k)
minValues <- sort(values(fctAll3), decreasing=TRUE)[NumSitesGoal3]
fctAllSolutionBudget3 <- fctAll3
values(fctAllSolutionBudget3) <- 0
fctAllSolutionBudget3[Which(fctAll3 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal3 - cellStats(fctAllSolutionBudget3, sum)
fctAllSolutionBudget3[sample(Which(fctAll3 == minValues, cells=TRUE), overage)] <- 1

fctAll4 <- stackApply(All4, nlayers(All4), k)
minValues <- sort(values(fctAll4), decreasing=TRUE)[NumSitesGoal4]
fctAllSolutionBudget4 <- fctAll4
values(fctAllSolutionBudget4) <- 0
fctAllSolutionBudget4[Which(fctAll4 > minValues, cells=TRUE)] <- 1 
overage <-  NumSitesGoal4 - cellStats(fctAllSolutionBudget4, sum)
fctAllSolutionBudget4[sample(Which(fctAll4 == minValues, cells=TRUE), overage)] <- 1

# Final map, sum suit, all ecoregions + protected Alls
fctAllSol <- mosaic(fctAllSolutionBudget3, fctAllSolutionBudget4, fun="max", na.rm=TRUE) %>%
              mask(., naturalAreasBinaryFocal)
fctAllSolPA <- mosaic(fctAllSol, protectedAreasNA, fun="max", na.rm=TRUE)
final <- stack(fctAllSol, fctAllSolPA)
assign(paste0(k, "_FinalAll"), final)
rm(fctAll3, fctAll4, fctAllSol, fctAllSolPA, final)

} # End statistic choice loop  


## Scenario 2 -  using minimize_shortfall_objective in prioritizer----------------------------------------------

## Suitability
  # zone 3
minShortSuitProb3 <- problem(costLayer3, Suitability3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortSuitProb3) # < 30s
minShortSuitSol3 <- solve(minShortSuitProb3)
repSuit3 <- feature_representation(minShortSuitProb3, minShortSuitSol3)
freq(minShortSuitSol3)[2, "count"]/c(freq(minShortSuitSol3)[1, "count"] + freq(minShortSuitSol3)[2, "count"])

  # zone 4
minShortSuitProb4 <- problem(costLayer4, Suitability4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortSuitProb4) # < 30s
minShortSuitSol4 <- solve(minShortSuitProb4)

  # Final map, sum suit, all ecoregions + protected areas
minShortSuitSol <- merge(minShortSuitSol3, minShortSuitSol4)
minShortSuitSolPA <- merge(minShortSuitSol, protectedAreasNA)
minShort_FinalSuitability <- stack(minShortSuitSol, minShortSuitSolPA)

## Density  
  # zone 3
minShortdensityProb3 <- problem(costLayer3, Density3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortdensityProb3) # < 30s
minShortdensitySol3 <- solve(minShortdensityProb3)

  # zone 4
minShortdensityProb4 <- problem(costLayer4, Density4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortdensityProb4) # < 30s
minShortdensitySol4 <- solve(minShortdensityProb4)

  # Final map, sum suit, all ecoregions + protected areas
minShortDensitySol <- merge(minShortdensitySol3, minShortdensitySol4)
minShortDensitySolPA <- merge(minShortDensitySol, protectedAreasNA)
minShort_FinalDensity <- stack(minShortDensitySol, minShortDensitySolPA)

## Area
  # zone 3
minShortAreaProb3 <- problem(costLayer3, Area3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortAreaProb3) # < 30s
minShortAreaSol3 <- solve(minShortAreaProb3)

  # zone 4
minShortAreaProb4 <- problem(costLayer4, Area4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortAreaProb4) # < 30s
minShortAreaSol4 <- solve(minShortAreaProb4)

# Final map, sum suit, all ecoregions + protected areas
minShortAreaSol <- merge(minShortAreaSol3, minShortAreaSol4)
minShortAreaSolPA <- merge(minShortAreaSol, protectedAreasNA)
minShort_FinalArea <- stack(minShortAreaSol, minShortAreaSolPA)

## All
 # zone 3
minShortProb3 <- problem(costLayer3, All3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(minShortProb3)
minShortSol3 <- solve(minShortProb3)

 # zone 4
minShortProb4 <- problem(costLayer4, All4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(minShortProb4)
minShortSol4 <- solve(minShortProb4)

  # Final map, sum suit, all ecoregions + protected areas
minShortAllSol <- merge(minShortSol3, minShortSol4)
minShortAllSolPA <- merge(minShortAllSol, protectedAreasNA)
minShort_FinalAll <- stack(minShortAllSol, minShortAllSolPA)

##
## Scenario 3 -  using add_max_objective in prioritizer----------------------------------------------

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
maxUtilitySuitSol <- merge(maxUtilitySuitSol3, maxUtilitySuitSol4)
maxUtilitySuitSolPA <- merge(maxUtilitySuitSol, protectedAreasNA)
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
maxUtilityDensitySol <- merge(maxUtilitydensitySol3, maxUtilitydensitySol4)
maxUtilityDensitySolPA <- merge(maxUtilityDensitySol, protectedAreasNA)
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
maxUtilityAreaSol <- merge(maxUtilityAreaSol3, maxUtilityAreaSol4)
maxUtilityAreaSolPA <- merge(maxUtilityAreaSol, protectedAreasNA)
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
maxUtilityAllSol <- merge(maxUtilitySol3, maxUtilitySol4)
maxUtilityAllSolPA <- merge(maxUtilityAllSol, protectedAreasNA)
maxUtility_FinalAll <- stack(maxUtilityAllSol, maxUtilityAllSolPA)

save.image("4_prioritizationScenarioEcoregions_working.RData")


## Combine and summarize output files --------------------------------------------------------------------

numModels <- 20

# Raster stack all solns incl PA
outputAll <- stack(sum_FinalAll,
                   sum_FinalSuitability, 
                   sum_FinalDensity, 
                   sum_FinalArea,
                   max_FinalAll,
                   max_FinalSuitability,
                   max_FinalDensity,
                   max_FinalArea,
                   mean_FinalAll,
                   mean_FinalSuitability,
                   mean_FinalDensity,
                   mean_FinalArea,
                   minShort_FinalAll,
                   minShort_FinalSuitability,
                   minShort_FinalDensity,  
                   minShort_FinalArea,
                   maxUtility_FinalAll,
                   maxUtility_FinalSuitability,
                   maxUtility_FinalDensity,  
                   maxUtility_FinalArea)

mapNames <- paste0(rep(c("SumAll",
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
  "MSAll",
  "MSSuit",
  "MSDensity",
  "MSArea",
  "MUAll",
  "MUSuit",
  "MUDensity",
  "MUArea"), each=2), c("", "_PA"))
names(outputAll) <- mapNames

  # Reorder stack 
outputAll <- subset(outputAll, c(seq(from=1, to=39, by=2), seq(from=2, to=40, by=2)))


## Calculate  representation

  # Matrix in which to save representation values, col=models, rows=features
allModelRep <- matrix(NA, nrow=(length(specieslist) * 3 + 1), ncol=numModels)
colnames(allModelRep) <- names(outputAll)[1:numModels]
rownames(allModelRep) <- c("Size", paste(rep(specieslist, each=3), c("Suit", "Density", "Area"), sep="_"))

# Per model output, extract total amount of each feature included in solution


for(k in 1:numModels){
  modelSol <- outputAll[[k]]
  allModelRep[1, k] <- freq(modelSol)[2, "count"]
  
  for(l in 1:length(specieslist)){
        sp <- specieslist[l]
    #Suitability
  coverage <- overlay(x=modelSol, y=Suitability[[l]], fun=function(x, y){(x * y)})
  suitcover <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(Suitability[[l]], "sum", na.rm=TRUE)
    # Density
  coverage <- overlay(x=modelSol, y=Density[[l]], fun=function(x, y){(x * y)})
  densecover <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(Density[[l]], "sum", na.rm=TRUE)
    #Area
  coverage <- overlay(x=modelSol, y=Area[[l]], fun=function(x, y){(x * y)})
  areacover <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(Area[[l]], "sum", na.rm=TRUE)
  
  allModelRep[((l*3)-1):((l*3)+ 1), k] <- c(suitcover, densecover, areacover)
         }
  }

barplot(allModelRep[1, ], cex.names=0.75, angle=30, beside=TRUE)
barplot(allModelRep[-1, ], cex.names=0.75, angle=30, beside=TRUE)



par(mfrow=c(1,2))
plot(outputAll[["MSAll"]])
plot(outputAll[["MUAll"]])
## Calculate raster correlations

#cors <- layerStats(subset(outputAll, 1:20), "pearson", na.rm=TRUE)
set.seed(10)
subSamp <- sampleRandom(subset(outputAll, 1:20), size= ncell(outputAll[[1]]) * 0.25)
cors <- cor(subSamp)

#                SumAll     SumSuit SumDensity    SumArea      MaxAll     MaxSuit
#SumAll      1.00000000  0.18031173 0.39492364 0.85946678  0.02376495  0.02152637
#SumSuit     0.18031173  1.00000000 0.09935234 0.15752819  0.02224282  0.02319036
#SumDensity  0.39492364  0.09935234 1.00000000 0.27171082  0.02370270  0.02189052
#SumArea     0.85946678  0.15752819 0.27171082 1.00000000  0.02351016  0.02247970
#MaxAll      0.02376495  0.02224282 0.02370270 0.02351016  1.00000000  0.02334434
#MaxSuit     0.02152637  0.02319036 0.02189052 0.02247970  0.02334434  1.00000000
#MaxDensity  0.28527213  0.05557887 0.79551261 0.17142847  0.01991553  0.01886139
#MaxArea     0.65934464  0.02441032 0.24309495 0.72321101  0.01961339  0.01730375
#MeanAll     0.80266011  0.11150790 0.45632494 0.73537352  0.02396630  0.02115919
#MeanSuit    0.17975504  0.87297009 0.09948263 0.15706626  0.02073859  0.02250339
#MeanDensity 0.39492777  0.09935559 0.99999408 0.27171458  0.02369388  0.02189354
#MeanArea    0.56072858 -0.01531539 0.19019824 0.61677831  0.01878429  0.01797888
#MSAll       0.63769316  0.13761795 0.29308980 0.62265307 -0.04505638 -0.04624081
#MSSuit      0.12229810  0.68322472 0.03861469 0.10632536 -0.04626450 -0.04631188
#MSDensity   0.29013255  0.02927832 0.74297243 0.17412202 -0.04505638 -0.04671458
#MSArea      0.60205351  0.08546759 0.20654339 0.63242461 -0.04590917 -0.04652507
#MUAll       0.74470687  0.07564868 0.30769385 0.69155125 -0.04403777 -0.04627634
#MUSuit      0.08577012  0.69011809 0.02981436 0.06042878 -0.04660798 -0.04587364
#MUDensity   0.29052341  0.02109392 0.74470170 0.17415756 -0.04414437 -0.04691593
#MUArea      0.69244960  0.05834418 0.17733531 0.74469654 -0.04502085 -0.04602761
#              MaxDensity      MaxArea    MeanAll    MeanSuit MeanDensity    MeanArea
#SumAll       0.285272129  0.659344637 0.80266011  0.17975504  0.39492777  0.56072858
#SumSuit      0.055578874  0.024410321 0.11150790  0.87297009  0.09935559 -0.01531539
#SumDensity   0.795512607  0.243094949 0.45632494  0.09948263  0.99999408  0.19019824
#SumArea      0.171428466  0.723211014 0.73537352  0.15706626  0.27171458  0.61677831
#MaxAll       0.019915535  0.019613388 0.02396630  0.02073859  0.02369388  0.01878429
#MaxSuit      0.018861388  0.017303753 0.02115919  0.02250339  0.02189354  0.01797888
#MaxDensity   1.000000000  0.174614596 0.38394800  0.05546043  0.79551792  0.14914924
#MaxArea      0.174614596  1.000000000 0.63679299  0.02459983  0.24309863  0.75863726
#MeanAll      0.383947997  0.636792994 1.00000000  0.11163818  0.45632925  0.62872699
#MeanSuit     0.055460431  0.024599830 0.11163818  1.00000000  0.09948588 -0.01507851
#MeanDensity  0.795517923  0.243098630 0.45632925  0.09948588  1.00000000  0.19020176
#MeanArea     0.149149243  0.758637256 0.62872699 -0.01507851  0.19020176  1.00000000
#MSAll        0.189727987  0.544421603 0.53713463  0.13719156  0.29309363  0.43678078
#MSSuit      -0.005644573  0.002249673 0.06212880  0.68259698  0.03861776 -0.04838463
#MSDensity    0.644503666  0.153489288 0.33756917  0.02959812  0.74297760  0.09939646
#MSArea       0.113177950  0.597116806 0.50881479  0.08562156  0.20654696  0.49037615
#MUAll        0.209898916  0.576282717 0.68033283  0.07461823  0.30769772  0.48741508
#MUSuit      -0.008309552 -0.053382912 0.02910675  0.69069846  0.02981741 -0.08565858
#MUDensity    0.664544307  0.150729571 0.34882130  0.02149663  0.74470687  0.10270101
#MUArea       0.090981638  0.628646280 0.61161191  0.05860476  0.17733880  0.53762299
#               MSAll       MSSuit   MSDensity      MSArea       MUAll       MUSuit
#SumAll       0.63769316  0.122298100  0.29013255  0.60205351  0.74470687  0.085770123
#SumSuit      0.13761795  0.683224724  0.02927832  0.08546759  0.07564868  0.690118095
#SumDensity   0.29308980  0.038614689  0.74297243  0.20654339  0.30769385  0.029814363
#SumArea      0.62265307  0.106325364  0.17412202  0.63242461  0.69155125  0.060428778
#MaxAll      -0.04505638 -0.046264499 -0.04505638 -0.04590917 -0.04403777 -0.046607983
#MaxSuit     -0.04624081 -0.046311876 -0.04671458 -0.04652507 -0.04627634 -0.045873638
#MaxDensity   0.18972799 -0.005644573  0.64450367  0.11317795  0.20989892 -0.008309552
#MaxArea      0.54442160  0.002249673  0.15348929  0.59711681  0.57628272 -0.053382912
#MeanAll      0.53713463  0.062128800  0.33756917  0.50881479  0.68033283  0.029106751
#MeanSuit     0.13719156  0.682596977  0.02959812  0.08562156  0.07461823  0.690698465
#MeanDensity  0.29309363  0.038617764  0.74297760  0.20654696  0.30769772  0.029817412
#MeanArea     0.43678078 -0.048384625  0.09939646  0.49037615  0.48741508 -0.085658575
#MSAll        1.00000000  0.336874280  0.46784833  0.88000560  0.84359628  0.305001321
#MSSuit       0.33687428  1.000000000  0.23604392  0.29209106  0.28561223  0.789017835
#MSDensity    0.46784833  0.236043922  1.00000000  0.38981821  0.47150821  0.228842599
#MSArea       0.88000560  0.292091056  0.38981821  1.00000000  0.79362526  0.258204570
#MUAll        0.84359628  0.285612235  0.47150821  0.79362526  1.00000000  0.254627597
#MUSuit       0.30500132  0.789017835  0.22884260  0.25820457  0.25462760  1.000000000
#MUDensity    0.45788729  0.228155631  0.96312876  0.37914652  0.47221887  0.220255496
#MUArea       0.79253559  0.280033579  0.35179807  0.81665054  0.85945577  0.240935609
#MUDensity      MUArea
#SumAll       0.29052341  0.69244960
#SumSuit      0.02109392  0.05834418
#SumDensity   0.74470170  0.17733531
#SumArea      0.17415756  0.74469654
#MaxAll      -0.04414437 -0.04502085
#MaxSuit     -0.04691593 -0.04602761
#MaxDensity   0.66454431  0.09098164
#MaxArea      0.15072957  0.62864628
#MeanAll      0.34882130  0.61161191
#MeanSuit     0.02149663  0.05860476
#MeanDensity  0.74470687  0.17733880
#MeanArea     0.10270101  0.53762299
#MSAll        0.45788729  0.79253559
#MSSuit       0.22815563  0.28003358
#MSDensity    0.96312876  0.35179807
#MSArea       0.37914652  0.81665054
#MUAll        0.47221887  0.85945577
#MUSuit       0.22025550  0.24093561
#MUDensity    1.00000000  0.35199943
#MUArea       0.35199943  1.00000000

corrplot(cors, "ellipse", type="upper", diag=FALSE)

corsAll <- cor(sampleRandom(subset(outputAll, c(1,5,9,13,17)), size= ncell(outputAll[[1]]) * 0.25 ))
corrplot(corsAll, "ellipse", type="upper", diag=FALSE)

# jaccard
jacs <- cross_jaccard(subset(outputAll, 1:20), thresholds = 0.05)
jacsAll <- cross_jaccard(subset(outputAll, c(1,5,9,13,17)), thresholds = 0.05)
corrplot(as.matrix(jacsAll$`0.05`), "ellipse", type="upper", diag=FALSE, is.corr=FALSE)


## Export solution rasters ----------------------------------------------------------------------


writeRaster(outputAll, 
            filename=file.path(outDir, 
                               paste0(names(outputAll), "_prioritizationMap.tif")), 
            bylayer=TRUE,
            overwrite=TRUE)

