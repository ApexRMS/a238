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

# set separate graphics device as the default and cause it to be created/set
.rs.addFunction( "initGraphicsDevice", function()
{
  # options(device="RStudioGD")
  # grDevices::deviceIsInteractive("RStudioGD")
  grDevices::deviceIsInteractive()
})

## Functions - move to source file eventually
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

dev.new()

#setwd("~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity/")
setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

## Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Results"

## Load files ---------------------------------------------------------
  # Focal area
naturalAreasFocal <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(procDataDir, "LULCbinary_FocalArea.tif"))

  # Protected areas
protectedAreas <- raster(file.path(procDataDir, "protectedAreasNatural_FocalArea.tif"))
#values so 0 = protected, 1 = unprotected, reflecting cost of inclusion in network
protectedAreasBinary <-  reclassify(protectedAreas, rcl=matrix(c(0, 1, 1, 0), ncol=2, byrow=T))
#values so 1 = protected, 0 = unprotected, reflecting cost of inclusion in network
protectedAreasNA <-  reclassify(protectedAreasBinary, rcl=matrix(c(0, NA, 1, 1), ncol=2, byrow=T))

  # Ecoregions
ecoregions <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
              calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
              crop(., naturalAreasFocal) %>%
              mask(., naturalAreasFocal) %>%
              calc(., fun=function(x){ifelse(x==1, 3, x)}) #(zone 1 is combined with zone 3)

  # Set up 2 zones 
zone3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
zone4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})

naturalAreasBinaryFocal3 <-  mask(naturalAreasBinaryFocal, zone3)
naturalAreasBinaryFocal4 <-  mask(naturalAreasBinaryFocal, zone4)

protectedAreasNA3 <-  mask(protectedAreasNA, zone3)
protectedAreasNA4 <-  mask(protectedAreasNA, zone4)

#speciesID <- read.csv(
#					file.path(
#					paste0(rawDataDir, "/Focal Species"), 
#					"Species.csv"), 
#					stringsAsFactors = FALSE)

#specieslist <- speciesID$Code

specieslist <- c("PLCI", "RASY", "URAM", "MAAM", "BLBR")

## Using for loop over all species 5 for 3 layers each
  # all layers normalized
  # log scaling for density values (many orders of magnitude var)
  # Range scaling for input into prioritizr
  
for(i in specieslist){

species <- i	

  # Load Focal species data 
density <- raster(file.path(procDataDir, paste0(species, "_curmap_FocalArea.tif"))) %>%
				calc(., fun = log) %>%
				scale(.) %>%
				calc(., fun = rescaleR) %>%
        calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) %>%
        crop(., naturalAreasBinaryFocal) %>%
				mask(., naturalAreasBinaryFocal)
nam2 <- paste0(species, "_density")
assign(nam2, density)		

habitatSuitability <- raster(file.path(procDataDir, paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
				scale(.) %>%
				#calc(., fun = function(x){rescaleR(x, new.min=0.6, new.max=1)}) %>%
				calc(., rescaleR) %>%
        calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) %>%
				crop(., naturalAreasBinaryFocal) %>%
				mask(., naturalAreasBinaryFocal)
nam3 <- paste0(species, "_habitatSuitability")
assign(nam3, habitatSuitability)

habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
  			scale(.) %>%
  			calc(., fun = rescaleR) %>%
        calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) %>%
				crop(., naturalAreasBinaryFocal) %>%
				mask(., naturalAreasBinaryFocal)
nam4 <- paste0(species, "_habitatArea")
assign(nam4, habitatArea)

	} #Ignore warnings


## Combine features ---------------------------------------------------------

  # One stack with four features per species
for(j in specieslist){
	
	species <- j
	nam <- paste0(species, "stack")
	spStack <- stack(evaltext2(species, c("_density", "_habitatSuitability", "_habitatArea")))
	names(spStack) <- c("density", "habitatSuitability", "habitatArea")
	assign(nam, spStack)
}

  # One stack per feature, all species
suitabilityStack <- stack(evaltext(specieslist, "_habitatSuitability"))  
names(suitabilityStack) <- specieslist
areaStack <- stack(evaltext(specieslist, "_habitatArea"))  
names(areaStack) <- specieslist
densityStack <- stack(evaltext(specieslist, "_density"))  
names(densityStack) <- specieslist

  # Combine for all features into complete stack
allSuitabilities <- c(evaltext(specieslist, "_habitatSuitability"), 
					evaltext(specieslist, "_habitatArea"),
					evaltext(specieslist, "_density"))
allSuitabilities <- stack(allSuitabilities)

names(allSuitabilities) <- c(paste0(specieslist, "_habitatSuitability"), 
							paste0(specieslist, "_habitatArea"),
							paste0(specieslist, "_density"))


## Set up separate ecoregions ------------------------------------------------------------------

suit3 <- mask(suitabilityStack, zone3)
suit4 <- mask(suitabilityStack, zone4)

density3 <- mask(densityStack, zone3)
density4 <- mask(densityStack, zone4)

area3 <- mask(areaStack, zone3)
area4 <- mask(areaStack, zone4)

all3 <- mask(allSuitabilities, zone3)
all4 <- mask(allSuitabilities, zone4)


## Prioritization scenarios evaluating ecoregions separately----------------------------------------

## Set parameter values 
  # Budget
Budget <- 0.25  
costLayer3 <- zone3
costLayer4 <- zone4

protectedAreas3 <- mask(protectedAreasBinary, zone3)
protectedAreas4 <- mask(protectedAreasBinary, zone4)

  # New budgets, by ecoregion size
NumSitesGoal3 <- round(Budget*276727, 0)
NumSitesGoal4 <- round(Budget*173560, 0)


## Scenario 1 - take features from multiple species and kinds of data and summarize into a single layer. 
    # Here, we take the sum or median value per pixel across 5 sp x 3 layers 
    # Rescale for prioritizer

# Scenario 1 A - use  single summary layer approach, but include only 1 type of input feature per species

  # Habitat Suitability

sumSuitability3 <- sum(suit3, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal3) %>%
  mask(., protectedAreasNA3, inv=TRUE) %>%
  calc(., rescaleR)

minValues <- sort(values(sumSuitability3), decreasing=TRUE)[NumSitesGoal3]
sumSuitSolutionBudget3 <- calc(sumSuitability3, fun=function(x){ifelse(x < minValues[1], 0, 1)})

sumSuitability4 <- sum(suit4, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal4) %>%
  mask(., protectedAreasNA4, inv=TRUE) %>%
  calc(., rescaleR)

minValues <- sort(values(sumSuitability4), decreasing=TRUE)[NumSitesGoal4]
sumSuitSolutionBudget4 <- calc(sumSuitability4, fun=function(x){ifelse(x < minValues[1], 0, 1)})

  # Final map, sum suit, all ecoregions + protected areas
sumSuitSol <- merge(sumSuitSolutionBudget3, sumSuitSolutionBudget4)
sumSuitSolPA <- merge(sumSuitSol, protectedAreasNA)

  # Calculate overall representation
sumSuitSolRep <- matrix(NA, ncol=1, nrow=dim(suitabilityStack)[3])
rownames(sumSuitSolRep) <- paste0(names(suitabilityStack), "_habitatSuitability")
for(k in 1:dim(suitabilityStack)[3]){
  coverage <- overlay(x=sumSuitSol, y=suitabilityStack[[k]], fun=function(x, y){(x*y)})
  sumSuitSolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(suitabilityStack[[k]], "sum", na.rm=TRUE)
}

  # Density

sumDensity3 <- sum(density3, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal3) %>%
  mask(., protectedAreasNA3, inv=TRUE) %>%
  calc(., rescaleR)

minValues <- sort(values(sumDensity3), decreasing=TRUE)[NumSitesGoal3]
sumDensitySolutionBudget3 <- calc(sumDensity3, fun=function(x){ifelse(x < minValues[1], 0, 1)})

sumDensity4 <- sum(density4, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal4) %>%
  mask(., protectedAreasNA4, inv=TRUE) %>%
  calc(., rescaleR)

minValues <- sort(values(sumDensity4), decreasing=TRUE)[NumSitesGoal4]
sumDensitySolutionBudget4 <- calc(sumDensity4, fun=function(x){ifelse(x < minValues[1], 0, 1)})

  # Final map, sum density, all ecoregions + protected areas
sumDensitySol <- merge(sumDensitySolutionBudget3, sumDensitySolutionBudget4)
sumDensitySolPA <- merge(sumDensitySol, protectedAreasNA)

  # Calculate overall representation
sumDensitySolRep <- matrix(NA, ncol=1, nrow=dim(densityStack)[3])
rownames(sumDensitySolRep) <- paste0(names(densityStack), "_density")
for(k in 1:dim(densityStack)[3]){
    coverage <- overlay(x=sumDensitySol, y=densityStack[[k]], fun=function(x, y){(x*y)})
    sumDensitySolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(densityStack[[k]], "sum", na.rm=TRUE)
  }  
  
  # Habitat area

sumarea3 <- sum(area3, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal3) %>%
  mask(., protectedAreasNA3, inv=TRUE) %>%
  calc(., rescaleR)

minValues <- sort(values(sumarea3), decreasing=TRUE)[NumSitesGoal3]
sumareaSolutionBudget3 <- calc(sumarea3, fun=function(x){ifelse(x < minValues[1], 0, 1)})

sumarea4 <- sum(area4, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal4) %>%
  mask(., protectedAreasNA4, inv=TRUE) %>%
  calc(., rescaleR)

minValues <- sort(values(sumarea4), decreasing=TRUE)[NumSitesGoal4]
sumareaSolutionBudget4 <- calc(sumarea4, fun=function(x){ifelse(x < minValues[1], 0, 1)})

  # Final map, sum area, all ecoregions + protected areas
  sumAreaSol <- merge(sumareaSolutionBudget3, sumareaSolutionBudget4)
  sumAreaSolPA <- merge(sumAreaSol, protectedAreasNA)

  # Calculate overall representation
  sumAreaSolRep <- matrix(NA, ncol=1, nrow=dim(areaStack)[3])
  rownames(sumAreaSolRep) <- paste0(names(areaStack), "_habitatArea")
   rownames(sumAreaSolRep) <- names(areaStack)
  for(k in 1:dim(areaStack)[3]){
    coverage <- overlay(x=sumAreaSol, y=areaStack[[k]], fun=function(x, y){(x*y)})
    sumAreaSolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(areaStack[[k]], "sum", na.rm=TRUE)
  }  
  
  
# Scenario 1 B - use  single summary layer approach incorporating all 3 layers for all 5 species 

   # Habitat
   
sumall3 <- sum(all3, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal3) %>%
  mask(., protectedAreasNA3, inv=TRUE) %>%
  calc(., rescaleR)

minValues <- sort(values(sumall3), decreasing=TRUE)[NumSitesGoal3]
sumallSolutionBudget3 <- calc(sumall3, fun=function(x){ifelse(x < minValues[1], 0, 1)})

sumall4 <- sum(all4, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal4) %>%
  mask(., protectedAreasNA4, inv=TRUE) %>%
  calc(., rescaleR)
minValues <- sort(values(sumall4), decreasing=TRUE)[NumSitesGoal4]
sumallSolutionBudget4 <- calc(sumall4, fun=function(x){ifelse(x < minValues[1], 0, 1)})

  # Final map, sum suit, all ecoregions + protected areas
sumAllSol <- merge(sumallSolutionBudget3, sumallSolutionBudget4)
sumAllSolPA <- merge(sumAllSol, protectedAreasNA)

  # Calculate overall representation
sumAllSolRep <- matrix(NA, ncol=1, nrow=dim(allSuitabilities)[3])
rownames(sumAllSolRep) <- names(allSuitabilities)
for(k in 1:dim(allSuitabilities)[3]){
    coverage <- overlay(x=sumAllSol, y=allSuitabilities[[k]], fun=function(x, y){(x*y)})
    sumAllSolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(allSuitabilities[[k]], "sum", na.rm=TRUE)
  }  
  
  
  
## Scenario 2 - Zones + Multi features for multi species----------------------------------------------
  # Using minimize_shortfall_objective in prioritizer

## Scenario 2A - single type of feature stack per solution 

  # Suitability

  # zone 3
minShortSuitProb3 <- problem(costLayer3, suit3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortSuitProb3) # < 30s
minShortSuitSol3 <- solve(minShortSuitProb3)
repSuit3 <- feature_representation(minShortSuitProb3, minShortSuitSol3)
freq(minShortSuitSol3)[2, "count"]/c(freq(minShortSuitSol3)[1, "count"] + freq(minShortSuitSol3)[2, "count"])

  # zone 4
minShortSuitProb4 <- problem(costLayer4, suit4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortSuitProb4) # < 30s
minShortSuitSol4 <- solve(minShortSuitProb4)
repSuit4 <- feature_representation(minShortSuitProb4, minShortSuitSol4)
freq(minShortSuitSol4)[2, "count"]/c(freq(minShortSuitSol4)[1, "count"] + freq(minShortSuitSol4)[2, "count"])

  # Final map, sum suit, all ecoregions + protected areas
minShortSuitSol <- merge(minShortSuitSol3, minShortSuitSol4)
minShortSuitSolPA <- merge(minShortSuitSol, protectedAreasNA)

  # Calculate overall representation
minShortSuitSolRep <- matrix(NA, ncol=1, nrow=dim(suitabilityStack)[3])
rownames(minShortSuitSolRep) <- paste0(names(suitabilityStack), "_habitatSuitability")
for(k in 1:dim(suitabilityStack)[3]){
    coverage <- overlay(x=minShortSuitSol, y=suitabilityStack[[k]], fun=function(x, y){(x*y)})
    minShortSuitSolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(suitabilityStack[[k]], "sum", na.rm=TRUE)
  }
  
  # Density  

  # zone 3
minShortdensityProb3 <- problem(costLayer3, density3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortdensityProb3) # < 30s
minShortdensitySol3 <- solve(minShortdensityProb3)
repdensity3 <- feature_representation(minShortdensityProb3, minShortdensitySol3)
freq(minShortdensitySol3)[2, "count"]/c(freq(minShortdensitySol3)[1, "count"] + freq(minShortdensitySol3)[2, "count"])

  # zone 4
minShortdensityProb4 <- problem(costLayer4, density4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortdensityProb4) # < 30s
minShortdensitySol4 <- solve(minShortdensityProb4)
repdensity4 <- feature_representation(minShortdensityProb4, minShortdensitySol4)
freq(minShortdensitySol4)[2, "count"]/c(freq(minShortdensitySol4)[1, "count"] + freq(minShortdensitySol4)[2, "count"])

  # Final map, sum suit, all ecoregions + protected areas
minShortDensitySol <- merge(minShortdensitySol3, minShortdensitySol4)
minShortDensitySolPA <- merge(minShortDensitySol, protectedAreasNA)

  # Calculate overall representation
minShortDensitySolRep <- matrix(NA, ncol=1, nrow=dim(densityStack)[3])
rownames(minShortDensitySolRep) <- paste0(names(densityStack), "_density")
for(k in 1:dim(densityStack)[3]){
    coverage <- overlay(x=minShortDensitySol, y=densityStack[[k]], fun=function(x, y){(x*y)})
    minShortDensitySolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(densityStack[[k]], "sum", na.rm=TRUE)
  }  
  
  
  # Area

  # zone 3
minShortareaProb3 <- problem(costLayer3, area3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortareaProb3) # < 30s
minShortareaSol3 <- solve(minShortareaProb3)
reparea3 <- feature_representation(minShortareaProb3, minShortareaSol3)
freq(minShortareaSol3)[2, "count"]/c(freq(minShortareaSol3)[1, "count"] + freq(minShortareaSol3)[2, "count"])

  # zone 4
minShortareaProb4 <- problem(costLayer4, area4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.8) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortareaProb4) # < 30s
minShortareaSol4 <- solve(minShortareaProb4)
reparea4 <- feature_representation(minShortareaProb4, minShortareaSol4)
freq(minShortareaSol4)[2, "count"]/c(freq(minShortareaSol4)[1, "count"] + freq(minShortareaSol4)[2, "count"])

# Final map, sum suit, all ecoregions + protected areas
minShortAreaSol <- merge(minShortareaSol3, minShortareaSol4)
minShortAreaSolPA <- merge(minShortAreaSol, protectedAreasNA)

# Calculate overall representation
minShortAreaSolRep <- matrix(NA, ncol=1, nrow=dim(areaStack)[3])
rownames(minShortAreaSolRep) <- paste0(names(areaStack), "_habitatArea")
for(k in 1:dim(areaStack)[3]){
    coverage <- overlay(x=minShortAreaSol, y=areaStack[[k]], fun=function(x, y){(x*y)})
    minShortAreaSolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(areaStack[[k]], "sum", na.rm=TRUE)
  }  
  
  
## Scenario 2 B - use  all 3 feature layers for all species
 
 # zone 3
minShortProb3 <- problem(costLayer3, all3) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal3) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA3) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(minShortProb3)
minShortSol3 <- solve(minShortProb3)
repAllSuit3 <- feature_representation(minShortProb3, minShortSol3)
freq(minShortSol3)[2, "count"]/c(freq(minShortSol3)[1, "count"] + freq(minShortSol3)[2, "count"])

 # zone 4
minShortProb4 <- problem(costLayer4, all4) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal4) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasNA4) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(minShortProb4)
minShortSol4 <- solve(minShortProb4)
repAllSuit4 <- feature_representation(minShortProb4, minShortSol4)
freq(minShortSol4)[2, "count"]/c(freq(minShortSol4)[1, "count"] + freq(minShortSol4)[2, "count"])

  # Final map, sum suit, all ecoregions + protected areas
minShortAllSol <- merge(minShortSol3, minShortSol4)
minShortAllSolPA <- merge(minShortAllSol, protectedAreasNA)

  # Calculate overall representation
minShortAllSolRep <- matrix(NA, ncol=1, nrow=dim(allSuitabilities)[3])
rownames(minShortAllSolRep) <- names(allSuitabilities)
for(k in 1:dim(allSuitabilities)[3]){
    coverage <- overlay(x=minShortAllSol, y=allSuitabilities[[k]], fun=function(x, y){(x*y)})
    minShortAllSolRep[k, ] <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(allSuitabilities[[k]], "sum", na.rm=TRUE)
  }  



## Combine representation files --------------------------------------------------------------------

  # Confirmed that num cells selected match budget

sumSuitSolRep
sumDensitySolRep
sumAreaSolRep
sumAllSolRep

minShortSuitSolRep
minShortDensitySolRep
minShortAreaSolRep
minShortAllSolRep

  # Combine
allRep <- matrix(NA, ncol=8, nrow=nrow(sumAllSolRep))
rownames(allRep) <- rownames(sumAllSolRep)
colnames(allRep) <- c("SumAll", "SumSuit", "SumDensity", "SumArea", "MSAll", "MSSuit", "MSDensity", "MSArea")
allRep[, 1] <- sumAllSolRep
allRep[1:5, 2] <- sumSuitSolRep
allRep[11:15, 3] <- sumDensitySolRep
allRep[6:10, 4] <- sumAreaSolRep
allRep[, 5] <- minShortAllSolRep
allRep[1:5, 6] <- minShortSuitSolRep
allRep[11:15, 7] <- minShortDensitySolRep
allRep[6:10, 8] <- minShortAreaSolRep



  # Raster correlations

library(corrplot)


outputAll <- stack(sumAllSol,
              sumSuitSol, 
              sumDensitySol, 
              sumAreaSol,  
              minShortAllSol,
              minShortSuitSol,  
              minShortDensitySol,
              minShortAreaSol)
names(outputAll) <- c("SumAll", "SumSuit", "SumDensity", "SumArea", "MSAll", "MSSuit", "MSDensity", "MSArea")

cors <- layerStats(outputAll, "pearson", na.rm=TRUE)

# pearson correlation coefficient`
#            SumAll    SumSuit   SumDensity    SumArea     MSAll    MSSuit
#SumAll     1.0000000 0.20361908 0.44122686 0.74487334 0.6593783 0.1344801
#SumSuit    0.2036191 1.00000000 0.09491711 0.09333594 0.1410388 0.4068946
#SumDensity 0.4412269 0.09491711 1.00000000 0.24968335 0.2953564 0.0918679
#SumArea    0.7448733 0.09333594 0.24968335 1.00000000 0.7485384 0.0539615
#MSAll      0.6593783 0.14103877 0.29535637 0.74853839 1.0000000 0.2758525
#MSSuit     0.1344801 0.40689458 0.09186790 0.05396150 0.2758525 1.0000000
#MSDensity  0.4065938 0.08998305 0.81861420 0.22649380 0.4180313 0.2406869
#MSArea     0.6160193 0.09767758 0.24525240 0.75969434 0.9294792 0.2318155

 #           MSDensity     MSArea
#SumAll     0.40659380 0.61601935
#SumSuit    0.08998305 0.09767758
#SumDensity 0.81861420 0.24525240
#SumArea    0.22649380 0.75969434
#MSAll      0.41803129 0.92947916
#MSSuit     0.24068688 0.23181551
#MSDensity  1.00000000 0.37472860
#MSArea     0.37472860 1.00000000

corrplot(cors$'pearson correlation coefficient', is.corr=FALSE)






save.image("4_prioritizationScenarioEcoregions_working.RData")
load("4_prioritizationScenarioEcoregions_working.RData")