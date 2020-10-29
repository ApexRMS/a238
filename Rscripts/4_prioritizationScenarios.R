#####################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region for 5 focal species 
# 10-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    -habitat suitability, habitat area, current density layers
#	- protected areas layer
#	- natural areas layer 
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
#library("lpsymphony")

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

## Combine for all features - omit HabitatArea for now
allSuitabilities <- c(evaltext(specieslist, "_habitatSuitability"), 
					evaltext(specieslist, "_habitatArea"),
					evaltext(specieslist, "_density"))
allSuitabilities <- stack(allSuitabilities)

names(allSuitabilities) <- c(paste0(specieslist, "_habitatSuitability"), 
							paste0(specieslist, "_habitatArea"),
							paste0(specieslist, "_density"))


plot(allSuitabilities)

  # Prelim plotting
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "YlOrRd")
spplot(allSuitabilities, col.regions= pal, cuts = 8, par.settings = list(fontsize = list(text = 6)))


## Run prioritization scenarios ---------------------------------------------------------

  # Inputs for solver
Budget <- 0.25 #c(0.15, 0.225, 0.30) #% of landscape
LSize <- 450287
NumSitesGoal <- round(Budget * LSize, 0)
costLayer <- naturalAreasBinaryFocal
  

## Scenario 1 - take features from multiple species and kinds of data and summarize into a single layer. 
    # Here, we take the sum or median value per pixel across 5 sp x 3 layers 
    # Rescale for prioritizer


## Scenario 1 A - use  single summary layer approach, but include only 1 type of input feature per species

  # Habitat Suitability
sumSuitability <- sum(suitabilityStack, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal) %>%
  mask(., protectedAreasNA, inv=TRUE) %>%
  calc(., rescaleR)

plot(sort(values(sumSuitability), decreasing=TRUE))

minValues <- sort(values(sumSuitability), decreasing=TRUE)[NumSitesGoal]
sumSuitSolutionBudgetL <- calc(sumSuitability, fun=function(x){ifelse(x < minValues[1], 0, 1)})

  # Density
sumDensity <- sum(densityStack, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal) %>%
  mask(., protectedAreasNA, inv=TRUE) %>%
  calc(., rescaleR)
plot(sort(values(sumDensity), decreasing=TRUE))

minValues <- sort(values(sumDensity), decreasing=TRUE)[NumSitesGoal]
sumDensitySolutionBudgetL <- calc(sumDensity, fun=function(x){ifelse(x < minValues[1], 0, 1)})
plot(sumDensitySolutionBudgetL)


  # Habitat area
sumArea <- sum(areaStack, na.rm=TRUE) %>%
  mask(., naturalAreasBinaryFocal) %>%
  mask(., protectedAreasNA, inv=TRUE) %>%
  calc(., rescaleR)
plot(sort(values(sumArea), decreasing=TRUE))

minValues <- sort(values(sumArea), decreasing=TRUE)[NumSitesGoal]
sumAreaSolutionBudgetL <- calc(sumArea, fun=function(x){ifelse(x < minValues[1], 0, 1)})
plot(sumAreaSolutionBudgetL)


## Scenario 1 B - use  single summary layer approach incorporating all 3 layers for all 5 species 

  # Sum
sumAll <- sum(allSuitabilities, na.rm=TRUE) %>%
			mask(., naturalAreasBinaryFocal) %>%
      mask(., protectedAreasNA, inv=TRUE) %>%
			calc(., rescaleR)
plot(sumAll)
hist(values(sumAll))

plot(sort(values(sumAll), decreasing=TRUE))
minValues <- sort(values(sumAll), decreasing=TRUE)[NumSitesGoal]

sumAllsolutionBudgetL <- calc(sumAll, fun=function(x){ifelse(x < minValues[1], 0, 1)})
plot(sumAllsolutionBudgetL)

  # Median
medianAll <- calc(allSuitabilities, median, na.rm=TRUE) %>%
			mask(., naturalAreasBinaryFocal) %>%
      mask(., protectedAreasNA, inv=TRUE) %>%
      calc(., rescaleR)

minValues <- sort(values(medianAll), decreasing=TRUE)[NumSitesGoal]
medianAllsolutionBudgetL <- calc(medianAll, fun=function(x){ifelse(x < minValues[1], 0, 1)})


## Scenario 2 - Multi features for multi species -------------------------------------------------------
  # Using minimize_shortfall_objective in prioritizer

## Scenario 2A - single type of feature stack per solution---------------------------------------------

  # Suitability
minShortSuitProb <- problem(costLayer, suitabilityStack) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal[1]) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortSuitProb) # < 30s
minShortSuitSol <- solve(minShortSuitProb)
repSuit <- feature_representation(minShortSuitProb, minShortSuitSol)
freq(minShortSuitSol)[2, "count"]/c(freq(minShortSuitSol)[1, "count"] + freq(minShortSuitSol)[2, "count"])
plot(minShortSuitSol)

  # Density  
minShortDensityProb <- problem(costLayer, densityStack) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal[1]) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion

presolve_check(minShortDensityProb)  # 256 s 
minShortDensitySol <- solve(minShortDensityProb)
plot(minShortDensitySol)
repDens <- feature_representation(minShortDensityProb, minShortDensitySol)
freq(minShortDensitySol)[2, "count"]/c(freq(minShortDensitySol)[1, "count"] + freq(minShortDensitySol)[2, "count"])

  # Area
minShortAreaProb <- problem(costLayer, areaStack) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal[1]) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
  add_binary_decisions() #inclusion vs no-inclusion	

presolve_check(minShortAreaProb) #<30 s
minShortAreaSol <- solve(minShortAreaProb)
plot(minShortAreaSol)
repArea <- feature_representation(minShortAreaProb, minShortAreaSol)
freq(minShortAreaSol)[2, "count"]/c(freq(minShortAreaSol)[1, "count"] + freq(minShortAreaSol)[2, "count"])


## Scenario 2 B - use  all 3 layers for all 5 species ---------------------------------------------------------

minShortProb <- problem(costLayer, allSuitabilities) %>% #input is the cost surface + features 
  add_min_shortfall_objective(NumSitesGoal[1]) %>% #minimize cost surface
  add_relative_targets(0.80) %>%  # set high to avoid being a constraint??
  add_binary_decisions() %>% #inclusion vs no-inclusion	
  add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
  add_default_solver()
presolve_check(minShortProb)
minShortSol <- solve(minShortProb)
plot(minShortSol)
repAllSuit <- feature_representation(minShortProb, minShortSol)
freq(minShortSol)[2, "count"]/c(freq(minShortSol)[1, "count"] + freq(minShortSol)[2, "count"])

#save.image("4_prioritizationScenario_working.RData")


rasterCorrelation()
cor(values(minShortSol), values(sumAllsolutionBudgetM), use="complete.obs")
#0.81
cor(values(minShortSol), values(minShortDensitySol), use="complete.obs")
