#####################################################################
# a238 Multispecies prioritization for the Monteregie       
# Explore approaches to prioritizing the region for 5 focal species 
# 10-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    
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
library(sf)
library(prioritizr)
library(vegan)
#library("lpsymphony")
library(Rsymphony)
set.seed(10)

setwd("~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity/")

  # Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Results"

  # Functions
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
naturalAreasFocal <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(procDataDir, "LULCbinary_FocalArea.tif"))

  # Protected areas
protectedAreas <- raster(file.path(procDataDir, "protectedAreasNatural_FocalArea.tif"))
protectedAreasBinary <-  reclassify(protectedAreas, rcl=matrix(c(0, 1, 1, 0), ncol=2, byrow=T))
#values so 0 = protected, 1 = unprotected, reflecting cost of inclusion in network

#speciesID <- read.csv(
#					file.path(
#					paste0(rawDataDir, "/Focal Species"), 
#					"Species.csv"), 
#					stringsAsFactors = FALSE)

#specieslist <- speciesID$Code
specieslist <- c("PLCI", "RASY", "URAM", "MAAM", "BLBR")

## Using for loop over all species 5 for 4 layers each
  # Notes - too little variation in most of these individual layers to use normalization, so have tried range scaling only.
  # log scaling for density values (many orders of magnitude var)
  # take inverse of resistance
  
for(i in specieslist){

species <- i	

  # Load Focal species data 
resistance <- raster(file.path(procDataDir, paste0(species, "_Resistance.tif"))) %>%
				calc(., fun = function(x){33 - x}) %>% #take inverse
				scale(.) %>%
				calc(., fun = rescaleR) %>%
				crop(., naturalAreasBinaryFocal) %>%
				mask(., naturalAreasBinaryFocal)
nam1 <- paste0(species, "_ResistanceInv")
assign(nam1, resistance)		

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
	spStack <- stack(evaltext2(species, c("_ResistanceInv", "_density", "_habitatSuitability", "_habitatArea")))
	names(spStack) <- c("ResistanceInv", "density", "habitatSuitability", "habitatArea")
	assign(nam, spStack)
}

  # One stack per feature, all species
resistanceStack <- stack(evaltext(specieslist, "_ResistanceInv"))  
names(resistanceStack) <- specieslist
suitabilityStack <- stack(evaltext(specieslist, "_habitatSuitability"))  
names(suitabilityStack) <- specieslist
areaStack <- stack(evaltext(specieslist, "_habitatArea"))  
names(areaStack) <- specieslist
densityStack <- stack(evaltext(specieslist, "_density"))  
names(densityStack) <- specieslist

## Combine for all features - omit HabitatArea for now
allSuitabilities <- c(evaltext(specieslist, "_habitatSuitability"), 
					evaltext(specieslist, "_ResistanceInv"),
					evaltext(specieslist, "_density"))
allSuitabilities <- stack(allSuitabilities)

names(allSuitabilities) <- c(paste0(specieslist, "_habitatSuitability"), 
							paste0(specieslist, "_ResistanceInv"),
							paste0(specieslist, "_density"))

plot(allSuitabilities)

  # Prelim plotting
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "YlOrRd")
spplot(allSuitabilities, col.regions= pal, cuts = 8, par.settings = list(fontsize = list(text = 6)))


## Run prioritization scenarios ---------------------------------------------------------

Budget <- 0.25 #% of landscape
LSize <- 450287
NumSitesGoal <- round(Budget * LSize, 0)

  # Scenario 1 - consider only the sum/mean/median of feature values
#Sum
sumAll <- sum(allSuitabilities, na.rm=TRUE) %>%
			mask(., naturalAreasBinaryFocal) %>%
			calc(., rescaleR)
			#scale(.)

sumProblem <- problem(naturalAreasBinaryFocal, sumAll) %>% #input is the cost surface + features 
      add_min_set_objective() %>% #minimize cost surface
      add_relative_targets(Budget) %>%
	  add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
	  add_binary_decisions() #inclusion vs no-inclusion	

#presolve_check(p2rMax)
sumSoln <- solve(sumProblem)
plot(sumSoln)
feature_representation(sumProblem, sumSoln)

fr <- as.data.frame(freq(sumSoln))
fr$count[2]/fr$count[1]
#numbers don't match






#Max - not really useful!
#maxAll <- max(allSuitabilities, na.rm=TRUE) 

medianAll <- sum(allSuitabilities, na.rm=TRUE) %>%
			mask(., naturalAreasBinaryFocal) %>%
			scale(.)





