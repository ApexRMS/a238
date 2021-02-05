#####################################################################
# a238 Multispecies prioritization for the Monteregie       
# Calculate feature representation for different budgets
# 11-2020                                       					
#                             
#	  1.Inputs:                                   
#    - solution raster layers 
#    - 1 raster stack per target budget 
#    - Suitability, Density, and Patch raster stacks
#
#   2. Outputs:
#    - representation of suitability, density, and patch data, per model, per budget (.csv)
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

## Workspace ---------------------------------------------------------

options(warn = -1)

# Packages
library(tidyverse)
library(raster)
library(sp)
library(zonator)
library(corrplot)
library(purrr)

setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

## Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Results/PrioritizationSolutions"

## Load files

# Load species list
speciesID <- read.csv(
  file.path(
    paste0(rawDataDir, "/Focal Species"), 
    "Species.csv"), 
  stringsAsFactors = FALSE)
specieslist <- speciesID$Code

## Load feature files
naturalAreasBinaryFocal <- raster(file.path(procDataDir, "LULCbinary_FocalArea.tif"))

  # Use the feature layers without random variation to calculate actual representation
Suitability <- stack(file.path(procDataDir, "PrioritizationInputs",  "Suitability.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
Density <- stack(file.path(procDataDir,"PrioritizationInputs",   "Density.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
Area <- stack(file.path(procDataDir, "PrioritizationInputs",  "Area.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 
names(Suitability) <- names(Density) <- names(Area) <- specieslist

All <- stack(file.path(procDataDir, "PrioritizationInputs",  "All.tif")) %>%
  crop(., naturalAreasBinaryFocal) %>%
  mask(., naturalAreasBinaryFocal) 

## Define parameters
properNames <- c("Generic-Resistance",
                 "Ecoprofile-Taxon",
                 "Ecoprofile-Trophic",
                 "Sum-Species-Resistance",
                 "Mean-Species-Resistance",
                 "Sum-Species-Suitability",
                 "Sum-Species-Density",
                 "Sum-Species-Area",
                 "Sum-Species-All",
                 "Mean-Species-Suitability",
                 "Mean-Species-Density",
                 "Mean-Species-Area",
                 "Mean-Species-All",
                 "Minimize-Shortfall-Species-Suitability",
                 "Minimize-Shortfall-Species-Density",
                 "Minimize-Shortfall-Species-Area",
                 "Minimize-Shortfall-Species-All",
                 "Maximum-Utility-Species-Suitability",
                 "Maximum-Utility-Species-Density",
                 "Maximum-Utility-Species-Area",
                 "Maximum-Utility-Species-All")


numModels = length(properNames)
numSpecies <- length(specieslist)

## Calculate representation of solution maps ----------------------------------------------------------------


# Load solution maps for each budget of interest (for loop starts)

budgets <- c(0.05, 0.1, 0.17) 
for(i in budgets){
  
  path <- file.path(outDir, paste0("Allsolutions_", i, ".tif"))
  
  outputAll <- stack(path)
  names(outputAll) <- properNames
  
  # Matrix in which to save representation values, col=models, rows=features
allModelRep <- matrix(NA, nrow=(length(specieslist) * 3), ncol=numModels)
colnames(allModelRep) <- properNames
rownames(allModelRep) <- paste(specieslist, rep(c("Suit", "Density", "Area"), each=14), sep="_")

  # Per model output, extract total amount of each feature included in solution
for(m in 1:numModels){
  
  modelSol <- outputAll[[m]]
      #1:nspecies = suitability, (nspecies + 1):2*nspecies = density, (2*nspecies +1):3*nspecies
  
    #Suitability
    Suit_sum <- cellStats(Suitability, "sum", na.rm=TRUE)
    Sover <- overlay(x=modelSol, y=Suitability, fun=function(x, y){(x * y)})
    Sover_sum <- cellStats(Sover, "sum", na.rm=TRUE)
    allModelRep[1:numSpecies, m] <- Sover_sum/Suit_sum
    
    # Density
    Dense_sum <- cellStats(Density, "sum", na.rm=TRUE)
    Dover <- overlay(x=modelSol, y=Density, fun=function(x, y){(x * y)})
    Dover_sum <- cellStats(Dover, "sum", na.rm=TRUE)
    allModelRep[(numSpecies+1):(2*numSpecies), m] <- Dover_sum/Dense_sum

    #Area
    Area_sum <- cellStats(Area, "sum", na.rm=TRUE)
    Aover <- overlay(x=modelSol, y=Area, fun=function(x, y){(x * y)})
    Aover_sum <- cellStats(Aover, "sum", na.rm=TRUE)
    allModelRep[(2*numSpecies+1):(3*numSpecies), m] <- Aover_sum/Area_sum

    rm(Suit_sum, Sover, Sover_sum, Dense_sum, Dover, Dover_sum, Area_sum, Aover, Aover_sum)
    
  } #End model loop

assign(paste0("allModelRep", i), allModelRep)  # for outputs


## Save representation table
write.csv(allModelRep, file.path(outDir, paste0("AllRepresentation_", i, ".csv")))

} # End budget loop

