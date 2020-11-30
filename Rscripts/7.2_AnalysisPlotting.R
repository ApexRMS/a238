#####################################################################
# a238 Multispecies prioritization for the Monteregie       
# Plotting and correlations from model prioritization
# 11-2020                                       					
#                             
#	  1.Inputs:                                   
#    -solution raster layers 
#    - 1 raster stack per target budget 
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

 
## Define parameters
mapNames <- c("genericResSol",
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


numModels = length(mapNames)
numSpecies <- length(specieslist)


## Analysis of solution maps ----------------------------------------------------------------


## Load solution maps per budget (for loop starts)

budgets <- c(0.05, 0.1, 0.17)

for(i in budgets){
  
  path <- file.path(outDir, paste0("Allsolutions_", i, ".tif"))
  outputAll <- stack(path)
  names(outputAll) <- properNames
  assign(paste0("outputAll", i), outputAll)

#outputAll0.05
#outputAll0.1
#outputAll0.17


## Calculate correlations ----------------------------------------------------------------------    
 
    set.seed(1)
    # All models
    subSamp <- sampleRandom(outputAll, size= ncell(outputAll[[1]]) * 0.5)
    cors <- cor(subSamp)
    assign(paste0("cors", i), cors)
    
    # Models with all types of inputs
    corsAll <- cor(sampleRandom(subset(outputAll, c(9, 13, 17, 21)), size= ncell(outputAll[[1]]) * 0.5 ))
    assign(paste0("corsAll", i), corsAll)
    
    # Models with density only
    corsDensity <- cor(sampleRandom(subset(outputAll, c(1:5,7,11,15,19)), size= ncell(outputAll[[1]]) * 0.5 ))
    assign(paste0("corsDensity", i), corsDensity)
    
   # Calculate beta diversity (jaccard) - very time consuming
   # jacs <- cross_jaccard(outputAll, thresholds = 0.05)
   # assign(paste0("jacs", i), cors)

} # end budget loop

corrplot(cors0.17, "ellipse", "upper", tl.cex=0.65, diag=FALSE)
corrplot(cors0.17, "number", "upper", number.cex=0.65, tl.cex=0.65, diag=FALSE)

corrplot(corsDensity0.17, "number", "upper", number.cex=0.65, tl.cex=0.65, diag=FALSE)

