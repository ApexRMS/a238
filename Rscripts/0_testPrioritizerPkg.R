#####################################################################
# a238 Test Prioritizer library  for 5 test species
# 09-2020                                       					
#                                                                                                        
#  Inputs (per species)
#    - HabitatSuitability, HabitatPatch, and curmap 
#    - Natural areas raster
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

# Workspace ---------------------------------------------------------

## Packages
  #To initially install
install.packages("prioritizr", repos = "https://cran.rstudio.com/")
library(prioritizr)

#CITE: Hanson JO, Schuster R, Morrell N, Strimas-Mackey M, Watts ME, Arcese P, Bennett J, Possingham HP #(2020). prioritizr: Systematic Conservation Prioritization in R. R package version 5.0.2. Available #at https://CRAN.R-project.org/package=prioritizr.

  # Load Additional packages
library(tidyverse)
library(raster)
library(sf)

## Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
dataDir <- file.path(projectDir, "Data/Raw/BTSL_extent")
studyAreaDir <- file.path(projectDir, "Data/Raw/StudyArea/SHP")
outDir <- file.path(projectDir, "Data/Processed")

## Input parameters
species <- "BLBR"		#BLBR  MAAM  URAM  RANA  PLCI
polygonBufferWidth <- 20 # In km
suitabilityThreshold <- 60


## Load files and inputs ------

naturalAreas <- raster(file.path(outDir, "LULC_naturalAreas.tif"))
  # Regio
focalArea <- st_read(file.path(studyAreaDir, paste0("regio_l.shp"))) #regios
  # Protected areas
protectedAreas <- raster(file.path(dataDir, paste0("spatialMultiplier_ProtectedAreas.tif")))
  # update -9999 to NA
protectedAreasNA <- protectedAreas  %>%
					calc(protectedAreas, fun = function(x){ifelse(x==-9999, NA, x)})

  # Focal species data 
habitatSuitability <- raster(file.path(dataDir, paste0("HabitatSuitability.", species, ".it1.ts2010.tif")))
habitatPatch <- raster(file.path(dataDir, paste0("HabitatPatch.", species, ".it1.ts2010.tif")))
curmap <- raster(file.path(dataDir, paste0("curmap_BAUNC_Resistance.", species, ".it1.ts2010.tif")))


