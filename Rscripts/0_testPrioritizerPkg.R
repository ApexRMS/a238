#####################################################################
# a238 Test Prioritizer library  for 5 test species
# 09-2020                                       					
#                                                                                                        
#  Inputs (per species)
#    - HabitatSuitability, HabitatPatch, and curmap 
#    - Natural areas raster, protected areas raster
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

# Workspace ---------------------------------------------------------

## Packages
  #To initially install
install.packages("prioritizr", repos = "https://cran.rstudio.com/")
library(prioritizr)

#CITE: Hanson JO, Schuster R, Morrell N, Strimas-Mackey M, Watts ME, Arcese P, Bennett J, Possingham HP #(2020). prioritizr: Systematic Conservation Prioritization in R. R package version 5.0.2. Available #at https://CRAN.R-project.org/package=prioritizr.

data(sim_pu_polygons)
head(sim_pu_polygons@data)
spplot(sim_pu_polygons, "locked_in", main = "Planning units in protected areas",xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
data(sim_features)

# plot the distribution of suitable habitat for each feature
plot(sim_features, main = paste("Feature", seq_len(nlayers(sim_features))),
     nr = 2)
p1 <- problem(sim_pu_polygons, features = sim_features,
              cost_column = "cost") %>%
      add_min_set_objective() %>%
      add_relative_targets(0.15) %>%
      add_binary_decisions() %>%
      add_default_solver(gap = 0)
      
      
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
  # Protected areas
protectedAreas <- raster(file.path(dataDir, paste0("spatialMultiplier_ProtectedAreas.tif")))
  # update -9999 to NA

  # Focal species data 
habitatSuitability <- raster(file.path(dataDir, paste0("HabitatSuitability.", species, ".it1.ts2010.tif")))
habitatPatch <- raster(file.path(dataDir, paste0("HabitatPatch.", species, ".it1.ts2010.tif")))
curmap <- raster(file.path(dataDir, paste0("curmap_BAUNC_Resistance.", species, ".it1.ts2010.tif")))

