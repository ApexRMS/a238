#####################################################################
# a238 ECCC Multi-species connectivity analysis            
# Create resistance maps for 14 species              
#                                                                   
# Inputs (per species):                                                           
#    - The natural areas LULC layer         
#    - The habitat resistance suitability crosswalk                       
# Outputs:                                                                                                                                        
#    - Resistance maps                                              
#                                                                   
# Script created by Caroline Tucker  for ApexRMS                    
#####################################################################

# Workspace ---------------------------------------------------------
# Packages
library(tidyverse)
library(raster)
library(sf)

# Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
dataDir <- file.path(projectDir, "Data/Raw")
studyAreaDir <- file.path(projectDir, "Data/Raw/StudyArea/SHP")
outDir <- file.path(projectDir, "Data/Processed")

# Input parameters
suitabilityThreshold <- 60

## Load data
  # Landuse/landsclass map
LULC <- raster(file.path(outDir, "LULC_FocalArea.tif"))
#naturalAreas <- raster(file.path(outDir, "LULCnatural_FocalArea.tif"))

  # Species characteristics
minPatchSize <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "Habitat Patch.csv"))
crosswalkHabSuit <- read_csv(file.path(paste0(dataDir, "/Focal Species"),"Habitat Suitability.csv"))
crosswalkHabSuit <- crosswalkHabSuit[,1:4]
crosswalkResist <- read_csv(file.path(paste0(dataDir, "/Focal Species"),"Resistance.csv"))
species <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "Species.csv"))
species <- species[1:14, 1:3]


## Generate habitat resistance files for all species ----------------------------------

  # for loop over all species in species list
specieslist <- species$Code

for(i in specieslist){

	species <- i		

  # Reclassify  to generate resistance map
resistanceRasterReclass <- LULC %>%
  calc(., fun = function(x){ifelse(x == -9999, NA, x)}) %>%
  reclassify(., rcl=crosswalkResist[which(crosswalkResist$SpeciesID == species), c("StateClassID", "Value")])

  # Overlay habitat patches where habitat suitability is >60%
patchRaster <- raster(file.path(outDir, paste0(species, "_HabitatPatch.tif")))

  # Reclass to assign habitat patches a resistance value = 1 
resistanceRaster <- overlay(patchRaster, resistanceRasterReclass,  
							fun = function(x, y){return(ifelse(x==1, 1, y))})
resistanceRaster[resistanceRaster < 0] <- NA

## Save outputs ---------------------------------------------------------

writeRaster(resistanceRaster, file.path(outDir, paste0(species, "_Resistance.tif")), overwrite=TRUE)


				} # end loop