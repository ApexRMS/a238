#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Prepare data layers for input into Zonation prioritization 
# 12-2020                                       					
# 
#	  Inputs (for focal species):                                   
#    - habitat suitability, habitat area, current density layers
#    - protected areas layer
#	   - natural areas layer 
#	   - ecoregions layer
#
#   - Suitability and Area data were generated with a small amount of random variation 
#       - Names include "R" to identify these. 
#       - Generated w rnorm(x, 0, 0.001) - less than expected d/t meaningful variation
#
#   Outputs: Prepared data for scenarios to be used in Rscripts 7.1 & 7.2
#   - Prepared data separate by type, sp, and ecoregion to be using in Zonation scenarios 
#
#
# Script by Bronwyn Rayfield for ApexRMS 									
#############################################################################################

## Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(raster)
library(sp)

## Directories
rawDataDir <- paste0("Data/Raw")
procDataDir <- paste0("Data/Processed")
prioritizrInputDir <- paste0("Data/Processed/PrioritizationInputs/prioritizR")
zonationInputDir <- paste0("Data/Processed/PrioritizationInputs/Zonation")

##Load files ---------------------------------------------------------

# Load species list
speciesID <- read.csv(
  file.path(
    paste0(rawDataDir, "/Focal Species"), 
    "Species.csv"), 
  stringsAsFactors = FALSE)
speciesList <- speciesID$Code

# Natural areas
naturalAreasFocal <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- Which(naturalAreasFocal>0)

# Protected areas
protectedAreas <- raster(file.path(procDataDir, "protectedAreasTerrestrial.tif")) %>%
                  crop(., naturalAreasFocal)
protectedAreas[is.na(protectedAreas)]<-0

# Ecoregions 
ecoregions <- raster(file.path(procDataDir, "ecoregions_FocalArea.tif")) %>%
  crop(., naturalAreasFocal)
ecoregions[ecoregions==-9999]<-NA
ecoregions[ecoregions==0]<-NA

# Make a Monteregie study area
studyArea <- Which(ecoregions>0)
studyArea[studyArea==0] <- NA

# Unprotected natural areas
unprotectedNaturalAreas<- naturalAreasBinaryFocal - protectedAreas
unprotectedNaturalAreas[unprotectedNaturalAreas<0]<-0
unprotectedNaturalAreasFocal <- mask(unprotectedNaturalAreas, studyArea)
writeRaster(unprotectedNaturalAreasFocal, file.path(procDataDir, "unprotectedNaturalAreas_Focal.tif"), overwrite=T)

# Ecoregions - zone1 Adirondacks, zone 3 = StL lowlands, Zone 4 = appalachians
ecoregionList <- c(1, 4, 3)

# Conservation criteria to be included in prioritization analyses
criteriaList <- c("SuitabilityR", "Density", "AreaR")

for(ecoregionID in ecoregionList){ 
  ecoregion <- calc(ecoregions, fun=function(x){ifelse(x==ecoregionID, 1, NA)})
  # Create a new extent object that adds a buffer of 1 pixel around ecoregion extent
  zonationExtent <- extent(xmin(ecoregion)-res(ecoregion)[1],
                           xmax(ecoregion)+res(ecoregion)[1],
                           ymin(ecoregion)-res(ecoregion)[1],
                           ymax(ecoregion)+res(ecoregion)[1])
  
  # Crop protected areas and natural areas rasters to ecoregion
  protectedAreasEcoregion <- mask(protectedAreas, ecoregion)
  naturalAreasEcoregion <- mask(naturalAreasBinaryFocal, ecoregion)
  unprotectedNaturalAreasEcoregion <- naturalAreasEcoregion - protectedAreasEcoregion
  
  # Extent protected areas and natural areas rasters to new extent
  protectedAreasEcoregionExtend <- extend(protectedAreasEcoregion, zonationExtent, value=NA)
  naturalAreasEcoregionExtend <- extend(naturalAreasEcoregion, zonationExtent, value=NA)
  unprotectedNaturalAreasEcoregionExtend <- extend(unprotectedNaturalAreasEcoregion, zonationExtent, value=NA)
  
  # Save extended protected areas and natural areas rasters
  writeRaster(protectedAreasEcoregionExtend, file.path(zonationInputDir, paste0("zonation_ProtectedAreas_ecoregion", ecoregionID, ".tif")), overwrite=T)
  writeRaster(naturalAreasEcoregionExtend, file.path(zonationInputDir, paste0("zonation_NaturalAreas_ecoregion", ecoregionID, ".tif")), overwrite=T)
  writeRaster(unprotectedNaturalAreasEcoregionExtend, file.path(zonationInputDir, paste0("zonation_UnprotectedNaturalAreas_ecoregion", ecoregionID, ".tif")), overwrite=T)
  
  for(species in speciesList){ 
    for(criterion in criteriaList){
      
      # Raster name with species, conservation criterion, and ecoregion
      rasterName <- paste0(species, "_", criterion, "_ecoregion", ecoregionID, ".tif")      

      # Read in raster
      inRaster <- raster(file.path(prioritizrInputDir, rasterName))
      
      # Extend raster to new extent
      extendedRaster <- extend(inRaster, zonationExtent, value=NA)
      
      # Save extended rasters
      writeRaster(extendedRaster, file.path(zonationInputDir, rasterName), overwrite=T)
    }
  }
}


## End script