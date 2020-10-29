#####################################################################
# a238 For ECCC Multi-species connectivity analyses
# Create habitat suitability & habitat patch for 14 species
# 
# 09-2020                                       					
#                                                                                                        
#  Inputs (per species)
#    - HabitatSuitability cross walk 
#    - LULC 
#    - Natural areas raster, protected areas raster
#	 - Habitat patch cross walk
#  Outputs: (all cropped to focal area)
#	 - Habitat Suitability, patch ID, patch Area rasters
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################


## Workspace ---------------------------------------------------------

## Packages
library(tidyverse)
library(raster)
library(sp)

## Directories

#setwd("~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity/")
setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
studyDataDir <- "Data/Raw/StudyArea/SHP"

## Input parameters
suitabilityThreshold <- 60

## Load files and inputs ---------------------------------------------------------------

  # Landuse/landsclass map
LULC <- raster(file.path(procDataDir, "LULC_FocalArea.tif"))
naturalAreas <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
 
  # Species characteristics
minPatchSize <- read.csv(file.path(
              paste0(rawDataDir, "/Focal Species"), 
              "Habitat Patch.csv"), 
              stringsAsFactor = FALSE)
crosswalkHabSuit <- read.csv(file.path(paste0(rawDataDir, "/Focal Species"),
                                       "Habitat Suitability.csv"), 
                                       stringsAsFactor = FALSE)
species <- read.csv(file.path(
                    paste0(rawDataDir, "/Focal Species"), 
                    "Species.csv"),
                    stringsAsFactor = FALSE)


## Generate habitat suitability files for all species --------------------------

  # for loop over all species in species list
specieslist <- species$Code

for(i in specieslist){

species <- i		

  # Generate habitat suitability layer using from crosswalk
suitabilityRaster <- LULC %>%
  calc(., fun = function(x){ifelse(x == -9999, NA, x)}) %>%
  reclassify(., rcl = crosswalkHabSuit[which(crosswalkHabSuit$SpeciesID == species), c("Code", "Value")]) 
  
  # Generate habitat patch layer by constraining patches to those with >suitability threshold
habitatRaster <- Which(suitabilityRaster >= suitabilityThreshold) #habitats above 
conversionFromHa <- res(habitatRaster)[1] * res(habitatRaster)[2] * (1/10000) # Convert from hectares to m

# Create habitat patch layer by constraining to species min patch size
patchSizeThreshold <- minPatchSize$MinimumHabitatPatchSize[minPatchSize$SpeciesID == species]

  # Generate patch IDs by combine neighbouring patches of like suitability
habitatClump <- clump(habitatRaster)
habitatClumpID <- data.frame(freq(habitatClump))
  # Identify patches less than minimum
habitatClumpID <- habitatClumpID[habitatClumpID$count < patchSizeThreshold/conversionFromHa, ]

  # Remove clump observations smaller than minimum habitat patch size (ha)
habitatRaster[Which(habitatClump %in% habitatClumpID$value)] <- 0

  # Create raster with ID's for all patches > min threshold & > suitability threshold
habitatRasterCont <- clump(habitatRaster) 

  # Generate raster with patch area as value by combining & IDing neighbouring patches of like suitability
habitatClumpArea <- clump(habitatRaster)
habitatClumpAreaID <- data.frame(freq(habitatClumpArea))
habitatClumpAreaID[nrow(habitatClumpAreaID), "count"] <- NA 
habitatArea <- calc(habitatClumpArea, 
				fun = function(x){
					habitatClumpAreaID[x, 2] })
		

## Crop to Focal Area extent ---------------------------------------------------------
  # Habitat suitability only
suitabilityRasterFocal <-  suitabilityRaster %>%
						   crop(., naturalAreas) %>%
						   mask(., naturalAreas)
  # Binary, suitable habitat by size and threshold
habitatRasterFocal <- habitatRaster %>%
						   crop(., naturalAreas) %>%
						   mask(., naturalAreas)
  # Habitat area value per pixel
habitatAreaFocal <- habitatArea %>%
						   crop(., naturalAreas) %>%
						   mask(., naturalAreas)	
  # Patch ID
habitatRasterContFocal <-  habitatRasterCont %>%
						   crop(., naturalAreas) %>%
						   mask(., naturalAreas)						   	
		
				
## Save outputs ---------------------------------------------------------

  # All LULC area
writeRaster(suitabilityRaster, 
            file.path(procDataDir, 
            paste0(species, "_HabitatSuitability.tif")), 
            overwrite=TRUE)
writeRaster(habitatRaster, 
            file.path(procDataDir, 
            paste0(species, "_HabitatPatch.tif")), 
            overwrite=TRUE)
writeRaster(habitatArea, 
            file.path(procDataDir, 
            paste0(species, "_HabitatArea.tif")), 
            overwrite=TRUE)
writeRaster(habitatRasterCont, 
            file.path(procDataDir, 
            paste0(species, "_HabitatID.tif")), 
            overwrite=TRUE)

  # Focal area
writeRaster(suitabilityRasterFocal, 
            file.path(procDataDir, 
            paste0(species, "_HabitatSuitability_FocalArea.tif")), 
            overwrite=TRUE)
writeRaster(habitatRasterFocal, 
            file.path(procDataDir, 
            paste0(species, "_HabitatPatch_FocalArea.tif")), 
             overwrite=TRUE)
writeRaster(habitatAreaFocal, 
            file.path(procDataDir, 
            paste0(species, "_HabitatArea_FocalArea.tif")), 
            overwrite=TRUE)
writeRaster(habitatRasterContFocal, 
            file.path(procDataDir, 
            paste0(species, "_HabitatID_FocalArea.tif")), 
            overwrite=TRUE)

} # End loop

## End script
