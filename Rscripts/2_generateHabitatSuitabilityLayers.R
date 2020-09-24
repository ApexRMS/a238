#####################################################################
# a238 Create habitat suitability, patch + resistance layers for 14 species
# 09-2020                                       					
#                                                                                                        
#  Inputs (per species)
#    - HabitatSuitability cross walk 
#    - LULC 
#    - Natural areas raster, protected areas raster
#  Outputs: cropped to focal area 
#	 - Habitat Suitability, patch ID, patch Area
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################


# Workspace ---------------------------------------------------------

## Packages
library(tidyverse)
library(raster)
library(sf)


## Directories

projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
dataDir <- file.path(projectDir, "Data/Raw")
studyAreaDir <- file.path(projectDir, "Data/Raw/StudyArea/SHP")
outDir <- file.path(projectDir, "Data/Processed")

## Load files and inputs ------

  # Landuse/landsclass map
naturalAreas <- raster(file.path(outDir, "LULCnatural.tif"))
naturalAreasBinary <- raster(file.path(outDir, "LULCbinary.tif"))
naturalAreasFocal <- raster(file.path(outDir, "LULCnaturalFocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(outDir, "LULCbinaryFocalArea.tif"))
 
  # Protected areas
protectedAreasNatural <- raster(file.path(outDir, "protectedAreasNatural.tif"))
protectedAreasNaturalFocal <- raster(file.path(outDir, "protectedAreasNaturalFocalArea.tif"))
#values so 0 = protected, 1 = unprotected 

  # Species characteristics
minPatchSize <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "Habitat Patch.csv"))
crosswalkHabSuit <- read_csv(file.path(paste0(dataDir, "/Focal Species"),"Habitat Suitability.csv"))
crosswalkHabSuit <- crosswalkHabSuit[,1:4]
species <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "Species.csv"))
species <- species[1:14, 1:3]

## Input parameters
suitabilityThreshold <- 60


## Generate habitat suitability files for all species------
  # for loop
specieslist <- species$Code

for(i in specieslist){

species <- i		

## Generate habitat suitability from crosswalk
suitabilityRaster <- naturalAreas %>%
  reclassify(., rcl = crosswalkHabSuit[which(crosswalkHabSuit$SpeciesID == species), c("Code", "Value")])

  # Create habitat patch layer by species min patch size ---------------------------------
patchSizeThreshold <- minPatchSize$MinimumHabitatPatchSize[minPatchSize$SpeciesID == species]

## Generate habitat patch raster
habitatRaster <- Which(suitabilityRaster >= suitabilityThreshold) #habitats above suitability threshold
  # Convert from hectares to m
conversionFromHa <- res(habitatRaster)[1] * res(habitatRaster)[2] * (1/10000)
  # Combine neighbouring patches of like suitability
habitatClump <- clump(habitatRaster)
habitatClumpID <- data.frame(freq(habitatClump))
  # Remove clump observations with frequency smaller than minimum habitat patch size (ha)
habitatClumpID <- habitatClumpID[habitatClumpID$count < patchSizeThreshold/conversionFromHa, ]
habitatRaster[Which(habitatClump %in% habitatClumpID$value)] <- 0

  # Create raster with ID's for all patches > min threshold & > suitability threshold
habitatRasterCont <- clump(habitatRaster) 

  #Combine neighbouring patches of like suitability
habitatClumpArea <- clump(habitatRaster)
habitatClumpAreaID <- data.frame(freq(habitatClumpArea))
habitatClumpAreaID[nrow(habitatClumpAreaID), "count"] <- NA 
habitatArea <- calc(habitatClumpArea, 
				fun = function(x){
					habitatClumpAreaID[x, 2] })
				
					
## Crop focal species features to Monteregie extent

habitatSuitabilityFocal <- suitabilityRaster %>%
  crop(., extent(naturalAreasFocal), snap="out") %>% 
  mask(., mask= naturalAreasFocal) 

habitatPatchFocal <- habitatRaster %>%
  crop(., extent(naturalAreasFocal), snap="out") %>% 
  mask(., mask= naturalAreasFocal) 

habitatAreaFocal <- habitatArea %>%
  crop(., extent(naturalAreasFocal), snap="out") %>% 
  mask(., mask= naturalAreasFocal) 


## Save outputs ---------------------------------------------------------

writeRaster(habitatSuitabilityFocal, file.path(outDir, paste0(species, "_HabitatSuitability_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatPatchFocal, file.path(outDir, paste0(species, "_HabitatPatch_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatAreaFocal, file.path(outDir, paste0(species, "_HabitatArea_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatRasterCont, file.path(outDir, paste0(species, "_HabitatID_FocalArea.tif")), overwrite=TRUE)

}


## Circuit maps tbd
#curmap <- raster(file.path(dataDir, paste0("curmap_BAUNC_Resistance.", species, ".it1.ts2010.tif")))

#curmapFocal <- curmap %>%
#  crop(., extent(naturalAreasFocal), snap="out") %>% 
#  mask(., mask= naturalAreasFocal) 



