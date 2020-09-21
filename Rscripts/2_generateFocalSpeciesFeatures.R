#####################################################################
# a238 Create features for 5 test species for prioritizer
# 09-2020                                       					
#                                                                                                        
#  Inputs (per species)
#    - HabitatSuitability, HabitatPatch, and curmap 
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
dataDir <- file.path(projectDir, "Data/Raw/BTSL_extent")
studyAreaDir <- file.path(projectDir, "Data/Raw/StudyArea/SHP")
outDir <- file.path(projectDir, "Data/Processed")

## Input parameters
suitabilityThreshold <- 60

  # min habitat area

## Load files and inputs ------

naturalAreasFocal <- raster(file.path(outDir, "LULCnaturalFocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(outDir, "LULCbinaryFocalArea.tif"))
 
  # Protected areas
protectedAreas <- raster(file.path(outDir, "protectedAreasFocalArea.tif"))
#values so 0 = protected, 1 = unprotected, reflecting cost of inclusion in network
  # Mask protected areas to naturalAreasFocal
protectedAreasFocal <-  mask(protectedAreas, naturalAreasFocal)  
  # Min patch size
minPatchSize <- read_csv(file.path(paste0(projectDir, "/Data/Raw/Focal Species"),"FocalSpeciesMinPatchSize.csv"))

## Generate files for all species using for loop
specieslist <- c("BLBR", "MAAM", "URAM", "RANA", "PLCI")

for(i in specieslist){

species <- i		#BLBR  MAAM  URAM  RANA  PLCI

  # Focal species data 
habitatSuitability <- raster(file.path(dataDir, paste0("HabitatSuitability.", species, ".it1.ts2010.tif")))
#habitatPatch <- raster(file.path(dataDir, paste0("HabitatPatch.", species, ".it1.ts2010.tif")))
curmap <- raster(file.path(dataDir, paste0("curmap_BAUNC_Resistance.", species, ".it1.ts2010.tif")))

  # Generate additional patch info
habitatRaster <- Which(habitatSuitability >= suitabilityThreshold) #habitats above suitability threshold
  # Create habitat patch layer by species min patch size ---------------------------------
patchSizeThreshold <- minPatchSize$MinPatchSizeHa[minPatchSize$Species == species]
  #Convert from hectares to m
conversionFromHa <- res(habitatRaster)[1]*res(habitatRaster)[2]*(1/10000)
  #Combine neighbouring patches of like suitability
habitatClump <- clump(habitatRaster)
habitatClumpID <- data.frame(freq(habitatClump))
  # Remove clump observations with frequency smaller than minimum habitat patch size (ha)
habitatClumpID <- habitatClumpID[habitatClumpID$count < patchSizeThreshold/conversionFromHa,]
habitatRaster[Which(habitatClump %in% habitatClumpID$value)] <- 0


  #Combine neighbouring patches of like suitability
habitatClump <- clump(habitatRaster)
habitatClumpID <- data.frame(freq(habitatClump))

habitatClumpID[40913, "count"] <- NA 
habitatArea <- calc(habitatClump, 
				fun = function(x){
					habitatClumpID[x, 2] })
				
					
## Crop focal species features to Monteregie extent

habitatSuitabilityFocal <- habitatSuitability %>%
  crop(., extent(naturalAreasFocal), snap="out") %>% 
  mask(., mask= naturalAreasFocal) 

habitatPatchFocal <- habitatRaster %>%
  crop(., extent(naturalAreasFocal), snap="out") %>% 
  mask(., mask= naturalAreasFocal) 

habitatAreaFocal <- habitatArea %>%
  crop(., extent(naturalAreasFocal), snap="out") %>% 
  mask(., mask= naturalAreasFocal) 

curmapFocal <- curmap %>%
  crop(., extent(naturalAreasFocal), snap="out") %>% 
  mask(., mask= naturalAreasFocal) 


## Save outputs ---------------------------------------------------------

writeRaster(habitatSuitabilityFocal, file.path(outDir, paste0(species, "_HabitatSuitability_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatPatchFocal, file.path(outDir, paste0(species, "_HabitatPatch_FocalArea.tif")), overwrite=TRUE)
writeRaster(habitatAreaFocal, file.path(outDir, paste0(species, "_HabitatArea_FocalArea.tif")), overwrite=TRUE)
writeRaster(curmapFocal, file.path(outDir, paste0(species, "_curmap_FocalArea.tif")), overwrite=TRUE)


}




