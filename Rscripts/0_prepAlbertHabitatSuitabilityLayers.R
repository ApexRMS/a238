#####################################################################
# a238 ECCC Multi-species connectivity analysis            
# Compare original Albert et al. habitat suitabilities to new layers             
#                                                                   
# Inputs (per species):                                                           
#    - The natural areas LULC layer         
#    - Albert et al. habitat suitability layers
#	 - Albert et al patch ID layers
#	 - Albert projection file                 
# Outputs:                                                                                                                                        
#    - Resistance maps                                              
#                                                                   
# Script created by Caroline Tucker  for ApexRMS                    
#####################################################################

## Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(raster)
library(sf)
library(rgdal)

# Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
outDir <- file.path(projectDir, "Data")
outDirAlbert <- file.path(projectDir, "Data/Processed/Albert")
AlbertDir <- file.path(projectDir, "Albert2017/HabitatSuitability")
AlbertDir2 <- file.path(projectDir, "Albert2017/PatchID")

# Landuse/landsclass map
naturalAreas <- raster(file.path(paste0(outDir, "/Processed"),
                                 "LULC_FocalArea.tif"))

# Load species names
species <- read_csv(file.path(
  paste0(outDir, "/Raw/Focal Species"), 
  "Species.csv"))
species <- species[1:14, 1:3]
specieslist <- species$Code

# Input parameters
habitatSuitability <- 60
# Define original Albert projection (EPSG 42104)
crsAlbert <- " +proj=tmerc +lat_0=0 +lon_0=-73.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs"

## Reproject and crop original raters ----------------------------------------------
# For loop over all species

for(i in specieslist){
  
  species <- i
  
  ## Load files
  # Albert files
  AlbertHabSuit <- raster(file.path(AlbertDir, paste0(species, "_pixelquality.asc")))
  AlbertPatchID <- raster(file.path(AlbertDir2, paste0(species, "_PatchId_merged.asc")))
  
  # Set CRS Albert
  crs(AlbertPatchID) <- crs(AlbertHabSuit) <- crsAlbert
  
  # Reproject Albert
  AlbertHabSuitProj <- projectRaster(AlbertHabSuit, crs=crs(naturalAreas))
  
  # Match Updated resolution
  AlbertHabSuitProjAgg <- aggregate(AlbertHabSuitProj, fact=3)
  AlbertHabSuitSamp <- resample(AlbertHabSuitProj, naturalAreas)
  AlbertHabSuitSamp[Which(AlbertHabSuitSamp < (habitatSuitability/100))] <- NA #habitats above 
  
  AlbertPatchProj <- projectRaster(AlbertPatchID, crs=crs(naturalAreas))
  AlbertPatchSamp <- resample(AlbertPatchProj, naturalAreas)
  
  AlbertHabSuitSamp <- mask(AlbertHabSuitSamp, AlbertPatchSamp)
  
  
  # Updated files
  UpdatedHabSuit <- raster(file.path(
    paste0(outDir, "/Processed"), 
    paste0(species, "_HabitatSuitability_FocalArea.tif")))
  UpdatedPatchID <- raster(file.path(
    paste0(outDir, "/Processed"), 
    paste0(species, "_HabitatPatch_FocalArea.tif")))
  
  # Reformat Updated files to match data in Albert
  # Reduce updated habitat suitability layer to include only patches larger than min & habitat suitability > 60
  UpdatedHabSuit[Which(UpdatedHabSuit < habitatSuitability)] <- NA #habitats above 
  UpdatedPatchID <- calc(UpdatedPatchID, fun=function(x){ifelse(x < 1, NA, x)})
  UpdatedHabSuitPatch <- mask(UpdatedHabSuit, UpdatedPatchID)
  
  
  ## Crop and Mask Updated and Albert files to share extent
  AlbertHabSuitCrop <- AlbertHabSuitSamp %>%
    crop(., extent(naturalAreas), snap = "out") %>% # Crop to monteregie extent
    mask(., mask = naturalAreas) %>% # Clip to focal area
    trim(.)  %>% 
    calc(., function(x){ifelse(x < 0, 0, ifelse(x > 1, 1, x)) * 100})
  
  AlbertPatchCrop <- AlbertPatchSamp %>%
    crop(., extent(naturalAreas), snap = "out") %>% # Crop to monteregie extent
    mask(., mask = naturalAreas) %>% # Clip to focal area
    trim(.) # Trim extra white spaces
  
  UpdatedHabSuitCrop <- UpdatedHabSuitPatch %>%
    crop(., AlbertHabSuitCrop, snap = "out") # Crop to Albert extent
  
  ## Export updated raster files ----------------------------------------------
  
  writeRaster(AlbertHabSuitCrop, 
              file.path(outDirAlbert,
                        paste0(species, "_FocalArea.tif")), 
              overwrite=TRUE)
  
  writeRaster(AlbertPatchCrop, 
              file.path(outDirAlbert,
                        paste0(species, "_PatchID.tif")), 
              overwrite=TRUE)
  
  writeRaster(UpdatedHabSuitCrop, 
              file.path(outDirAlbert,
                        paste0(species, "_UpdatedHabitatSuit.tif")), 
              overwrite=TRUE)
  
} #end loop


