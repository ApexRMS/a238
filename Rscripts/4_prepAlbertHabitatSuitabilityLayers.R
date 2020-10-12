#####################################################################
# a238 ECCC Multi-species connectivity analysis            
# Compare original Albert et al. habitat suitabilities to new layers             
#                                                                   
# Inputs (per species):                                                           
#    - The natural areas LULC layer         
#    - Albert et al. habitat suitability layers                 
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
outDir <- file.path(projectDir, "Data/Processed")
outDirAlbert <- file.path(projectDir, "Data/Processed/Albert")
AlbertDir <- file.path(projectDir, "Albert2017/HabitatSuitability")

  # Landuse/landsclass map
naturalAreas <- raster(file.path(outDir, "LULC_FocalArea.tif"))

  # Load species names
species <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "Species.csv"))
species <- species[1:14, 1:3]
specieslist <- species$Code

  # Define original Albert projection (EPSG 42104)
crsAlbert <- " +proj=tmerc +lat_0=0 +lon_0=-73.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs"

## Reproject and crop original raters ----------------------------------------------
  # For loop over all species
  
for(i in specieslist){

species <- i
AlbertHabSuit <- raster(file.path(AlbertDir, paste0(species, "_pixelquality.asc")))


crs(AlbertHabSuit) <- crsAlbert

  # Change resolution
AlbertHabSuitAgg <- aggregate(AlbertHabSuit, fact=3)

  # Reproject
AlbertHabSuitProj <- projectRaster(AlbertHabSuitAgg, crs=crs(naturalAreas))
AlbertHabSuitSamp <- resample(AlbertHabSuitProj, naturalAreas, "bilinear")

AlbertHabSuitCrop <- AlbertHabSuitSamp %>%
  crop(., extent(naturalAreas), snap="out") %>% # Crop to monteregie extent
  mask(., mask= naturalAreas) %>% # Clip to focal area
  trim(.)  %>% 
  calc(., function(x){x*100})# Trim extra white spaces


writeRaster(AlbertHabSuitCrop, 
				file.path(outDirAlbert,
				paste0(species, "_FocalArea.tif")), 
				overwrite=TRUE)

} #end loop

