#####################################################################
# a238 Generate LULC layer for the Monteregie, identify and crop to natural areas
# For ECCC Multi-species connectivity analyses
# 09-2020                                       					
#                                                                                                        
#  Inputs 
#    - Monteregie Study area shpfile
#	 - LULC input raster
#    
#  Outputs
#    - Raster file with non-natural areas coded as NA, cropped to Monteregie focal area
#    - Raster file with protected areas cropped to natural areas in Monteregie focal area
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

# Workspace ---------------------------------------------------------

## Load packages
library(tidyverse)
library(raster)
library(sf)

## Directories ------
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
dataDir <- file.path(projectDir, "Data/Raw/BTSL_extent")
outDir <- file.path(projectDir, "Data/Processed")

## Load files and inputs ------
  # Protected areas
protectedAreas <- raster(file.path(dataDir, "spatialMultiplier_ProtectedAreas.tif"))

  # Study area (Monteregie)
focalArea <- st_read(file.path(outDir, "regioMonteregie.shp")) #political boundary
  
  # Natural areas
LULC <- raster(file.path(dataDir, "InitialStateClass_AgeMean.tif"))
    # There are 9 forest classes: Decid, Mixed, Conif x Young, Med, Old. 
    # Decid Y, Med, Old (511, 512, 513)
    # Mixed Y, Med, Old (521, 522, 523)
    # Conf Y, Med, Old (531, 532, 533)
    # Urban (400)
    # Roads (410)

## Reproject shapefile to lcc
  # Create spatial polygon from points
monteregiePolygon <- focalArea %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4269) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
monteregieProj <- st_transform(monteregiePolygon, crs(LULC)) # Reproject polygon to the same crs as LULC


## Generate natural areas rasters ---------------------------------------------------------
    
# Natural area classes (codes 511, 512, 513, 522, 523, 531, 532, 533) 

  # Reclassify non-natural areas to NA 
naturalClass <- data.frame(ID = c(-9999, 100, 400, 410, 511, 512, 513, 700, 521, 522, 523, 531, 532, 533, 800, 810), type = c(NA, NA, NA, NA, 511, 512, 513, NA, 521, 522, 523, 531, 532, 533, 800, 810))

  # Generate LULC raster including only natural areas
LULCnatural <- subs(x=LULC, y= naturalClass, by="ID") %>%
					calc(., fun = function(x){ifelse(x == -9999, NA, x)}) 
LULCbinary <- calc(LULCnatural, fun = function(x){ifelse(is.na(x), NA, 1)}) 

## Generate protected area raster including only those occurring in natural areas
protectedAreasNatural <- protectedAreas %>%
							mask(., mask = LULCnatural) 
						 	

## Crop areas to Monteregie extent ---------------------------------------------------------
  # Crop LULC
LULCFocalArea <- LULC %>%
  crop(., extent(monteregieProj), snap="out") %>% # Crop to monteregie extent
  mask(., mask= monteregieProj) %>% # Clip to focal area
  trim(.)  # Trim extra white spaces

  # Crop LULC natural
LULCnaturalFocalArea <- LULCnatural %>%
  crop(., extent(monteregieProj), snap="out") %>% # Crop to monteregie extent
  mask(., mask= monteregieProj) %>% # Clip to focal area
  trim(.)  # Trim extra white spaces
  
  # Crop LULC binary
LULCbinaryFocalArea <- LULCbinary %>%
  crop(., extent(monteregieProj), snap="out") %>% # Crop to monteregie extent
  mask(., mask= monteregieProj) %>% # Clip to focal area
  trim(.) # Trim extra white spaces

  # Crop protected areas 
protectedAreasNaturalFocalArea <- protectedAreasNatural %>%
  crop(., extent(monteregieProj), snap="out") %>% 
  mask(., mask= monteregieProj) %>% # Clip to focal areas
  trim(.) %>% # Trim extra white spaces
  calc(., fun = function(x){ifelse(x == -9999, NA, x)}) 

				
## Save natural areas rasters
   # Full extent  
writeRaster(LULCnatural, file.path(outDir, "LULCnatural.tif"), overwrite=TRUE)
writeRaster(LULCbinary, file.path(outDir, "LULCbinary.tif"), overwrite=TRUE)
writeRaster(protectedAreasNatural, file.path(outDir, "protectedAreasNatural.tif"), overwrite=TRUE)

  # Focal area  
writeRaster(LULCFocalArea, file.path(outDir, "LULC_FocalArea.tif"), overwrite=TRUE)
writeRaster(LULCnaturalFocalArea, file.path(outDir, "LULCnatural_FocalArea.tif"), overwrite=TRUE)
writeRaster(LULCbinaryFocalArea, file.path(outDir, "LULCbinary_FocalArea.tif"), overwrite=TRUE)
writeRaster(protectedAreasNaturalFocalArea, file.path(outDir, "protectedAreasNatural_FocalArea.tif"), overwrite=TRUE)



