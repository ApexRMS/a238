#####################################################################
# a238 Generate LULC layer with natural areas only
# 09-2020                                       					
#                                                                                                        
#  Inputs 
#    - Monteregie Study area shpfile
#	 - LULC file for Monteregie
#    
#  Outputs
#    - Raster file with non-natural areas coded as NA, cropped to Monteregie, study area, buffer
#    - Raster file with protected areas, cropped to Monteregie, study area, buffer
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

## Load packages
library(tidyverse)
library(raster)
library(sf)

## Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
dataDir <- file.path(projectDir, "Data/Raw/BTSL_extent")
outDir <- file.path(projectDir, "Data/Processed")

## Load data
  # Protected areas
  protectedAreas <- raster(file.path(dataDir, "spatialMultiplier_ProtectedAreas.tif"))

  # Study area (Monteregie)
  focalArea <- st_read(file.path(outDir, "regioMonteregie.shp")) #regios
  
  # Natural areas
LULC <- raster(file.path(dataDir, "InitialStateClass_AgeMean.tif"))
    # There are 9 forest classes: Decid, Mixed, Conif x Young, Med, Old. 
    # Decid Y, Med, Old (511, 512, 513)
    # Mixed Y, Med, Old (521, 522, 523)
    # Conf Y, Med, Old (531, 532, 533)
    # Urban (400, 410)
      # Subset to natural areas (codes 511, 512, 513, 522, 523, 531, 532, 533, 700) 
#labelClass <- data.frame(ID = c(-9999, 100, 400, 410, 511, 512, 513, 700, 521, 522, 523, 531, 532, 533, 800, 810), type = c("NA", "Agriculture", "Urban", "Urban", "DeciduousForest", "DeciduousForest", "DeciduousForest", "Eau", "MixedForest", "MixedForest", "MixedForest", "ConiferForest", "ConiferForest", "ConiferForest", "Wetland", "Wetland"))

  # Input parameters
polygonBufferWidth <- 20 #km #following C-E Ecopark

## Reproject shapefile to lcc
  # Create spatial polygon from points
polygon <- focalArea %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4269) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
# Reproject polygon to the same crs as lulcRaw
polygonProjected <- st_transform(polygon, crs(LULC))


## Crop data to Monteregie shapefile and standardize rasters ---------------------

## LULC
LULCFocalArea <- LULC %>%
  crop(., extent(polygonProjected), snap="out") %>% # Crop SOLRIS to study area buffer extent
  mask(., mask= polygonProjected) %>% # Clip to buffered study area
  trim(.) %>% # Trim extra white spaces
  calc(., fun = function(x){ifelse(x==-9999, NA, x)}) 
# unique(LULC_buffer)

## Reclassify natural/non-natural areas to NA 
naturalClass <- data.frame(ID = c(-9999, 100, 400, 410, 511, 512, 513, 700, 521, 522, 523, 531, 532, 533, 800, 810), type = c(NA, NA, NA, NA, 511, 512, 513, NA, 521, 522, 523, 531, 532, 533, 800, 810))
LULCnaturalFocalArea <- LULCFocalArea %>% subs(x=., y= naturalClass, by="ID")
LULCbinaryFocalArea <- calc(LULCnaturalFocalArea, fun = function(x){ifelse(is.na(x), NA, 1)}) 

## Crop protected areas to studyAreaBuffer 
protectedAreasFocalArea <- protectedAreas %>%
  crop(., extent(polygonProjected), snap="out") %>% 
  mask(., mask= polygonProjected) %>% # Clip to buffered study area
  trim(.) %>% # Trim extra white spaces
  calc(., fun = function(x){ifelse(x==-9999, NA, x)}) 

				
## Save natural areas raster
  #Focal area  
writeRaster(LULCnaturalFocalArea, file.path(outDir, "LULCnaturalFocalArea.tif"), overwrite=TRUE)
writeRaster(LULCbinaryFocalArea, file.path(outDir, "LULCbinaryFocalArea.tif"), overwrite=TRUE)
writeRaster(protectedAreasFocalArea, file.path(outDir, "protectedAreasFocalArea.tif"), overwrite=TRUE)


