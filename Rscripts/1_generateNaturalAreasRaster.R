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
  protectedAreas <- raster(file.path(dataDir, paste0("spatialMultiplier_ProtectedAreas.tif")))

  # Study area (Monteregie)
  focalArea <- st_read(file.path(outDir, paste0("regioMonteregie.shp"))) #regios
  
  
  # Natural areas
LULC <- raster(file.path(dataDir, "InitialStateClass_AgeMean.tif"))
    # There are 9 forest classes: Decid, Mixed, Conif x Young, Med, Old. 
    # Decid Y, Med, Old (511, 512, 513)
    # Mixed Y, Med, Old (521, 522, 523)
    # Conf Y, Med, Old (531, 532, 533)
    # Urban (400, 410)
      # Subset to natural areas (codes 511, 512, 513, 522, 523, 531, 532, 533, 700) 
#labelClass <- data.frame(ID = c(-9999, 100, 400, 410, 511, 512, 513, 700, 521, 522, 523, 531, 532, 533, 800, 810), type = c("NA", "Agriculture", "Urban", "Urban", "DeciduousForest", "DeciduousForest", "DeciduousForest", "Eau", "MixedForest", "MixedForest", "MixedForest", "ConiferForest", "ConiferForest", "ConiferForest", "Wetland", "Wetland"))


## Reproject shapefile to lcc
  # Create spatial polygon from points
polygon <- focalArea %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4269) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
# Reproject polygon to the same crs as lulcRaw
polygonProjected <- st_transform(polygon, crs(LULC))

## Create buffer around polygon for the study area using polygon buffer width
studyArea <- st_buffer(polygonProjected, (polygonBufferWidth*1000))

  # Create bounding box around study area
  # Calculate 20% of the average of the max length and width of the bounding box to be the study area buffer width
  # Connectivity analyses will be run in this buffered study area to minimize edge effects
studyAreaBBox <- st_bbox(studyArea)
studyAreaBufferWidth <- round(((studyAreaBBox$ymax - studyAreaBBox$ymin) + (studyAreaBBox$xmax - studyAreaBBox$xmin))/2 * 0.2, digits=0) #in m
studyAreaBuffer <- st_buffer(studyArea, studyAreaBufferWidth)


## LULC - Crop data to STUDY AREA BUFFER and standardize rasters ---------------------
  # Crop to buffered study area
LULCbuffer <- LULC %>%
  crop(., extent(studyAreaBuffer), snap="out") %>% # Crop SOLRIS to study area buffer extent
  mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  trim(.) %>% # Trim extra white spaces
  calc(., fun = function(x){ifelse(x==-9999, NA, x)}) 
# unique(LULC_buffer)

## Reclassify natural/non-natural areas to NA 
naturalClass <- data.frame(ID = c(-9999, 100, 400, 410, 511, 512, 513, 700, 521, 522, 523, 531, 532, 533, 800, 810), type = c(NA, NA, NA, NA, 511, 512, 513, 700, 521, 522, 523, 531, 532, 533, 800, 810))
LULCnatural <- LULCbuffer %>% subs(x=., y= naturalClass, by="ID")
LULCbinary <- calc(LULCnatural, fun = function(x){ifelse(is.na(x), NA, 1)}) 

## Crop protected areas to studyAreaBuffer 
protectedAreasBuffer <- protectedAreas %>%
  crop(., extent(studyAreaBuffer), snap="out") %>% 
  mask(., mask=studyAreaBuffer) %>% # Clip to buffered study area
  trim(.) %>% # Trim extra white spaces
  calc(., fun = function(x){ifelse(x==-9999, NA, x)}) 


## Crop to study area
LULCnaturalStudyArea <- LULCnatural %>%
  crop(., extent(studyArea), snap="out") 

LULCbinaryStudyArea <- LULCnatural %>%
  crop(., extent(studyArea), snap="out") 

protectedAreasStudyArea <- protectedAreasBuffer %>%
  crop(., extent(studyArea), snap="out") 
  	
## Crop to focal area
LULCnaturalFocalArea <- LULCnatural %>%
  crop(., extent(polygonProjected), snap="out") 

LULCbinaryFocalArea <- LULCnatural %>%
  crop(., extent(polygonProjected), snap="out") 

protectedAreasFocalArea <- protectedAreasBuffer %>%
  crop(., extent(polygonProjected), snap="out") 

				
## Save natural areas raster
  # With buffer
writeRaster(LULCnatural, file.path(outDir, "LULCnaturalAreasBuffer.tif"), overwrite=TRUE)
writeRaster(LULCbinary, file.path(outDir, "LULCbinaryBuffer.tif"), overwrite=TRUE)
writeRaster(protectedAreaBuffer, file.path(outDir, "protectedAreasBuffer.tif"), overwrite=TRUE)

  # 20 km Study Area
writeRaster(LULCnaturalStudyArea, file.path(outDir, "LULCnaturalStudyArea.tif"), overwrite=TRUE)
writeRaster(LULCbinaryStudyArea, file.path(outDir, "LULCbinaryStudyArea.tif"), overwrite=TRUE)
writeRaster(protectedAreasStudyArea, file.path(outDir, "protectedAreasStudyArea.tif"), overwrite=TRUE)
  
  #Focal area  
writeRaster(LULCnaturalFocalArea, file.path(outDir, "LULCnaturalFocalArea.tif"), overwrite=TRUE)
writeRaster(LULCbinaryFocalArea, file.path(outDir, "LULCbinaryFocalArea.tif"), overwrite=TRUE)
writeRaster(protectedAreasFocalArea, file.path(outDir, "protectedAreasFocalArea.tif"), overwrite=TRUE)


