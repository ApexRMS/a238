
library(raster)
library(sf)

primaryStratum <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/PrimaryStratum_NAchanged.tif")
studyExtent <- raster("Data/Processed/LULC_FocalAreaBuffer.tif")
onlyMonteregie <- st_transform(st_read("Data/Processed/regioMonteregie.shp"),
                               crs(studyExtent))

primaryStratumCropped <- 
  extend(mask(crop(primaryStratum, studyExtent), onlyMonteregie), studyExtent)
primaryStratumCropped[is.na(primaryStratumCropped)] <- 0

writeRaster(primaryStratumCropped, 
            "Data/stsim/primary_stratum_FocalAreaBuffer.tif",
            overwrite = TRUE)
