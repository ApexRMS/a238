
library(raster)
library(sf)

primaryStratum <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/PrimaryStratum_NAchanged.tif")
secondaryStratum <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/SecondaryStratum_NAchanged.tif")
tertiaryStratum <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/TertiaryStratum_NAchanged.tif")
stateclasses <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/InitialStateClass_AgeMean_NAchanged.tif")


studyExtent <- raster("Data/Processed/LULC_FocalAreaBuffer.tif")
onlyMonteregie <- st_transform(st_read("Data/Processed/regioMonteregie.shp"),
                               crs(studyExtent))

primaryStratumCropped <- 
  extend(mask(crop(primaryStratum, studyExtent), onlyMonteregie), studyExtent)
primaryStratumCropped[is.na(primaryStratumCropped)] <- 0

secondaryStratumCropped<- 
  extend(mask(crop(secondaryStratum, studyExtent), onlyMonteregie), studyExtent)
tertiaryStratumCropped<- 
  extend(mask(crop(tertiaryStratum, studyExtent), onlyMonteregie), studyExtent)
stateclassesCropped<- 
  extend(mask(crop(stateclasses, studyExtent), onlyMonteregie), studyExtent)


writeRaster(primaryStratumCropped, 
            "Data/stsim/primary_stratum_FocalAreaBuffer.tif",
            overwrite = TRUE)
writeRaster(secondaryStratumCropped, 
            "Data/stsim/secondary_stratum_FocalAreaBuffer.tif",
            overwrite = TRUE)
writeRaster(tertiaryStratumCropped, 
            "Data/stsim/tertiary_stratum_FocalAreaBuffer.tif",
            overwrite = TRUE)
writeRaster(stateclassesCropped, 
            "Data/stsim/stateclasses_FocalAreaBuffer.tif",
            overwrite = TRUE)
