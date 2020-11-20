
library(raster)
library(sf)

primaryStratum <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/PrimaryStratum_NAchanged.tif")
secondaryStratum <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/SecondaryStratum_NAchanged.tif")
tertiaryStratum <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-127/stsim_InitialConditionsSpatial/TertiaryStratum_NAchanged.tif")
protectedAreas <- 
  raster("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim.input/Scenario-148/stsim_TransitionSpatialMultiplier/protectedAreas_FocalAreaBuffer.tif")

studyExtent <- raster("Data/Processed/LULC_FocalAreaBuffer.tif")
onlyMonteregie <- st_transform(st_read("Data/Processed/regioMonteregie.shp"),
                               crs(studyExtent))

primaryStratumCropped <- 
  extend(mask(crop(primaryStratum, studyExtent), onlyMonteregie), studyExtent)
primaryStratumCropped[is.na(primaryStratumCropped)] <- 0
primaryStratumCropped[is.na(studyExtent)] <- NA
secondaryStratumCropped<- 
  extend(mask(crop(secondaryStratum, studyExtent, snap="out"), onlyMonteregie), studyExtent)
secondaryStratumCropped[is.na(secondaryStratumCropped)] <- 0
secondaryStratumCropped[is.na(studyExtent)] <- NA
tertiaryStratumCropped<- 
  extend(mask(crop(tertiaryStratum, studyExtent), onlyMonteregie), studyExtent)
tertiaryStratumCropped[is.na(tertiaryStratumCropped)] <- 0
tertiaryStratumCropped[is.na(studyExtent)] <- NA
protectedAreasCropped <- 
  extend(mask(crop(protectedAreas, studyExtent), onlyMonteregie), studyExtent)
protectedAreasCropped[is.na(protectedAreasCropped)] <- 0
protectedAreasCropped[is.na(studyExtent)] <- NA

writeRaster(primaryStratumCropped, 
            "Data/stsim/primary_stratum_FocalAreaBuffer.tif",
            overwrite = TRUE)
writeRaster(secondaryStratumCropped, 
            "Data/stsim/secondary_stratum_FocalAreaBuffer.tif",
            overwrite = TRUE)
writeRaster(tertiaryStratumCropped, 
            "Data/stsim/tertiary_stratum_FocalAreaBuffer.tif",
            overwrite = TRUE)
writeRaster(studyExtent, 
            "Data/stsim/stateclasses_FocalAreaBuffer.tif",
            overwrite = TRUE)
writeRaster(protectedAreasCropped, 
            "Data/stsim/protectedAreas_FocalAreaBuffer.tif",
            overwrite = TRUE)

## Filter the targets & multipliers

library(rsyncrosim)
mylib <- ssimLibrary("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim")

sceTar <- scenario(mylib, 7)
myproj <- project(mylib, "Definitions")

secondaryStratumDatasheet <- datasheet(myproj, "stsim_SecondaryStratum")
goodValues <- 
  as.numeric(names(which(table(values(secondaryStratumCropped))>1000)))
secondaryStratumDatasheetFiltered <- secondaryStratumDatasheet %>% 
  dplyr::filter(ID %in% goodValues)

# saveDatasheet(ssimObject = myproj, data = secondaryStratumDatasheetFiltered, 
#               name = "stsim_SecondaryStratum")

targets <- datasheet(sceTar, "stsim_TransitionTarget")
targetsFiltered <- targets %>% 
  dplyr::filter(SecondaryStratumID %in% secondaryStratumDatasheetFiltered$Name)

saveDatasheet(sceTar, targetsFiltered, "stsim_TransitionTarget")

sceMult <- scenario(mylib, 9)
sceTrans <- scenario(mylib, 4)

tertiaryStratumDatasheet <- datasheet(sceMult, "stsim_TertiaryStratum")
tertiaryStratumDatasheetFiltered <- tertiaryStratumDatasheet %>% 
  dplyr::filter(ID %in% unique(values(tertiaryStratumCropped)))

multipliers <- datasheet(sceMult, "stsim_TransitionMultiplierValue")

trans <- datasheet(sceTrans, "stsim_Transition")
uniqueTrans <- paste0(trans$TransitionTypeID, " [Type]")

multipliersFiltered <- multipliers %>% 
  dplyr::filter(TertiaryStratumID %in% tertiaryStratumDatasheetFiltered$Name)
addMultipliers <- data.frame(Amount = 0, StratumID= "NotStudyRegion", 
                             TransitionGroupID = uniqueTrans)
multipliersFilteredAugmented <- dplyr::bind_rows(multipliersFiltered, addMultipliers)

saveDatasheet(sceMult, multipliersFilteredAugmented, "stsim_TransitionMultiplierValue")
