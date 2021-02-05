## A238 Preparing Inputs for Monteregie library
## V.L.

library(raster)
library(sf)
library(rsyncrosim)

# Prepare spatial inputs --------------------------------------------------

# Load the original data
primaryStratum <- 
  raster("Data/stsim/PrimaryStratum_NAchanged.tif")
secondaryStratum <- 
  raster("Data/stsim/SecondaryStratum_NAchanged.tif")
tertiaryStratum <- 
  raster("Data/stsim/TertiaryStratum_NAchanged.tif")
protectedAreas <- 
  raster("Data/stsim/protectedAreas_FocalAreaBuffer.tif")

# Load the study region data
studyExtent <- raster("Data/Processed/LULC_FocalAreaBuffer.tif")
onlyMonteregie <- st_transform(st_read("Data/Processed/regioMonteregie.shp"),
                               crs(studyExtent))

# Resize and clean secondary stratum
secondaryStratumCropped<- 
  extend(mask(crop(secondaryStratum, studyExtent, snap="out"), onlyMonteregie), studyExtent)
secondaryStratumCropped[is.na(secondaryStratumCropped)] <- 0
goodValues <- 
  as.numeric(names(which(table(values(secondaryStratumCropped))>1000)))
secondaryStratumCropped[!(secondaryStratumCropped %in% goodValues)] <- 0
secondaryStratumCropped[is.na(studyExtent)] <- NA

# Cleaning
secondaryStratumCroppedClean <- secondaryStratumCropped
secondaryStratumCroppedClean[secondaryStratumCroppedClean!=0] <- 1
secondaryStratumCroppedClean[secondaryStratumCroppedClean==0] <- NA

# Create mask
writeRaster(secondaryStratumCroppedClean, "Data/Processed/Monteregie_mask.tif")

# Do the same with primary stratum
primaryStratumCropped <- 
  extend(mask(crop(primaryStratum, studyExtent), onlyMonteregie), secondaryStratumCroppedClean)
primaryStratumCropped[is.na(secondaryStratumCroppedClean)] <- 0
primaryStratumCropped[is.na(studyExtent)] <- NA

# And with the tertiary stratum
tertiaryStratumCropped<- 
  extend(mask(crop(tertiaryStratum, studyExtent), onlyMonteregie), secondaryStratumCroppedClean)
tertiaryStratumCropped[is.na(secondaryStratumCroppedClean)] <- 0
tertiaryStratumCropped[is.na(studyExtent)] <- NA
protectedAreasCropped <- 
  extend(mask(crop(protectedAreas, studyExtent), onlyMonteregie), secondaryStratumCroppedClean)
protectedAreasCropped[is.na(secondaryStratumCroppedClean)] <- 0
protectedAreasCropped[is.na(studyExtent)] <- NA

# Write out the files
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

# Filter Targets ----------------------------------------------------------

mylib <- ssimLibrary("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim")
myproj <- project(mylib, "Definitions")

# Collect which secondary stratums are in Monteregie + filter tertiary
secondaryStratumDatasheet <- datasheet(myproj, "stsim_SecondaryStratum")
secondaryStratumDatasheetFiltered <- secondaryStratumDatasheet %>% 
  dplyr::filter(ID %in% goodValues)

tertiaryStratumDatasheet <- datasheet(myproj, "stsim_TertiaryStratum")
tertiaryStratumDatasheetFiltered <- tertiaryStratumDatasheet %>% 
  dplyr::filter(ID %in% unique(values(tertiaryStratumCropped)))

# Historic
sceTar <- scenario(mylib, 7)

targets <- datasheet(sceTar, "stsim_TransitionTarget")
targetsFiltered <- targets %>% 
  dplyr::filter(SecondaryStratumID %in% secondaryStratumDatasheetFiltered$Name)

saveDatasheet(sceTar, targetsFiltered, "stsim_TransitionTarget")

#  Filter Multipliers And Targets -----------------------------------------

# Pick up unique transitions
sceTrans <- scenario(mylib, 4)
trans <- datasheet(sceTrans, "stsim_Transition")
uniqueTrans <- paste0(trans$TransitionTypeID, " [Type]")

# Generate multipliers for the buffer
addMultipliers <- data.frame(Amount = 0, StratumID= "NotStudyRegion", 
                             TransitionGroupID = uniqueTrans)

# Baseline
sceMult <- scenario(mylib, 9)
multipliers <- datasheet(sceMult, "stsim_TransitionMultiplierValue")

multipliersFiltered <- multipliers %>% 
  dplyr::filter(TertiaryStratumID %in% tertiaryStratumDatasheetFiltered$Name)
multipliersFilteredAugmented <- dplyr::bind_rows(multipliersFiltered, addMultipliers)

saveDatasheet(sceMult, multipliersFilteredAugmented, "stsim_TransitionMultiplierValue")

# RCP 4.5
sceMult45 <- scenario(mylib, 97)
# multipliers45 <- datasheet(sceMult45, "stsim_TransitionMultiplierValue") # This doesnt work!
multipliers45 <- read_csv("config/stsim/TransitionMultipliers45.csv")

multipliers45Filtered <- multipliers45 %>% 
  dplyr::filter(TertiaryStratumID %in% tertiaryStratumDatasheetFiltered$Name)
multipliers45FilteredAugmented <- dplyr::bind_rows(multipliers45Filtered, addMultipliers)

saveDatasheet(sceMult45, multipliers45FilteredAugmented, "stsim_TransitionMultiplierValue")

write.csv(multipliers45FilteredAugmented, "config/stsim/TransitionMultipliers45_new.csv")

# RCP 8.5
sceMult85 <- scenario(mylib, 98)
# multipliers85 <- datasheet(sceMult85, "stsim_TransitionMultiplierValue") # This doesnt work!
multipliers85 <- read_csv("config/stsim/TransitionMultipliers85.csv")

multipliers85Filtered <- multipliers85 %>% 
  dplyr::filter(TertiaryStratumID %in% tertiaryStratumDatasheetFiltered$Name)
multipliers85FilteredAugmented <- dplyr::bind_rows(multipliers85Filtered, addMultipliers)

saveDatasheet(sceMult85, multipliers85FilteredAugmented, "stsim_TransitionMultiplierValue")

write.csv(multipliers85FilteredAugmented, "config/stsim/TransitionMultipliers85_new.csv")
