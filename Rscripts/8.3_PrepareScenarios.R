# Program the spatial multipliers scenarios

# Packages
library(tidyverse)
library(rsyncrosim)
library(raster)
library(stringr)

# Lib
mylib <- ssimLibrary("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim")
myproj <- project(mylib, "Definitions")

# Input folders for priorization files
inputFolders <- c("Results/MapsBudget0.17/", "Results/MapsBudget0.1/")

# Read list of files 
listOfFiles <- unlist(lapply(inputFolders, list.files, full.names=T), recursive = F)

# Reclassficfication matrix
reclassMatrix <- matrix(c(0,1,1,0), 2, 2)

# Read the current multpliers
protectedAreas <- raster("Data/stsim/protectedAreas_FocalAreaBuffer.tif")

# Process the multipliers -------------------------------------------------

# Function that process multipliers
processMultipliers <- function(file_name, 
                               reclassMatrix = reclassMatrix,
                               match_combine = protectedAreas){
  # Load
  multiplier <- raster(file_name)
  # Reformat name
  names(multiplier) <- tools::file_path_sans_ext(basename(file_name)) %>% 
    str_replace_all("-", "_")
  # reclassify, extend
  multiplierUpdated <- 
    extend(reclassify(multiplier, reclassMatrix), protectedAreas)
  # Where there are protected areas, and in the buffer, set multiplier to 0
  multiplierUpdated[protectedAreas==0] <- 0
  # Everywhere else set the value to 1
  multiplierUpdated[!is.na(protectedAreas) & is.na(multiplierUpdated)] <- 1
  multiplierUpdated
}

# Process all files
listOfProcessedMultpliers <- lapply(listOfFiles, processMultipliers, 
                                    reclassMatrix = reclassMatrix,
                                    match_combine = protectedAreas)
names(listOfProcessedMultpliers) <- unlist(lapply(listOfProcessedMultpliers , names))

# Write the files out
listOfFilesToWrite <- paste0("Results/MultipliersProcessed/", names(listOfProcessedMultpliers), ".tif")
mapply(listOfProcessedMultpliers, listOfFilesToWrite, FUN=writeRaster, overwrite=TRUE)

# Create subscenario for each, then combine these with the other scenarios
# Here we only do 2 as a test
listOfProcessedMultpliersSUBSET <- listOfProcessedMultpliers[1:2]

for (mult in listOfProcessedMultpliersSUBSET){
  mySce <- Scenario(myproj, names(mult))
  theDatasheet <- 
    data.frame(TransitionGroupID = c("Urbanisation"), 
               MultiplierFileName = paste0(getwd(),"Results/MultipliersProcessed/", 
                                           names(mult), ".tif"))
  saveDatasheet(ssimObject = mySce, data = theDatasheet, name = 'stsim_TransitionSpatialMultiplier')
}

# Now combine these subscenario with the climate and land use in order to produce full scenarios
## TODO