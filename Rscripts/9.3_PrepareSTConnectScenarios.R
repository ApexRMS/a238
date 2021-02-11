# Program the spatial multipliers scenarios

# Packages
library(tidyverse)
library(rsyncrosim)
library(raster)
library(stringr)

# Lib
mylib <- ssimLibrary("Simulations/Monteregie_stconnect.ssim")
myproj <- project(mylib, "Definitions")

# MSC methods
mscList <- c("Outputs/PrioritizationSolutions/MapsBudget0.05/CAZ_All_0.05.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.05/Sum-Species-All0.05.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.17/CAZ_All_0.17.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.17/Sum-Species-All0.17.tif")

# Reclassficfication matrix
reclassMatrix <- matrix(c(0,1,1,0), 2, 2)

# Read the current multpliers
protectedAreas <- raster("Data/Processed/spatialMultiplier_ProtectedAreas_FocalArea.tif")

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
listOfProcessedMultpliers <- lapply(mscList, processMultipliers, 
                                    reclassMatrix = reclassMatrix,
                                    match_combine = protectedAreas)
names(listOfProcessedMultpliers) <- unlist(lapply(listOfProcessedMultpliers , names))

# Write the files out
listOfFilesToWrite <- paste0("Outputs/SpatialMultipliers/", names(listOfProcessedMultpliers), ".tif")
mapply(listOfProcessedMultpliers, listOfFilesToWrite, FUN=writeRaster, overwrite=TRUE)

# Create subscenario for each, then combine these with the other scenarios
# Here we only do 2 as a test

for (mult in listOfProcessedMultpliers){
  mySce <- scenario(myproj, names(mult), overwrite = TRUE)
  theDatasheet <- 
    data.frame(TransitionGroupID = c("Agricultural Expansion and Urbanization"), 
               MultiplierFileName = paste0(getwd(),"/Outputs/SpatialMultipliers/", 
                                           names(mult), ".tif"))
  saveDatasheet(ssimObject = mySce, data = theDatasheet, name = 'stsim_TransitionSpatialMultiplier')
}

# Now combine these subscenario with the climate and land use in order to produce full scenarios

# First step: list scenarios that are common in all scnearios
# WARNING: This assumes we are using targets given the fact that we did not 
#          make the multipliers work.
commonSet <- list(RunControl = 220, 
                  IniCond = 247, 
                  OutOpt = 8, 
                  TrAdj = 13, 
                  StAtrr = 14, 
                  TrSiz = 15, 
                  TSTr = 16)
species14Set <- c(HabSuit = 201, HabPat = 203, Res = 204)

# Second, the 2 land use change scenarios
lucSet <- list(noLULC = 6, historicLULC = 7)

# Third, the climate change scenarios
ccSet <- list(baseline = 9, rcp85 = 98) #rcp45 = 97, 

# Loop and create the full scenarios

for (mult in listOfProcessedMultpliers){
  for (lu in names(lucSet)){
    for (cc in names(ccSet)){
      
      # Create name
      scenarioName <- paste(names(mult), lu, cc,  sep= "_")
      print(scenarioName)
      
      # Create the scenario
      tempSce <- scenario(myproj, scenarioName, overwrite = TRUE)
      
      # Set dependencies
      dependency(tempSce, unlist(c(commonSet, species14Set, 
                                   lucSet[[lu]], ccSet[[cc]])))
    }
  }
}
