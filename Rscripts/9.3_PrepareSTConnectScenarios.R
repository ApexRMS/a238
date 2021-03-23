# Program the spatial multipliers scenarios

# Packages
library(tidyverse)
library(rsyncrosim)
library(raster)
library(stringr)

# Lib
mySession <- "F:/2.2.7/SyncroSim.Console.exe"
mylib <- ssimLibrary("Simulations/Monteregie_stconnect.ssim", session = mySession)
myproj <- project(mylib, "Definitions")

# MSC methods
mscList <- c("Outputs/PrioritizationSolutions/MapsBudget0.1/Generic-Resistance0.1.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.1/Ecoprofile-Trophic0.1.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.1/CAZ_Density_0.1.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.1/Sum-Species-Density0.1.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.1/CAZ_All_0.1.tif",
             "Outputs/PrioritizationSolutions/MapsBudget0.1/Sum-Species-All0.1.tif")

"Outputs/PrioritizationSolutions/MapsBudget0.1/Minimize-Shortfall-Species-All0.1.tif", 

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
for (mult in names(listOfProcessedMultpliers)){
  mySce <- scenario(myproj, mult, overwrite = TRUE)
  theDatasheet <- 
    data.frame(TransitionGroupID = c("Agricultural Expansion and Urbanization"), 
               MultiplierFileName = paste0(getwd(),"/Outputs/SpatialMultipliers/", 
                                           mult, ".tif"))
  saveDatasheet(ssimObject = mySce, data = theDatasheet, name = 'stsim_TransitionSpatialMultiplier')
}

# Get scenario id's for multipliers
multSet <- list()
for (mult in names(listOfProcessedMultpliers)){
  multSet<-c(multSet, scenarioId(scenario(myproj, mult)))
}
names(multSet) <- names(listOfProcessedMultpliers)

# Now combine these subscenario with the climate and land use in order to produce full scenarios

# First step: list scenarios that are common in all scnearios
commonSet <- list(RunControl = 370, 
                  IniCond = 247, 
                  OutOpt = 8,
                  TrPath = 4,
                  TrAdj = 13, 
                  StAtrr = 14,
                  TrAttr = 321,
                  TrSiz = 15, 
                  TSTr = 16)
species14Set <- c(HabSuit = 201, HabPat = 203, Res = 204)

# Second, the 2 land use change scenarios
lucSet <- list(historicLULC = 7) #noLULC = 6, 

# Third, the climate change scenarios
ccSet <- list(baseline = 9) #rcp45 = 97, rcp85 = 98) 

# Loop and create the full scenarios

for (mult in names(multSet)){
  for (lu in names(lucSet)){
    for (cc in names(ccSet)){
      
      # Create name
      scenarioName <- paste(mult, lu, cc,  sep= "_")
      print(scenarioName)
      
      # Create the scenario
      tempSce <- scenario(myproj, scenarioName, overwrite = TRUE)
      
      # Set dependencies
      dependency(tempSce, unlist(c(commonSet, species14Set, 
                                   lucSet[lu], ccSet[cc], multSet[mult])))
    }
  }
}


# Loop and run the full scenarios
for (mult in names(listOfProcessedMultpliers)){
  for (lu in names(lucSet)){
    for (cc in names(ccSet)){
      
      # Create name
      scenarioName <- paste(mult, lu, cc,  sep= "_")
      print(scenarioName)

      run(myproj, scenario = scenarioName, jobs = 5)      
    }
  }
}


# Baseline scenario
# Create the scenario
scenarioName <- "PA2010_historicLULC_baseline"
tmpSce <- scenario(myproj, scenarioName, overwrite = TRUE)
PA2010scenarioID <- list(SpTrMult = 248)
# Set dependencies
dependency(tmpSce, unlist(c(commonSet, species14Set, 
                                  lucSet[lu], ccSet[cc], PA2010scenarioID)))
# Run
run(myproj, scenario = scenarioName, jobs = 3)

# # No longer used
# # Set up to run habitat sitability transformer separately
# for (mult in listOfProcessedMultpliers){
#   for (lu in names(lucSet)){
#     for (cc in names(ccSet)){
#       
#       # Create name
#       scenarioName <- paste("Habitat", names(mult), lu, cc,  sep= "_")
#       print(scenarioName)
#       
#       # Create the scenario
#       tempSce <- scenario(myproj, scenarioName, overwrite = TRUE)
#       
#       # Set dependencies
#       dependency(tempSce, c(scenarioName, "Habitat Suitability: Basic Landcover and Forest Age - 14 species"))
#     }
#   }
# }
# 
# for (mult in listOfProcessedMultpliers){
#   for (lu in names(lucSet)){
#     for (cc in names(ccSet)){
#       
#       # Create name
#       scenarioName <- paste("Habitat", names(mult), lu, cc,  sep= "_")
#       print(scenarioName)
#       
#       run(myproj, scenario = scenarioName, jobs = 3, transformerName = "stconnect_HabitatSuitability")      
#     }
#   }
# }
# 

run(myproj, scenario = "PA2010_historicLULC_baseline", jobs = 3, transformerName = "stconnect_STSimWrapper")
