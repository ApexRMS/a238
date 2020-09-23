library(raster)
library(rsyncrosim)

e = ssimEnvironment()

source(file.path(e$PackageDirectory, "common.R"))

#datasheets
stateClassIn = datasheet(GLOBAL_Scenario,"stsim_StateClass")
names(stateClassIn) <- c("StateClassID", "StateLabelXID", "StateLabelYID", "ID", "Color")
habitatSuitabilityIn = GetDataSheetExpectData("stconnect_HSHabitatSuitability", GLOBAL_Scenario)
habitatSuitabilityIn <- merge(habitatSuitabilityIn,stateClassIn,by="StateClassID")
habitatPatchIn = GetDataSheetExpectData("stconnect_HSHabitatPatch", GLOBAL_Scenario)
resistanceIn = GetDataSheetExpectData("stconnect_HSResistance", GLOBAL_Scenario)
resistanceIn <- merge(resistanceIn, stateClassIn, by="StateClassID")
habitatSuitabilityOut = datasheet(GLOBAL_Scenario, "stconnect_HSOutputHabitatSuitability")
habitatPatchOut = datasheet(GLOBAL_Scenario, "stconnect_HSOutputHabitatPatch")
resistanceOut = datasheet(GLOBAL_Scenario, "stconnect_HSOutputResistance")

#file paths
tempFolderPath = envTempFolder("HabitatSuitability")

# Temporal resolution of analysis in years
# e.g. analyse every 10 years
temporalRes <- 20
# Set of timesteps to analyse
timestepSet <- seq(GLOBAL_MinTimestep, GLOBAL_MaxTimestep, by=temporalRes)

#Simulation
envBeginSimulation(GLOBAL_TotalIterations * GLOBAL_TotalTimesteps)

for (iteration in GLOBAL_MinIteration:GLOBAL_MaxIteration) {

    for (timestep in timestepSet) {

        envReportProgress(iteration, timestep)

        for (sprow in 1:nrow(GLOBAL_Species)) {

            species = GLOBAL_Species[sprow, "Name"]
            speciesCode = GLOBAL_Species[sprow, "Code"]
            suitabilityValues <- habitatSuitabilityIn[which(habitatSuitabilityIn$SpeciesID == species), c("ID", "Value")]
            patchValues <- habitatPatchIn[which(habitatPatchIn$SpeciesID == species), c("HabitatSuitabilityThreshold", "MinimumHabitatPatchSize")]
            resistanceValues <- resistanceIn[which(resistanceIn$SpeciesID == species), c("ID", "Value")]

            if (nrow(suitabilityValues) == 0 || nrow(patchValues) == 0 || nrow(resistanceValues) == 0) {
                next
            }

            suitabilityThreshold <- patchValues$HabitatSuitabilityThreshold
            patchSizeThreshold <- patchValues$MinimumHabitatPatchSize

            #Input stateclass raster
            stateMap <- datasheetRaster(GLOBAL_Scenario, datasheet = "stsim_OutputSpatialState", iteration = iteration, timestep = timestep)

            #Habitat suitability map
            suitabilityRaster <- reclassify(stateMap, rcl = suitabilityValues)

            #Habitat patch map
            habitatRaster <- Which(suitabilityRaster >= suitabilityThreshold)
            #convert from hectares to m
            conversionFromHa <- res(habitatRaster)[1]*res(habitatRaster)[2]*(1/10000)
            habitatClump <- clump(habitatRaster)
            habitatClumpID <- data.frame(freq(habitatClump))
            #Remove clump observations with frequency smaller than minimum habitat patch size (ha)
            habitatClumpID <- habitatClumpID[habitatClumpID$count < patchSizeThreshold/conversionFromHa,]
            habitatRaster[Which(habitatClump %in% habitatClumpID$value)] <- 0

            #Resistance map
            #Reclassify based on resistance values
            resistanceRasterReclass <- reclassify(stateMap, rcl = resistanceValues)
            #Overlay habitat patches
            resistanceRasterOverlay <- overlay(resistanceRasterReclass, habitatRaster, fun = function(x,y){return(x+y)})
            #Reclass to assign habitat patches a resistance value = 1 (note that both overlaid values of both 3 and 5 correspond to habitat patches)
            resistanceRaster <- reclassify(resistanceRasterOverlay, rcl = matrix(c(3, 1, 5, 1), ncol=2, byrow = T))
            
            #Save rasters
            suitabilityName = file.path(tempFolderPath, CreateRasterFileName2("HabitatSuitability", speciesCode, iteration, timestep, "tif"))
            patchName = file.path(tempFolderPath, CreateRasterFileName2("HabitatPatch", speciesCode, iteration, timestep, "tif"))
            resistanceName = file.path(tempFolderPath, CreateRasterFileName2("Resistance", speciesCode, iteration, timestep, "tif"))

            writeRaster(suitabilityRaster, suitabilityName, overwrite = TRUE)
            df = data.frame(Iteration = iteration, Timestep = timestep, SpeciesID = species, Filename = suitabilityName)
            habitatSuitabilityOut = addRow(habitatSuitabilityOut, df)

            writeRaster(habitatRaster, patchName, overwrite = TRUE)
            df = data.frame(Iteration = iteration, Timestep = timestep, SpeciesID = species, Filename = patchName)
            habitatPatchOut = addRow(habitatPatchOut, df)

            writeRaster(resistanceRaster, resistanceName, overwrite = TRUE)
            df = data.frame(Iteration = iteration, Timestep = timestep, SpeciesID = species, Filename = resistanceName)
            resistanceOut = addRow(resistanceOut, df)
        }

        envStepSimulation()
    }
}

saveDatasheet(GLOBAL_Scenario, habitatSuitabilityOut, "stconnect_HSOutputHabitatSuitability")
saveDatasheet(GLOBAL_Scenario, habitatPatchOut, "stconnect_HSOutputHabitatPatch")
saveDatasheet(GLOBAL_Scenario, resistanceOut, "stconnect_HSOutputResistance")

envEndSimulation()
