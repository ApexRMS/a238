#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Run Zonation prioritizations 
# 12-2020                                       					
# 
#	  Inputs:                                
#    - Zonation run settings template
# 
#   Outputs:
#     - create Zonation input files:
#     - Run settings (.dat), Biodiversity feature list (.spp), batch file (.bat)
#     - Solution maps for 8 scenarios (2 cell removal rules * 4 conservation criteria combos)
#
# Script by Bronwyn Rayfield for ApexRMS 									
#############################################################################################

# Workspace ---------------------------------------------------------
# Packages
library(tidyverse)
library(raster)

# Directories
## Directories
rawDataDir <- paste0("Data/Raw")
procDataDir <- paste0("Data/Processed")
zonationInputDir <- paste0("Data/Processed/PrioritizationInputs/Zonation")
zonationOutputDir <- paste0("Outputs/PrioritizationSolutions/Zonation")
zonationPath<-"C:/Program Files/zonation 4.0.0rc1_compact/bin/zig4.exe"

# Cell removal rules
cellRemovalRuleNameList <- c("CAZ", "ABF")
cellRemovalRuleList <- c(1, 2)

# Ecoregions - zone1 Adirondacks, zone 3 = StL lowlands, Zone 4 = appalachians
ecoregionList <- c(1, 4, 3)

# Conservation criteria to be included in prioritization analyses
criteriaList <- c("SuitabilityR", "Density", "AreaR")

# Area targets to convert rank solution to binary solution
areaTargetList <- c(0.05, 0.1, 0.17)
targetSummary <- read_csv(file.path(procDataDir, "TargetAreaSummary - Caroline1.csv"))


# Load files ---------------------------------------------------------
# Load species list
speciesID <- read.csv(
  file.path(
    paste0(rawDataDir, "/Focal Species"), 
    "Species.csv"), 
  stringsAsFactors = FALSE)
speciesList <- speciesID$Code

# Unprotected natural areas - to be used to crop Zonation outputs to original extent
unprotectedNaturalAreasFocal <- raster(file.path(procDataDir, "unprotectedNaturalAreas_Focal.tif"))

# Zonation Run Settings template
zonationSetTemplate <- readLines(file.path("config", "zonation", "ALL_set_template.dat"))


# Run Zonation  ------------------------
# Loop over all cell removal rules
for(i in 1:length(cellRemovalRuleNameList)){
  cellRemovalRuleName <- cellRemovalRuleNameList[i]
  cellRemovalRule <- cellRemovalRuleList[i]
  
  # Loop over all ecoregions
  for(ecoregionID in ecoregionList){ 
    # Create run settings file  
    zonationSet <- zonationSetTemplate
    zonationSet[2] <- paste("removal rule =", cellRemovalRule)
    zonationSet[12] <- "use mask = 1"
    zonationSet[15] <- paste("mask file =", file.path(zonationInputDir, paste0("zonation_ProtectedAreas_ecoregion", ecoregionID, ".tif")))
    zonationSet[30] <- "mask missing areas = 1"
    zonationSet[31] <- paste("area mask file =", file.path(zonationInputDir, paste0("zonation_UnprotectedNaturalAreas_ecoregion", ecoregionID, ".tif")))
    # Name the run settings file tagged by ecoregion
    zonationSetName <- file.path(getwd(), zonationInputDir, paste0("RunSetting_", cellRemovalRuleName, "_ecoregion", ecoregionID, ".dat"))
    writeLines(zonationSet, zonationSetName)
    
    # Loop over all conservation criteria individually and then all together ("All") 
    # to create biodiversity feature list (.spp) and batch file (.bat) 
    for(criterion in c(criteriaList, "All")){
      # Create names of biodiversity feature maps for each species
      # Include only a single conservation criterion
      if(criterion != "All"){
        sppList <- paste0(speciesList, "_", criterion, "_ecoregion", ecoregionID, ".tif")
      }
      # Include all conservation criteria together
      if(criterion == "All"){
        sppList <- c(paste0(speciesList, "_", criteriaList[1], "_ecoregion", ecoregionID, ".tif"),
                    paste0(speciesList, "_", criteriaList[2], "_ecoregion", ecoregionID, ".tif"),
                    paste0(speciesList, "_", criteriaList[3], "_ecoregion", ecoregionID, ".tif"))
      }
      # If using CAZ or ABF cell removal rule then use the same parameters in the biodiversity feature list
      if(cellRemovalRule != 4){
        zonationSpp <- data.frame(1, 0, 1, 1, 1, file.path(getwd(), zonationInputDir, sppList))
      }
      # If using GBF cell removal rule then specify extra parameters in the biodiversity feature list
      if(cellRemovalRule == 4){
        zonationSpp <- data.frame(1, 0, 1, 1, 1, file.path(getwd(), zonationInputDir, sppList))
      }
      # Name the biodiversity feature list file tagged with cell removal rule, conservation criteria, and ecoregion
      zonationSppName <- file.path(getwd(), zonationInputDir, paste0("BiodiversityFeatureList_", cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".spp"))
      write.table(zonationSpp, zonationSppName, row.names = FALSE, col.names=FALSE)
      
      # Create batch file to run Zonation for this cell removal rule, criteria, and ecoregion
      zonationOutName <- file.path(getwd(), zonationOutputDir, paste0(cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".txt"))
      zonationScenario <- paste0("-r ", "\"", zonationSetName, "\" ", "\"", zonationSppName, "\" ", "\"", zonationOutName, "\" ", "0.0 0 1.0 0")
      zonationScenarioName <- file.path(getwd(), zonationInputDir, paste0("ZonationScenario_", cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".bat"))
      writeLines(zonationScenario, zonationScenarioName)
      
      # Run Zonation
      system(paste0("\"",zonationPath, "\" ", zonationScenario))
    } # conservation criteria
  } # ecoregion
} # cell removal rule


# Loop over all cell removal rules
for(i in 1:length(cellRemovalRuleNameList)){
  cellRemovalRuleName <- cellRemovalRuleNameList[i]
  cellRemovalRule <- cellRemovalRuleList[i]
  # Loop through all conservation criteria individually and all together
  # Loop through all area target levels
  for(criterion in c(criteriaList, "All")){
    # Read 3 ecoregion files
    er1 <- raster(file.path(zonationOutputDir, paste0(cellRemovalRuleName, "_", criterion, "_ecoregion1.rank.compressed.tif")))
    er3 <- raster(file.path(zonationOutputDir, paste0(cellRemovalRuleName, "_", criterion, "_ecoregion3.rank.compressed.tif")))
    er4 <- raster(file.path(zonationOutputDir, paste0(cellRemovalRuleName, "_", criterion, "_ecoregion4.rank.compressed.tif")))
    
    for(areaTarget in areaTargetList){  
      # Create binary solutions corresponding to area target
      er1Target <- Which(er1 >= (1-(pull(targetSummary %>% filter(ID==1) %>% dplyr::select(paste0('TargetNaturalAreaPercent', areaTarget)))/100)))
      er3Target <- Which(er3 >= (1-pull(targetSummary %>% filter(ID==3) %>% dplyr::select(paste0('TargetNaturalAreaPercent', areaTarget)))/100))
      er4Target <- Which(er4 >= (1-pull(targetSummary %>% filter(ID==4) %>% dplyr::select(paste0('TargetNaturalAreaPercent', areaTarget)))/100))
      
      # Combine all ecoregion solutions together
      erAllTarget <- mosaic(er3Target, er4Target, fun="max", na.rm=TRUE) %>%
        mosaic(., er1Target,  fun="max", na.rm=TRUE) %>%
        crop(., unprotectedNaturalAreasFocal) %>%
        mask(., unprotectedNaturalAreasFocal)
      
      # Save one solution for each cell removal rule, conservation criteria, and target level
      writeRaster(erAllTarget, file.path(zonationOutputDir, paste0(cellRemovalRuleName, "_", criterion, "_", areaTarget, ".tif")), overwrite=T)
    } # conservation criteria
  } # target level
} # cell removal rule


# Not used
# # Generalized benefit function 
# cellRemovalRuleName <- "GBF"
# cellRemovalRule <- 4
# relativeTargetPerFeature <- 0.75
# 
# for(ecoregionID in ecoregionList){ 
#   
#   zonationSet <- readLines(file.path("config", "zonation", "ALL_set_template.dat"))
#   zonationSet[12] <- paste("removal rule =", cellRemovalRule)
#   zonationSet[2] <- "use mask = 1"
#   zonationSet[15] <- paste("mask file =", file.path(zonationInputDir, paste0("zonation_ProtectedAreas_ecoregion", ecoregionID, ".tif")))
#   zonationSet[30] <- "mask missing areas = 1"
#   zonationSet[31] <- paste("area mask file =", file.path(zonationInputDir, paste0("zonation_NaturalAreas_ecoregion", ecoregionID, ".tif")))
#   
#   zonationSetName <- file.path(getwd(), "config", "zonation", paste0("RunSetting_ecoregion", ecoregionID, ".dat"))
#   writeLines(zonationSet, zonationSetName)
#   
#   # Loop through all conservation criteria individually
#   for(criterion in criteriaList){
#     sppList <- paste0(speciesList, "_", criterion, "_ecoregion", ecoregionID, ".tif")
#     zonationSpp <- data.frame(1, 0, 1, 1, 0, relativeTargetPerFeature, 1, 1, file.path(getwd(), zonationInputDir, sppList))
#     zonationSppName <- file.path(getwd(), "config", "zonation", paste0("BiodiversityFeatureList_", cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".spp"))
#     write.table(zonationSpp, zonationSppName, row.names = FALSE, col.names=FALSE)
#     
#     
#     zonationOutName <- file.path(getwd(), "config", "zonation", paste0(cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".txt"))
#     zonationScenario <- paste0("-r ", "\"", zonationSetName, "\" ", "\"", zonationSppName, "\" ", "\"", zonationOutName, "\" ", "0.0 0 1.0 0")
#     zonationScenarioName <- file.path("config", "zonation", paste0("ZonationScenario_", cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".bat"))
#     writeLines(zonationScenario, zonationScenarioName)
#     
#     #run Zonation
#     system(paste0("\"",zonationPath, "\" ", zonationScenario))
#   }
#   # Run all conservation criteria together
#   criterion <- "All"
#   sppList <- c(paste0(speciesList, "_", criteriaList[1], "_ecoregion", ecoregionID, ".tif"),
#                paste0(speciesList, "_", criteriaList[2], "_ecoregion", ecoregionID, ".tif"),
#                paste0(speciesList, "_", criteriaList[3], "_ecoregion", ecoregionID, ".tif"))
#   zonationSpp <- data.frame(1, 0, 1, 1, 0, relativeTargetPerFeature, 1, 1, file.path(getwd(), zonationInputDir, sppList))
#   zonationSppName <- file.path(getwd(), "config", "zonation", paste0("BiodiversityFeatureList_", cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".spp"))
#   write.table(zonationSpp, zonationSppName, row.names = FALSE, col.names=FALSE)
#   
#   
#   zonationOutName <- file.path(getwd(), "config", "zonation", paste0(cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".txt"))
#   zonationScenario <- paste0("-r ", "\"", zonationSetName, "\" ", "\"", zonationSppName, "\" ", "\"", zonationOutName, "\" ", "0.0 0 1.0 0")
#   zonationScenarioName <- file.path("config", "zonation", paste0("ZonationScenario_", cellRemovalRuleName, "_", criterion, "_ecoregion", ecoregionID, ".bat"))
#   writeLines(zonationScenario, zonationScenarioName)
#   
#   #run Zonation
#   system(paste0("\"",zonationPath, "\" ", zonationScenario))
# }


# End script