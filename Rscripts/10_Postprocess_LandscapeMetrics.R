library(tidyverse)
library(rsyncrosim)
library(raster)
library(sf)
library(viridis)
library(scales)
library(landscapemetrics)
library(cowplot)
library(ggthemes)

options(tibble.width = Inf, tibble.print_min = Inf)

#### Inputs ####
# Input parameters
resultTag <- c("Generic-Resistance0.1.tif",
               "Minimize-Shortfall-Species-All0.1.tif",
               "PA2010_historicLULC_baseline")

"Outputs/PrioritizationSolutions/MapsBudget0.1/Ecoprofile-Trophic0.1.tif",
"Outputs/PrioritizationSolutions/MapsBudget0.1/CAZ_Density_0.1.tif",
"Outputs/PrioritizationSolutions/MapsBudget0.1/Sum-Species-Density0.1.tif",
"Outputs/PrioritizationSolutions/MapsBudget0.1/CAZ_All_0.1.tif",
"Outputs/PrioritizationSolutions/MapsBudget0.1/Sum-Species-All0.1.tif",
resultScenarioNumber <- c(393, 394, 397)

numScenarios <- length(resultTag)
numIterations <- 3
timeInterval <- 25
timesteps <- seq(from=2010, to=2110, by=timeInterval)
numTimesteps = length(timesteps)
speciesCodes <- c("BLBR", "MAAM", "BUAM", "SCMI", "STVA", "PLCI", "RANA", "URAM", "SEAU", "DRPI", "SICA", "LEAM", "PELE", "ODVI")
numSpecies <- length(speciesCodes)


# Input files and folders
workingDir <- "Outputs/"

# Lib
mySession <- "F:/2.2.7/SyncroSim.Console.exe"
mylib <- ssimLibrary("Simulations/Monteregie_stconnect.ssim", session = mySession)
myproj <- project(mylib, "Definitions")


# Tabular Summary habitat patch per species per year ---------------------------------------------
# Create empty output table
landscapeSummary <- tibble(
  layer = as.integer(),
  level = as.character(), 
  class = as.integer(),
  id = as.integer(), 
  metric = as.character(),        
  value = as.double(), 
  Scenario = as.character(), 
  SpeciesCode = as.character(), 
  Timestep = as.integer(), 
  Iteration = as.integer()
)

# Process habitat patch maps
mapName <- "stconnect_HSOutputHabitatPatch"

# Calculate the following landscape metrics for habitat:
# lsm_c_ca - total area
# lsm_c_np - number of patches 
# lsm_c_area_mn - mean patch area
# lsm_c_area_cv - cv patch area
# lsm_c_para_mn - mean perimeter to area ratio
# lsm_c_para_cv - cv perimeter to area ratio

# Start with initial timestep (2010) and just 1 iteration 
year <- 2010  
i <- 1
for(scn in 1:numScenarios){
  # Setting up st-sim library, project, scenario
  myScenario <- scenario(myproj, resultScenarioNumber[scn])
  #datasheet(myScenario)
  
  # Get habitat patch map for year, iteration for all species
  habPatchStack <- datasheetRaster(myScenario, mapName, iteration = i, timestep = year)
  # Set area outside patches to NA so there is only 1 class in landscape (halves time for landscape metrics)
  habPatchStack[habPatchStack==0] <- NA
  
  for(speciesIndex in 1:numSpecies){
    print(paste("Working on scenario", resultTag[scn], "year", year, "iteration", i, "species", speciesCodes[speciesIndex]))
    # Get map for focal species
    habPatch <- habPatchStack[[speciesIndex]]
    # Check landscape for landscape metrics
    print(check_landscape(habPatch))
    # Calculate landscape metrics
    #start_time <- Sys.time()
    landscapeSummaryTemp <- calculate_lsm(habPatch,
                                          what = c("lsm_c_ca", "lsm_c_np", 
                                                   "lsm_c_area_mn", "lsm_c_area_cv", 
                                                   "lsm_c_para_mn", "lsm_c_para_cv")) %>%
      mutate(Scenario = resultTag[scn], 
             SpeciesCode = speciesCodes[speciesIndex], 
             Timestep = year,
             Iteration = i)
    
    #end_time <- Sys.time()
    #end_time - start_time
    landscapeSummary <- landscapeSummary %>%
      rbind(landscapeSummaryTemp)
  }
}
rm("year", "i")

# Repeat for all timesteps and numIterations 
for(scn in 1:numScenarios){
  # Setting up st-sim library, project, scenario
  myScenario <- scenario(myproj, resultScenarioNumber[scn])
  #datasheet(myScenario)
  
  for(year in timesteps){
    # Get habitat patch maps for year and all other iterations
#    for(i in 1:numIterations){
      # Get habitat patch map for year, iteration 1 for all species
      habPatchStack <- datasheetRaster(myScenario, mapName, iteration = i, timestep = year)
      # Set area outside patches to NA so there is only 1 class in landscape (halves time for landscape metrics)
      habPatchStack[habPatchStack==0] <- NA
      
      for(speciesIndex in 1:numSpecies){
        print(paste("Working on scenario", resultTag[scn], "year", year, "iteration", i, "species", speciesCodes[speciesIndex]))
        # Get map for focal species
        habPatch <- habPatchStack[[speciesIndex]]
        # Check landscape for landscape metrics
        print(check_landscape(habPatch))
        # Calculate landscape metrics
        #start_time <- Sys.time()
        landscapeSummaryTemp <- calculate_lsm(habPatch,
                        what = c("lsm_c_ca", "lsm_c_np", 
                                 "lsm_c_area_mn", "lsm_c_area_cv", 
                                 "lsm_c_para_mn", "lsm_c_para_cv")) %>%
          mutate(Scenario = resultTag[scn], 
                 SpeciesCode = speciesCodes[speciesIndex], 
                 Timestep = year,
                 Iteration = i)
    
        #end_time <- Sys.time()
        #end_time - start_time
        landscapeSummary <- landscapeSummary %>%
          rbind(landscapeSummaryTemp)
        
      }
      removeTmpFiles(h=2)
    #}    
  }
}

write_csv(landscapeSummary, paste0(workingDir, "landscapeMetrics_2021-02-19.csv"))


#Fixed/free scales
metricList <- c("Total Area", "Number of Patches", "Mean Patch Area", "Mean Perimeter-Area Ratio", "CV Perimeter-Area Ratio", "CV Patch Area")
for(m in metricList){
  #"NC_NC", "NC_45", "NC_85", "BAU_NC", "BAU_45", "BAU_85", "CON_NC", "CON_45","CON_85"
  speciesMetricSummary <- landscapeSummary %>%
    mutate(Scenario = factor(Scenario, levels = resultTag),
           Scenario = fct_recode(Scenario, "Generic Resistance" = "Generic_Resistance0.1_historicLULC_baseline",
                               "Ecoprofile (Trophic)" = "Ecoprofile_Trophic0.1_historicLULC_baseline",
                               "CAZ Current Density" = "CAZ_Density_0.1_historicLULC_baseline",
                               "Sum Current Density" = "Sum_Species_Density0.1_historicLULC_baseline",
                               "CAZ All" = "CAZ_All_0.1_historicLULC_baseline",
                               "Sum All" = "Sum_Species_All0.1_historicLULC_baseline",
                               "Minimize Shortfall All" = "Minimize_Shortfall_Species_All0.1_historicLULC_baseline"),
           metric = factor(metric, levels=c("ca", "np", "area_mn", "area_cv", "para_mn", "para_cv", "area_sd", "para_sd")),
           metric = fct_recode(metric, "Total Area" = "ca",
                               "Number of Patches" = "np",
                               "Mean Patch Area" = "area_mn",
                               "CV Patch Area" = "area_cv",
                               "Mean Perimeter-Area Ratio" = "para_mn",
                               "CV Perimeter-Area Ratio" = "para_cv",
                               "SD Patch Area" = "area_sd",
                               "SD Perimeter-area Ratio" = "para_sd")) %>%
    group_by(Scenario, SpeciesCode, Timestep, metric) %>%
    summarise(value_mn = mean(value), value_sd = sd(value)) %>%
    filter(metric == m)
  
  scnColours <- brewer.pal(n = 8, name = "Dark2")
  
  p <- ggplot(speciesMetricSummary, aes(Timestep, value_mn, group=Scenario, color=Scenario)) +
    geom_errorbar(aes(ymin=value_mn-value_sd, ymax=value_mn+value_sd), width=1) +
    geom_line(size=1.2) +
    geom_point() +
    ylab(m) +
    scale_color_manual(values = scnColours) +
    scale_x_continuous(breaks=c(2010, 2035, 2060)) +
    facet_wrap(~SpeciesCode, ncol = 4, scales="free") +
    theme_gdocs()+
    theme(axis.text.y = element_text( size=9),
          axis.text.x = element_text( size=8),
          legend.text = element_text( size=8))#,
          #legend.position = "none")
  x11(); plot(p)
  # p<-ggdraw(p) + draw_image("D:/Apex/Projects/A224/species_icons/BLBR.png", x = 0.65, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/MAAM.png", x = 0.81, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/PLCI.png", x = 0.975, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/RANA.png", x = 1.13, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/URAM.png", x = 1.3, y = 1.45, hjust = 1, vjust = 1, scale=0.05)
  
#  ggsave(paste0(workingDir, "LandscapeMetrics_,", m, "_yFree.png"), p, height=6, width=6)
  
}

