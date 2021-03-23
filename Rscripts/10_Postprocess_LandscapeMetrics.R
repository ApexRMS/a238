library(tidyverse)
library(rsyncrosim)
library(raster)
library(sf)
library(viridis)
library(scales)
library(landscapemetrics)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(metapop)

options(tibble.width = Inf, tibble.print_min = Inf)

#### Inputs ####
# Input parameters
resultTag <- c(
               #"Minimize-Shortfall-Species-All",
               #"PA2010_historicLULC_baseline",
              "Generic-Resistance",
              "Ecoprofile-Trophic",
               "CAZ_Density",
               "Sum-Species-Density",
               "CAZ_All",
               "Sum-Species-All")

resultScenarioNumber <- c(410, 411, 412, 413, 414, 415)
393, 394, 397, 

numScenarios <- length(resultTag)
numIterations <- 5
timeInterval <- 50
timesteps <- c(2110) #seq(from=2010, to=2110, by=timeInterval)
numTimesteps = length(timesteps)
#speciesCodes <- c("BLBR", "MAAM", "PLCI", "RANA", "URAM")
speciesCodes <- c("BLBR", "MAAM", "BUAM", "SCMI", "STVA", "PLCI", "RANA", "URAM", "SEAU", "DRPI", "SICA", "LEAM", "PELE", "ODVI")
numSpecies <- length(speciesCodes)


# Input files and folders
dataDir <- "Data"
workingDir <- "Outputs/"

# Dispersal distance
dispersal <- read_csv(file.path(dataDir, "FocalSpeciesDispersalDistance.csv"))

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
scn <- 1
# Setting up st-sim library, project, scenario
myScenario <- scenario(myproj, resultScenarioNumber[scn])
#datasheet(myScenario)

# Get habitat patch map for year, iteration for all species
habPatchStack <- datasheetRaster(myScenario, mapName, iteration = i, timestep = year)
# Set area outside patches to NA so there is only 1 class in landscape (halves time for landscape metrics)
habPatchStack[habPatchStack==0] <- NA

for(speciesIndex in 1:numSpecies){
  #    speciesIndex<-14
  print(paste("Working on scenario", resultTag[scn], "year", year, "iteration", i, "species", speciesCodes[speciesIndex]))
  # Get map for focal species
  habPatch <- habPatchStack[[speciesIndex]]
  # # Calculate metapop capacity
  # # Create patch id's
  # rc <- clump(habPatch, directions = 4)
  # clump_id <- getValues(rc)    
  # # Get xy coords
  # xy <- xyFromCell(rc,1:ncell(rc))
  # df <- data.frame(xy, clump_id, is_clump = rc[] %in% freq(rc, useNA = 'no')[,1])
  # # Take mean xy coordinate of each patch to be centroid
  # dfm <- df %>%
  #   filter(is_clump == T) %>%
  #   group_by(clump_id) %>%
  #   summarise(xm = mean(x), ym = mean(y), A=n()*8100)
  # # Distance between patches
  # d <- dist(dplyr::select(dfm, c('xm', 'ym')))
  # 
  # # Metapop capacity
  # sol<-metacapa(d, dfm$A, alpha=1/264)#dispersal$Lower[speciesIndex])
  # 
  # # Create full row
  # mc <- data.frame(layer=1, level=NA, class=NA, id=NA, metric="metacapa", value=sol$capacity)
  
  # Calculate landscape metrics
  # Check landscape for landscape metrics
  # print(check_landscape(habPatch))
  #start_time <- Sys.time()
  landscapeSummaryTemp <- calculate_lsm(habPatch,
                                        what = c("lsm_c_ca", "lsm_c_np", 
                                                 "lsm_c_area_mn", "lsm_c_area_cv", 
                                                 "lsm_c_para_mn", "lsm_c_para_cv")) %>%
    #     bind_rows(mc) %>%
    mutate(Scenario = resultTag[scn], 
           SpeciesCode = speciesCodes[speciesIndex], 
           Timestep = year,
           Iteration = i)
  
  #end_time <- Sys.time()
  #end_time - start_time
  landscapeSummary <- landscapeSummary %>%
    rbind(landscapeSummaryTemp)
}  
  
  # Copy values to all scenarios, timestep = 1, iteration = 1
# Copy values to all iterations, timestep = 1
landscapeSummaryCopy <- landscapeSummary
for(scn in 2:numScenarios){
    landscapeSummaryTemp <- landscapeSummaryCopy %>%
      mutate(Scenario = resultTag[scn])
    
    landscapeSummary <- landscapeSummary %>%
      rbind(landscapeSummaryTemp)
  }
  # Copy values to all iterations, timestep = 1
  landscapeSummaryCopy <- landscapeSummary
  for(i in 2:numIterations){
    landscapeSummaryTemp <- landscapeSummaryCopy %>%
      mutate(Iteration = i)
  
    landscapeSummary <- landscapeSummary %>%
      rbind(landscapeSummaryTemp)
  }

rm("year", "i")

# Repeat for all timesteps and numIterations 
for(scn in 1:numScenarios){
  # Setting up st-sim library, project, scenario
  myScenario <- scenario(myproj, resultScenarioNumber[scn])
  #datasheet(myScenario)
  
  for(year in timesteps){
    # Get habitat patch maps for year and all other iterations
    for(i in 1:numIterations){
      # Get habitat patch map for year, iteration 1 for all species
      habPatchStack <- datasheetRaster(myScenario, mapName, iteration = i, timestep = year)
      # Set area outside patches to NA so there is only 1 class in landscape (halves time for landscape metrics)
      habPatchStack[habPatchStack==0] <- NA
      
      for(speciesIndex in 1:numSpecies){
#        speciesIndex<-14
        print(paste("Working on scenario", resultTag[scn], "year", year, "iteration", i, "species", speciesCodes[speciesIndex]))
        # Get map for focal species
        habPatch <- habPatchStack[[speciesIndex]]
        # # Calculate metapop capacity
        # # Create patch id's
        # rc <- clump(habPatch, directions = 4)
        # clump_id <- getValues(rc)    
        # # Get xy coords
        # xy <- xyFromCell(rc,1:ncell(rc))
        # df <- data.frame(xy, clump_id, is_clump = rc[] %in% freq(rc, useNA = 'no')[,1])
        # # Take mean xy coordinate of each patch to be centroid
        # dfm <- df %>%
        #   filter(is_clump == T) %>%
        #   group_by(clump_id) %>%
        #   summarise(xm = mean(x), ym = mean(y), A=n()*8100)
        # # Distance between patches
        # d <- dist(dplyr::select(dfm, c('xm', 'ym')))
        # 
        # # Metapop capacity
        # sol<-metacapa(d, dfm$A, alpha=1/264)#dispersal$Lower[speciesIndex])
        # 
        # # Create full row
        # mc <- data.frame(layer=1, level=NA, class=NA, id=NA, metric="metacapa", value=sol$capacity)
        
        # Calculate landscape metrics
        # Check landscape for landscape metrics
        # print(check_landscape(habPatch))
        #start_time <- Sys.time()
        landscapeSummaryTemp <- calculate_lsm(habPatch,
                                              what = c("lsm_c_ca", "lsm_c_np", 
                                                       "lsm_c_area_mn", "lsm_c_area_cv", 
                                                       "lsm_c_para_mn", "lsm_c_para_cv")) %>%
#          bind_rows(mc) %>%
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
    }    
  }
}

write_csv(landscapeSummary, paste0(workingDir, "landscapeMetrics_2021-02-19.csv"))


#Fixed/free scales

#"Minimize Shortfall All" = "Minimize-Shortfall-Species-All0.1.tif",
#"Protected Areas 2010" = "PA2010_historicLULC_baseline"),

metricList <- c("Total Area", "Number of Patches", "Mean Patch Area", "Mean Perimeter-Area Ratio", "CV Perimeter-Area Ratio", "CV Patch Area")#, "Metapopulation Capacity")
for(m in metricList){
  #"NC_NC", "NC_45", "NC_85", "BAU_NC", "BAU_45", "BAU_85", "CON_NC", "CON_45","CON_85"
  speciesMetricSummary <- landscapeSummary %>%
    mutate(Scenario = factor(Scenario, levels = resultTag),
#           Scenario = fct_recode(Scenario, "Generic Resistance" = "Generic-Resistance",
#                               "Minimize Shortfall All" = "Minimize-Shortfall-Species-All0.1.tif",
#                               "Protected Areas 2010" = "PA2010_historicLULC_baseline"),
           metric = factor(metric, levels=c("ca", "np", "area_mn", "area_cv", "para_mn", "para_cv", "area_sd", "para_sd", "metacapa")),
           metric = fct_recode(metric, "Total Area" = "ca",
                               "Number of Patches" = "np",
                               "Mean Patch Area" = "area_mn",
                               "CV Patch Area" = "area_cv",
                               "Mean Perimeter-Area Ratio" = "para_mn",
                               "CV Perimeter-Area Ratio" = "para_cv",
                               "SD Patch Area" = "area_sd",
                               "SD Perimeter-area Ratio" = "para_sd",
                               "Metapopulation Capacity" = 'metacapa')) %>%
    dplyr::group_by(Scenario, SpeciesCode, Timestep, metric) %>%
    spread(Timestep, value) %>%
    mutate(pct_change = (`2110` - `2010`)/`2010`*100) %>%
    gather("Timestep", "value", 9:10 ) %>%
    dplyr::group_by(Scenario, SpeciesCode, Timestep, metric) %>%
    mutate(Timestep = as.numeric(Timestep),
           pct_change2=ifelse(Timestep==2010, 0, pct_change)) %>%
    
  
    # mutate(pct_change = (value - lag(value))/lag(value)*100,
    #        lagvalue=lag(value, order_by=c(metric, SpeciesCode, Scenario, Timestep))) %>%
    # replace_na(list(pct_change=0)) %>%

    summarise(value_mn = mean(value), 
              value_sd = sd(value), 
              pct_change_mn = mean(pct_change2), #, na.rm=T), 
              pct_change_sd = sd(pct_change2)) %>% #, na.rm=T)) %>%
    filter(metric == m)
  
    meanMetricSummary <- speciesMetricSummary %>%
      group_by(Scenario, Timestep, metric) %>%
      summarise(value_mn = mean(value_mn), 
                value_sd = sd(value_sd), 
                pct_change_mn = mean(pct_change_mn), #, na.rm=T), 
                pct_change_sd = sd(pct_change_sd)) %>% #, na.rm=T)) %>%
#      mutate(SpeciesCode = "Mean") %>%
      filter(metric == m)
    
    speciesMetricSummary <- speciesMetricSummary %>%
      bind_rows(meanMetricSummary)

    scnColours <- brewer.pal(n = 6, name = "Dark2")
  
  p <- ggplot(speciesMetricSummary, aes(Timestep, pct_change_mn, group=Scenario, color=Scenario)) +
    #geom_errorbar(aes(ymin=pct_change_mn-pct_change_sd, ymax=pct_change_mn+pct_change_sd), width=1) +
    geom_ribbon(aes(ymin=pct_change_mn-pct_change_sd, ymax=pct_change_mn+pct_change_sd, fill=Scenario), alpha=0.2, colour=NA) +
    geom_line(size=1.2) +
    geom_point() +
    ylab(m) +
    scale_color_manual(values = scnColours) +
    scale_fill_manual(values = scnColours) +
    scale_x_continuous(breaks=c(2010, 2110)) +
    facet_wrap(~SpeciesCode, ncol = 4, scales="free") +
    theme_gdocs() +
    theme(axis.text.y = element_text( size=9),
          axis.text.x = element_text( size=8),
          legend.text = element_text( size=8),
          strip.text.x = element_text(size=0),
          legend.position = "none")
  x11(); plot(p)
  # p<-ggdraw(p) + draw_image("D:/Apex/Projects/A224/species_icons/BLBR.png", x = 0.65, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/MAAM.png", x = 0.81, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/PLCI.png", x = 0.975, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/RANA.png", x = 1.13, y = 1.45, hjust = 1, vjust = 1, scale=0.05) +
  #   draw_image("D:/Apex/Projects/A224/species_icons/URAM.png", x = 1.3, y = 1.45, hjust = 1, vjust = 1, scale=0.05)
  
#  ggsave(paste0(workingDir, "LandscapeMetrics_,", m, "_yFree.png"), p, height=6, width=6)

  p <- ggplot(speciesMetricSummary, aes(Timestep, value_mn, group=Scenario, color=Scenario)) +
    #geom_errorbar(aes(ymin=value_mn-value_sd, ymax=value_mn+value_sd), width=1) +
    geom_ribbon(aes(ymin=value_mn-value_sd, ymax=value_mn+value_sd, fill=Scenario), alpha=0.2, colour=NA) +
    geom_line(size=1.2) +
    geom_point() +
    ylab(m) +
    scale_color_manual(values = scnColours) +
    scale_fill_manual(values = scnColours) +
    scale_x_continuous(breaks=c(2010, 2110)) +
    facet_wrap(~SpeciesCode, ncol = 4, scales="free") +
    theme_gdocs() +
    theme(axis.text.y = element_text( size=9),
          axis.text.x = element_text( size=8),
          legend.text = element_text( size=8),
          strip.text.x = element_text(size=8))#,
  #legend.position = "none")
  x11(); plot(p)
  
}


landscapeSummary1<- landscapeSummary %>%
  group_by(Scenario, SpeciesCode, Iteration, metric) %>%
  mutate(pct_change = (value - lag(value))/lag(value)*100,
         lagvalue=lag(value))




"Ecoprofile (Trophic)" = "Ecoprofile_Trophic0.1_historicLULC_baseline",
"CAZ Current Density" = "CAZ_Density_0.1_historicLULC_baseline",
"Sum Current Density" = "Sum_Species_Density0.1_historicLULC_baseline",
"CAZ All" = "CAZ_All_0.1_historicLULC_baseline",
"Sum All" = "Sum_Species_All0.1_historicLULC_baseline",
