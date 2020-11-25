#### a238
#### From a script by Chloé Debyser, Andy Kleinhesselink, Bronwyn Rayfield

########################################################################
# For each Stratum in a selected StratumMask, this code:               #
# 1. Computes LULC frequency for all years                             #
# 2. Produces transition targets of LULC change for all time periods   #
# 3. Produces transition multipliers                                   #
########################################################################

# Workspace ----------------------------------------------------------------------------------------------
# Packages
library(rgrass7)
library(tidyverse)
library(ggplot2)
library(scales)
library(rsyncrosim)

# Settings
Sys.setenv(TZ='GMT')
options(stringsAsFactors=FALSE)

# Input parameters
initializaNewMapset <- T
res <- 90 # Resolution for transition targets
# ecoregion <- 73

# Directories
gisBase <- "C:/Program Files/GRASS GIS 7.8"
gisDbase <- "libraries/GRASS"
spatialDataDir <- "Data/Processed"
tabularDataDir <- "Data/Processed/tabular"
resultsDir <- "outputs"
# subscenariosDir <- paste0("ER", ecoregion, "/STSim/Results/libraries/subscenarios")

# Create directories for initial conditions, transition targets, and transition multipliers
dir.create(subscenariosDir, recursive = T)
dir.create(file.path(subscenariosDir, "stsim_InitialConditionsNonSpatial"))
dir.create(file.path(subscenariosDir, "stsim_InitialConditionsNonSpatialDistribution"))
dir.create(file.path(subscenariosDir, "stsim_TransitionTarget"))
dir.create(file.path(subscenariosDir, "stsim_TransitionMultiplierValue"))

# Tabular data - Load
stateClasses <- read_csv(file=file.path(tabularDataDir, 'Raw', 'State Class.csv'))
stratumIDs <- ecoregion

# Function - Compute transition matrix

# Set up mapset
# Open DataPreProcessing
if(initializaNewMapset){
  initGRASS(gisBase = "C:/Program Files/GRASS GIS 7.8", gisDbase = "libraries/GRASS/",  
            location = "Mont", mapset = "PERMANENT", override = TRUE)
  execGRASS("g.proj", proj4 = projection(raster("Data/Processed/LULC_FocalAreaBuffer.tif")), flags="c")
  execGRASS("g.mapset", mapset="Monteregie", flags="c")
  execGRASS("r.in.gdal", input = "Data/Processed/LULC_FocalAreaBuffer.tif", 
            output = "focal", flag = c("o", "overwrite"))
  execGRASS("g.region", raster = "focal")
}

# ST-Sim Initial Conditions ----------------

# MRC raster
execGRASS('r.in.gdal', input = "Data/stsim/secondary_stratum_FocalAreaBuffer.tif", 
          output = "sec_stratum", flags = 'overwrite')

# 2010
# execGRASS('r.mapcalc', expression = "focalInt = round(focal)", flags = "overwrite")
# execGRASS('r.mapcalc', expression = "sec_stratumInt = round(sec_stratum)", flags = "overwrite")

execGRASS('r.stats', input = c("focal", "sec_stratum"), flags = c('a', 'i', 'overwrite'), 
          separator = ',', null_value = "NA",
          output = file.path(spatialDataDir, 'tabular', paste0('focal_90m.csv')))

# ST-Sim Initial Conditions Non-Spatial Distribution
# 1996
ic2010 <- read_csv("Data/Processed/tabular/focal_90m.csv", 
                   col_names = FALSE, col_types = cols(
                     X1 = col_double(),
                     X2 = col_double(),
                     X3 = col_double()
                   )
) %>% rename('FromStateClass' = X1, 'StratumID' = X2, 'TotalAmount' = X3) %>% 
  filter(!is.na(FromStateClass)) %>% 
  filter(StratumID != 0) %>% 
  mutate(Timestep = 2010)

# Transition Multipliers ----------------------------------------------------------------------
cellAreaSqKm <- (res*res)/1000000 # Area of an individual cell, in km2

# Get targets + stratums
mylib <- ssimLibrary("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim")
sceTar <- scenario(mylib, 7)
myproj <- project(mylib, "Definitions")

targets <- datasheet(sceTar, "stsim_TransitionTarget")
secondaryStratumDatasheet <- datasheet(myproj, "stsim_SecondaryStratum")

# Combine all years
stateClassTotalAmount <- bind_rows(ic2010) %>%
  select(Timestep, StratumID, FromStateClass, TotalAmount)

myDatasheet <- targets %>%
  mutate(FromStateClass = gsub("->.*","",TransitionGroupID),
         TargetAmount = ifelse(Amount < cellAreaSqKm, 0, Amount)) %>%
  dpyr::select(-Amount) %>%
  left_join(stateClassTotalAmount, by=c("Timestep", "StratumID", "FromStateClass")) %>%
  mutate(MultiplierAmount = ifelse(is.na(TotalAmount), 0, TargetAmount/TotalAmount)) %>%
  arrange(TargetAmount) %>%
  select(Timestep, StratumID, TransitionGroupID, "Amount"=MultiplierAmount)

datasheetName <- "stsim_TransitionMultiplierValue"
write_csv(myDatasheet, file.path(subscenariosDir, datasheetName, paste0(datasheetName, "_er", ecoregion, '_', res, "m.csv")))


# Export maps ------------------------------------------------------------------------
# Save Tidal CCAP simple 
# Max compression options
rastersToSave <- c("Tidal_CCAP1996_simple",
                   "Tidal_CCAP2001_simple",
                   "Tidal_CCAP2006_simple",
                   "Tidal_CCAP2010_simple",
                   "Tidal_CCAP2016_simple")

for(ras in rastersToSave){
  outfile <- file.path(resultsDir, "Spatial", paste0(ras, "_", res, "m.tif"))
  execGRASS('r.out.gdal', 
            input = ras, 
            output = outfile, 
            type = 'Byte', 
            flags = c('c', 'm', 'overwrite'),
            createopt = c("COMPRESS=DEFLATE")) #, "NUM_THREADS=4", "ZLEVEL=9", "SPARSE_OK=YES"))
}


