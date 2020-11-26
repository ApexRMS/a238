#### a238
#### From a script by Chloé Debyser, Andy Kleinhesselink, Bronwyn Rayfield

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

# Directories
gisBase <- "C:/Program Files/GRASS GIS 7.8"
gisDbase <- "libraries/GRASS"
spatialDataDir <- "Data/Processed"
tabularDataDir <- "Data/Processed/tabular"
resultsDir <- "outputs"

# Set up mapset
initGRASS(gisBase = "C:/Program Files/GRASS GIS 7.8", gisDbase = "libraries/GRASS/",  
          location = "Mont", mapset = "PERMANENT", override = TRUE)
execGRASS("g.proj", proj4 = projection(raster("Data/Processed/LULC_FocalAreaBuffer.tif")), flags="c")
execGRASS("g.mapset", mapset="Monteregie", flags="c")
execGRASS("r.in.gdal", input = "Data/Processed/LULC_FocalAreaBuffer.tif", 
          output = "focal", flag = c("o", "overwrite"))
execGRASS("g.region", raster = "focal")

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
) %>% rename('StateClassID' = X1, 'ID' = X2, 'TotalAmount' = X3) %>% 
  filter(!is.na(StateClassID)) %>% 
  filter(ID != 0) %>% 
  mutate(Timestep = 2000)

# Transition Multipliers ----------------------------------------------------------------------
cellAreaSqKm <- (res*res)/1000000 # Area of an individual cell, in km2

# Get targets + stratums
mylib <- ssimLibrary("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim")
sceTar <- scenario(mylib, 7)
myproj <- project(mylib, "Definitions")

targets <- datasheet(sceTar, "stsim_TransitionTarget")
secondaryStratumDatasheet <- datasheet(myproj, "stsim_SecondaryStratum")
StateClassDatasheet <- datasheet(myproj, "stsim_StateClass") %>% 
  dplyr::select(Name, ID) %>% 
  rename(ID_SC = ID)

# Combine all years
stateClassTotalAmount <- bind_rows(ic2010) %>%
  left_join(StateClassDatasheet, by = c("StateClassID"="ID_SC")) %>% 
  mutate(FromStateClass = gsub(":All*","", Name)) %>% 
  rename(StateClass_Name=Name) %>% 
  left_join(secondaryStratumDatasheet, by = "ID") %>% 
  dplyr::select(Timestep, 'SecondaryStratumID'=Name, FromStateClass, TotalAmount) %>% 
  mutate(TotalAmount = TotalAmount/10000)

myDatasheet <- targets %>%
  mutate(FromStateClass = gsub("->.*","",TransitionGroupID),
         TargetAmount = ifelse(Amount < cellAreaSqKm, 0, Amount)) %>%
  dplyr::select(-Amount) %>%
  left_join(stateClassTotalAmount, by=c("Timestep", "SecondaryStratumID", "FromStateClass")) %>%
  drop_na() %>% 
  mutate(MultiplierAmount = ifelse(is.na(TotalAmount), 0, TargetAmount/TotalAmount)) %>%
  arrange(SecondaryStratumID) %>%
  dplyr::select(Timestep, SecondaryStratumID, TransitionGroupID, "Amount"=MultiplierAmount) %>% 
  mutate(Timestep=2010)

datasheetName <- "stsim_TransitionMultiplierValue"
sceMult <- scenario(mylib, 9)
currentMultipliers <- datasheet(sceMult, "stsim_TransitionMultiplierValue")
newMultipliers <- bind_rows(currentMultipliers, myDatasheet)

mysce <- scenario(myproj, "Monteregie_Targets_as_multipliers_baseline")
saveDatasheet(mysce, newMultipliers, datasheetName)
