## A238 Convert targets to multipliers
## Adapted from a script by Chloé Debyser, Andy Kleinhesselink, Bronwyn Rayfield

# Workspace ----------------------------------------------------------------------------------------------
# Packages
library(rgrass7)
library(tidyverse)
library(ggplot2)
library(scales)
library(rsyncrosim)
library(stringr)

# Settings
Sys.setenv(TZ='GMT')
options(stringsAsFactors=FALSE)

# Input parameters
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
  mutate(Timestep = 2000) # not 2010 because we need to match the targets

# Transition Multipliers ----------------------------------------------------------------------
cellAreaHectares <- (res*res)/10000 # Area of an individual cell, in hectares

# Get targets + stratums
mylib <- ssimLibrary("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim")
sceTar <- scenario(mylib, 7)
myproj <- project(mylib, "Definitions")

targets <- datasheet(sceTar, "stsim_TransitionTarget") %>%
  # Drop the fct because we dont want TST
  mutate(TransitionGroupID = fct_drop(TransitionGroupID))
secondaryStratumDatasheet <- datasheet(myproj, "stsim_SecondaryStratum")
StateClassDatasheet <- datasheet(myproj, "stsim_StateClass") %>% 
  dplyr::select(Name, ID) %>% 
  rename(ID_SC = ID)

# Combine with stsim data
stateClassTotalAmount <- bind_rows(ic2010) %>%
  # Join with state class
  left_join(StateClassDatasheet, by = c("StateClassID"="ID_SC")) %>% 
  mutate(FromStateClass = gsub(":All*","", Name)) %>% 
  rename(StateClass_Name=Name) %>% 
  
  # Then with sec stratum
  left_join(secondaryStratumDatasheet, by = "ID") %>% 
  dplyr::select(Timestep, 'SecondaryStratumID'=Name, FromStateClass, TotalAmount) %>% 
  # Important to turn the values from m2 to hectares
  mutate(TotalAmount = TotalAmount/10000) %>% ungroup() %>% 
  
  # SUM of forest types because targets are specified at the level of forest only
  mutate(ForestOrNot = str_detect(FromStateClass, "Forest:")) %>% 
  mutate(FromStateClass=ifelse(ForestOrNot, "Forest", FromStateClass)) %>% 
  select(-ForestOrNot) %>% 
  group_by(Timestep, SecondaryStratumID, FromStateClass) %>% 
  summarise(TotalAmount = sum(TotalAmount))

myDatasheet <- targets %>%
  # mutate(FromStateClass = gsub("->.*","",TransitionGroupID),
  #        TargetAmount = ifelse(Amount < cellAreaHectares, 0, Amount)) %>%
  # Here we decide to keep targets smaller than one cell
  mutate(FromStateClass = gsub("->.*","",TransitionGroupID), TargetAmount = Amount) %>%
  dplyr::select(-Amount) %>%
  
  # Join with the total amount
  left_join(stateClassTotalAmount, by=c("Timestep", "SecondaryStratumID", "FromStateClass")) %>%
  drop_na() %>% # this drops row count by half due to unmatched 1990 timestep
  
  # Calculate the multiplier
  mutate(MultiplierAmount = ifelse(is.na(TotalAmount), 0, TargetAmount/TotalAmount)) %>%
  arrange(SecondaryStratumID) %>%
  
  # Select columns and change the timestep
  dplyr::select(Timestep, SecondaryStratumID, TransitionGroupID, "Amount"=MultiplierAmount) %>% 
  mutate(Timestep=2010) %>%
  # mutate(TransitionGroupID = fct_drop(TransitionGroupID)) %>% 
  
  # complete combinations
  complete(Timestep, SecondaryStratumID, TransitionGroupID, fill = list(Amount =0))

# Combine with the other multipliers --------------------------------------

datasheetName <- "stsim_TransitionMultiplierValue"

# Baseline
sceMult <- scenario(mylib, 9)
currentMultipliers <- datasheet(sceMult, "stsim_TransitionMultiplierValue")
newMultipliers <- bind_rows(currentMultipliers, myDatasheet)

mysce <- scenario(myproj, "Multipliers: Targets + Climate Baseline")
saveDatasheet(mysce, newMultipliers, datasheetName, append = FALSE)

# 4.5
sceMult45 <- scenario(mylib, 97)
# currentMultipliers45 <- datasheet(sceMult45, "stsim_TransitionMultiplierValue") # This doesnt't work!
currentMultipliers45 <- read.csv("config/stsim/TransitionMultipliers45_new.csv")
newMultipliers45 <- bind_rows(currentMultipliers45, myDatasheet)

mysce45 <- scenario(myproj, "Multipliers: Targets + Climate RCP4.5")
saveDatasheet(mysce45, newMultipliers45[,-1], datasheetName, append = FALSE)

# 8.5
sceMult85 <- scenario(mylib, 98)
# currentMultipliers85 <- datasheet(sceMult85, "stsim_TransitionMultiplierValue") # This doesnt't work!
currentMultipliers85 <- read.csv("config/stsim/TransitionMultipliers85_new.csv")
newMultipliers85 <- bind_rows(currentMultipliers85, myDatasheet)

mysce85 <- scenario(myproj, "Multipliers: Targets + Climate RCP8.5")
saveDatasheet(mysce85, newMultipliers85[,-1], datasheetName, append = FALSE)
