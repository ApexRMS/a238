## A238 Add species to bring the count to 14

# Packages
library(tidyverse)
library(rsyncrosim)
library(stringr)

# Lib
mylib <- ssimLibrary("libraries/BTSL_stconnect.ssim/BTSL_stconnect.ssim")
myproj <- project(mylib, "Definitions")

# Get species -------------------------------------------------------------

species <- datasheet(myproj, "stconnect_Species")

# Get State Class ---------------------------------------------------------

stateClasses <- datasheet(myproj, "stsim_StateClass") %>% 
  select(StateClassID=Name, CodeSC=ID)

# Suitability -------------------------------------------------------------

sceSuit <- scenario(mylib, 20)

suitabilityData <- read_csv("Data/Raw/Focal Species/Habitat Suitability.csv") %>% 
  left_join(species, by=c("SpeciesID"="Code")) %>% 
  select(Code=SpeciesID, CodeSC=Code, SpeciesID=Name, OldStateClassID=StateClassID) %>% 
  left_join(stateClasses, by="CodeSC") %>% 
  select(-c(OldStateClassID, Code, CodeSC))

newSceSuit <- scenario(mylib, "Habitat Suitability: Basic Landcover and Forest Age - 14 species")
saveDatasheet(ssimObject = newSceSuit, data = suitabilityData, name = "stconnect_HSHabitatSuitability") 
# write_csv(suitabilityData, "config/stsim/HabitatSuitability14Species.csv")

# Patch -------------------------------------------------------------------

patchData <- read_csv("Data/Raw/Focal Species/Habitat Patch.csv") %>% 
  rename(Code=SpeciesID) %>% 
  left_join(species, by = "Code") %>% 
  select(-Code, SpeciesID=Name)
newScePatch <- scenario(mylib, "Habitat Patch - 14 species")
saveDatasheet(ssimObject = newScePatch, data = patchData, name = "stconnect_HSHabitatPatch")

# Resistance --------------------------------------------------------------

resData <- read_csv("Data/Raw/Focal Species/Resistance.csv") %>%
  select(-c(`Class Name`)) %>% 
  left_join(species, by=c("SpeciesID"="Code")) %>% 
  select(-c(SpeciesID)) %>% 
  rename(SpeciesID=Name, CodeSC=StateClassID ) %>% 
  left_join(stateClasses, by=c("CodeSC")) %>% 
  select(-CodeSC)
newSceRes <- scenario(mylib, "Resistance - 14 species")
saveDatasheet(ssimObject = newSceRes, data = resData, name = "stconnect_HSResistance")
