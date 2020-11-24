#####################################################################
# a238 ECCC Multi-species connectivity analysis            
# Combine resistance maps by mean value per ecoprofile for calculation of 
#  current density maps
#                                                                   
# Inputs (per species):                                                           
#    - Resistance maps                       
# Outputs:                                                                                                                                        
#    - Resistance maps combined by ecoprofile scenarios                                             
#                                                                   
# Script created by Caroline Tucker  for ApexRMS                    
#####################################################################

## Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(raster)
library(sp)
library(vegan)


#setwd("~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity/")
setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

## Functions
evaltext <- function(x, y){
  sapply(x, 
         FUN = function(X){
           eval(parse(text = paste0(X, y)))}
  )}

rescaleR <- function(x, new.min = 0, new.max = 1){
  x.min = suppressWarnings(min(x, na.rm=TRUE))
  x.max = suppressWarnings(max(x, na.rm=TRUE))
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

## Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "outputs"

## Load species list
speciesID <- read.csv(
  file.path(
    paste0(rawDataDir, "/Focal Species"), 
    "Species.csv"), 
  stringsAsFactors = FALSE)
specieslist <- speciesID$Code
specieslist <- specieslist[specieslist != ""]


## Load resistance files for all species

for(i in specieslist){ # run for all species
  
  species <- i
  
  resistance <- raster(file.path(procDataDir, paste0(species, "_Resistance.tif")))
  nam2 <- paste0(species, "_resistance")
  assign(nam2, resistance)	
  
  # Ignore warnings - related to issue of sf/sp packages 
}  # End loop

# Combine resistance files
All <- evaltext(specieslist, "_resistance")
resistanceAll <- stack(All)  


## Ecoprofile Scenarios ----------------------------------------------------------------------

##  Scenario 1 - taxonomic grouping
birds <- c("SEAU",	"SICA",	"STVA",	"SCMI",	"DRPI")
birdStack <- subset(resistanceAll, birds)
mammals <- c("PELE",	"LEAM",	"ODVI",	"MAAM",	"URAM")
mammalStack <- subset(resistanceAll, mammals)
amphibians <- c("PLCI",	"RASY",	"BUAM")
amphibianStack <- subset(resistanceAll, amphibians)

#birds
combinedResistanceBird <- stackApply(birdStack, nlayers(birdStack), fun="mean", na.rm=TRUE) %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)})

#mammals
combinedResistanceMammal <- stackApply(mammalStack, nlayers(mammalStack), "mean", na.rm=TRUE) %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)})

#amphibians
combinedResistanceAmphibian <- stackApply(amphibianStack, nlayers(amphibianStack), "mean", na.rm=TRUE)  %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)}) 


## Scenario 2 - diet  groupings

carnivores <- c("MAAM", "STVA")
carnivoreStack <- subset(resistanceAll, carnivores)
herbivores <- c("LEAM", "ODVI")
herbivoreStack <- subset(resistanceAll, herbivores)
insectivores <- c("BLBR", "BUAM", "DRPI", "PLCI", "RASY", "SCMI", "SEAU")
insectivoreStack <- subset(resistanceAll, insectivores)
omnivores <- c("SICA", "PELE", "URAM")
omnivoreStack <- subset(resistanceAll, omnivores)

# carnivores
combinedResistanceCarnivores <- stackApply(carnivoreStack, nlayers(carnivoreStack), fun="mean", na.rm=TRUE) %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)})

# herbivores
combinedResistanceHerbivores <- stackApply(herbivoreStack, nlayers(herbivoreStack), "mean", na.rm=TRUE) %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)})

# insectivores
combinedResistanceInsectivores <- stackApply(insectivoreStack, nlayers(insectivoreStack), "mean", na.rm=TRUE)  %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)}) 

# omnivores
combinedResistanceOmnivores <- stackApply(omnivoreStack, nlayers(omnivoreStack), "mean", na.rm=TRUE)  %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)}) 



## Output summary resistance layers-----------------------------------------------------------

## Scenario 1 - taxonomy
writeRaster(combinedResistanceBird, 
            filename=file.path(procDataDir,  "combinedResistanceBird.tif"), 
            overwrite=TRUE)
writeRaster(combinedResistanceMammal, 
            filename=file.path(procDataDir,  "combinedResistanceMammal.tif"), 
            overwrite=TRUE)
writeRaster(combinedResistanceAmphibian, 
            filename=file.path(procDataDir,  "combinedResistanceAmphibian.tif"), 
            overwrite=TRUE)


## Scenario 2 - diet
writeRaster(combinedResistanceCarnivores, 
            filename=file.path(procDataDir,  "combinedResistanceCarnivores.tif"), 
            overwrite=TRUE)
writeRaster(combinedResistanceHerbivores, 
            filename=file.path(procDataDir,  "combinedResistanceHerbivores.tif"), 
            overwrite=TRUE)
writeRaster(combinedResistanceInsectivores, 
            filename=file.path(procDataDir,  "combinedResistanceInsectivores.tif"), 
            overwrite=TRUE)
writeRaster(combinedResistanceOmnivores, 
            filename=file.path(procDataDir,  "combinedResistanceOmnivores.tif"), 
            overwrite=TRUE)

# End script