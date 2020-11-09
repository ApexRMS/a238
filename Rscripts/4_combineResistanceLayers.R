#####################################################################
# a238 ECCC Multi-species connectivity analysis            
# Combine resistance maps            
#                                                                   
# Inputs (per species):                                                           
#    - Resistance maps                       
# Outputs:                                                                                                                                        
#    - Resistance maps combined across species                                              
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

#specieslist <- c("PLCI", "RASY", "URAM", "MAAM", "BLBR")


## Load resistance files for all species

for(i in specieslist){ # run for all species
  
  species <- i
  
  resistance <- raster(file.path(procDataDir, paste0(species, "_ResistanceBuffer.tif")))
  nam2 <- paste0(species, "_resistance")
  assign(nam2, resistance)	
  
  # Ignore warnings - related to issue of sf/sp packages 
}  # End loop

# Combine resistance files
All <- evaltext(specieslist, "_resistance")
resistanceAll <- stack(All)  


## Calculate summary files

#mean
combinedResistanceMean <- stackApply(resistanceAll, nlayers(resistanceAll), "mean", na.rm=TRUE) %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)})

#max
combinedResistanceMax <- stackApply(resistanceAll, nlayers(resistanceAll), "max", na.rm=TRUE) %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)})

#sum
combinedResistanceSum <- stackApply(resistanceAll, nlayers(resistanceAll), "sum", na.rm=TRUE)  %>%
  calc(., fun=function(x){rescaleR(x, new.min = 1, new.max = 32)}) %>%
  mask(., resistanceAll)


## Output summary resistance layers-----------------------------------------------------------


writeRaster(combinedResistanceMean, 
            filename=file.path(outDir,  "combinedResistanceRaster_Mean.tif"), 
            overwrite=TRUE)
writeRaster(combinedResistanceMax, 
            filename=file.path(outDir,  "combinedResistanceRaster_Max.tif"), 
            overwrite=TRUE)
writeRaster(combinedResistanceSum, 
            filename=file.path(outDir,  "combinedResistanceRaster_Sum.tif"), 
            overwrite=TRUE)


# End script