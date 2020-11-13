######################################################################
# a238 Multispecies connectivity           
#  Generate generic species resistance map      
#  11/2020                                  
#                                                                   
# Inputs:                                                           
#    - The combined LULC layer                                                                           
#    - The LULC/resistance crosswalk                                
#
#  Outputs:                                                          
#    - Generic Resistance raster layer                                            
#                                                                  
# Script created by  C Tucker    for ApexRMS 
#####################################################################


## Workspace ---------------------------------------------------------
  # Packages
library(tidyverse)
library(raster)
library(sf)

  # Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

setwd("~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

  # Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"

## Read in data
  # Combined LULC layers
LULC <- raster(file.path(procDataDir, "LULC_FocalArea.tif"))

 # Tabular data
   # Resistance crosswalk
crosswalk <- read_csv(
				file.path(rawDataDir, 
                "GenericResistanceCrosswalk.csv"))
                

## Create resistance layer  ---------------------------------------------------------
  # Reclassify
resistance <- LULC %>%
  reclassify(., rcl=crosswalk[, c("StateClassID", "Resistance")]) %>%
  calc(., fun=function(x){ifelse(x<0, NA, x)})

## Save outputs ---------------------------------------------------------

#geotif
writeRaster(resistance, 
			file.path(procDataDir, "Generic_Resistance_FocalArea.tif"), 
			overwrite=TRUE)

## End script-------------------------------------------------------
