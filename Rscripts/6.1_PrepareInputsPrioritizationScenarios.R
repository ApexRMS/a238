#########################################################################################
# a238 Multispecies prioritization for the Monteregie       
# Prepare data layers for input into prioritization solvers 
# 10-2020                                       					
# 
#	  Inputs (for focal species):                                   
#    -habitat suitability, habitat area, current density layers
#	   - generic resistance layer
#    - protected areas layer
#	   - natural areas layer 
#	   - ecoregions layer
#
#   - Suitability and Area data are also generated with a small amount of random variation 
#       - Names include "R" to identify these. 
#       - Generated w rnorm(x, 0, 0.001) - less than expected d/t meaningful variation
#
#   Outputs: Prepared data for scenarios to be used in Rscripts 7.1 & 7.2
#   - Prepared data separate by type, sp, and ecoregion to be using in Zonation scenarios 
#
#    * ignore message "Discarded datum Unknown based on GRS80 ellipsoid in CRS 
#       definition, but +towgs84= values preserved".
#       This is due to sp/sf compatibility bug but doesn't affect output
#     * Protected areas are omitted from input maps                                                                    
#
# Script by C Tucker for ApexRMS 									
#############################################################################################

## Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(raster)
library(sp)
library(vegan)
library(zonator)

options(warn = -1)

setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

## Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Data/Processed/PrioritizrInputs"

## Functions 
rescaleR <- function(x, new.min = 0, new.max = 1) {
  x.min = suppressWarnings(min(x, na.rm=TRUE))
  x.max = suppressWarnings(max(x, na.rm=TRUE))
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

evaltext <- function(x, y){
  sapply(x, 
         FUN = function(X){
           eval(parse(text = paste0(X, y)))}
  )}

evaltext2 <- function(x, y){
  sapply(y, 
         FUN = function(Y){
           eval(parse(text = paste0(x, Y)))}
  )}


##Load files ---------------------------------------------------------

# Load species list
speciesID <- read.csv(
  file.path(
    paste0(rawDataDir, "/Focal Species"), 
    "Species.csv"), 
  stringsAsFactors = FALSE)
specieslist <- speciesID$Code

# Focal area
LULC <- raster(file.path(procDataDir, "LULC_FocalArea.tif")) # 1465929 cells
naturalAreasFocal <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(procDataDir, "LULCbinary_FocalArea.tif"))

# Protected areas
protectedAreas <- raster(file.path(procDataDir, "protectedAreasTerrestrial_FocalArea.tif")) %>%
                  crop(., naturalAreasFocal)

# Ecoregions 
ecoregions <- raster(file.path(rawDataDir, "StudyArea/PrimaryStratum.tif")) %>%
  calc(., fun=function(x){ifelse(x==-9999, NA, x)}) %>%
  crop(., naturalAreasFocal) %>%
  mask(., naturalAreasFocal)
# Ecoregions - zone1 Adirondacks, zone 3 = StL lowlands, Zone 4 = appalachians

# Set up 3 ecoregion zones (1, 3 and 4)
zone1 <- calc(ecoregions, fun=function(x){ifelse(x==1, 1, NA)})
zone3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
zone4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})

naturalAreasBinaryFocal1 <-  mask(naturalAreasBinaryFocal, zone1)
naturalAreasBinaryFocal3 <-  mask(naturalAreasBinaryFocal, zone3)
naturalAreasBinaryFocal4 <-  mask(naturalAreasBinaryFocal, zone4)

protectedAreasNA1 <-  mask(protectedAreas, zone1)
protectedAreasNA3 <-  mask(protectedAreas, zone3)
protectedAreasNA4 <-  mask(protectedAreas, zone4)


## Load and generate feature files ---------------------------------------------------------

  # Using for-loop over ecoregion zones, and over i species
  # all layers normalized
  # log scaling for current density values only (many orders of magnitude var)
  # Range scaling reqd for input into prioritizr

ecozones <- c(1, 3, 4)

for(j in ecozones){ # multispecies inputs
  
  # ID  inputs for the correct zone
  zone <- eval(parse(text=paste0("zone", j))) 
  naturalAreasZ <- eval(parse(text=paste0("naturalAreasBinaryFocal", j))) 
  protectedAreasZ <- eval(parse(text=paste0("protectedAreasNA", j))) 
  
## Prepare species-level inputs-----------------------------------------------------------------
  for(i in specieslist){ 
    species <- i
    
    density <- raster(file.path(procDataDir, "CircuitscapeOutputs", paste0(species, "_Resistance_FocalAreaBuffer_out_cum_curmap.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
      calc(., fun = log) %>% #density has log normal distbn
      scale(.) %>%
      calc(., fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
    nam2 <- paste0(species, "_density")
    assign(nam2, density)		
    
    habitatSuitability <- raster(file.path(procDataDir, paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
      crop(., naturalAreasZ) %>%
      mask(., naturalAreasZ) %>%
      mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
      scale(.)
    habitatSuitabilityR <- habitatSuitability
    values(habitatSuitabilityR) <- values(habitatSuitability) + rnorm(length(values(habitatSuitability)), 0, 0.001) # add random var to avoid ties
    habitatSuitabilityR <- calc(habitatSuitabilityR, fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
    nam3 <- paste0(species, "_habitatSuitabilityR")
    assign(nam3, habitatSuitabilityR)
    
    habitatSuitability <- calc(habitatSuitability, fun = rescaleR) %>%
      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
    nam3 <- paste0(species, "_habitatSuitability")
    assign(nam3, habitatSuitability)
    
    if(i=="MAAM"&&j=="1"){
      
      # Replace MAAM Area raster due to single appropriate patch in zone1
      MAAM_habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
        crop(., naturalAreasZ) %>%
        mask(., naturalAreasZ) %>%
        mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
        calc(., fun=function(x){ifelse(x==271, 1, NA)})
      
      MAAM_habitatAreaR <- MAAM_habitatArea
      
    }else{
      
      habitatArea <- raster(file.path(procDataDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
        crop(., naturalAreasZ) %>%
        mask(., naturalAreasZ) %>%
        mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
        scale(.) 

      habitatAreaR <- habitatArea
      values(habitatAreaR) <- values(habitatArea) + rnorm(length(values(habitatArea)), 0, 0.001) # add random var to avoid ties
      habitatAreaR <- calc(habitatAreaR, fun = rescaleR) %>%
                      calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
      nam4R <- paste0(species, "_habitatAreaR")
      assign(nam4R, habitatAreaR) 
      
      habitatArea <- calc(habitatArea, fun = rescaleR) %>%
                    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
      nam4 <- paste0(species, "_habitatArea")
      assign(nam4, habitatArea) 
    }
    
    rm(density, habitatSuitability, habitatArea, habitatSuitabilityR, habitatAreaR)
    
  } # Finish species data loop
  

## Prepare generic resistance layer-----------------------------------------------------------------
genericRes <- raster(file.path(procDataDir, "CircuitscapeOutputs", paste0("Generic_Resistance_FocalAreaBuffer_out_cum_curmap.tif"))) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
namn3 <- paste0("genericRes", j)
assign(namn3, genericRes)
  
## Prepare summarized species resistance layers ----------------------------------------------------------
SumResDensity <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Sum_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("SumResDensity", j)
assign(nammm1, SumResDensity)
  
MeanResDensity <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Mean_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm2 <- paste0("MeanResDensity", j)
assign(nammm2, MeanResDensity)

  
## Prepare summarized ecoprofile resistances -----------------------------------------------------

  # Scenario 1
ecoprofileBird <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Bird_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileBird", j)
assign(nammm1, ecoprofileBird)
  
ecoprofileMammal <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Mammal_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm2 <- paste0("ecoprofileMammal", j)
assign(nammm2, ecoprofileMammal)
  
ecoprofileAmphibian <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Amphibian_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm3 <- paste0("ecoprofileAmphibian", j)
assign(nammm3, ecoprofileAmphibian)
  
  # Scenario2
ecoprofileOmnivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Omnivores_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileOmnivore", j)
assign(nammm1, ecoprofileOmnivore)
  
ecoprofileInsectivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Insectivores_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileInsectivore", j)
assign(nammm1, ecoprofileInsectivore)
  
ecoprofileCarnivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Carnivores_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileCarnivore", j)
assign(nammm1, ecoprofileCarnivore)
  
ecoprofileHerbivore <- raster(file.path(procDataDir, "CircuitscapeOutputs", "combined_Resistance_Herbivores_FocalAreaBuffer_out_cum_curmap.tif")) %>%
    crop(., naturalAreasZ) %>%
    mask(., naturalAreasZ) %>%
    mask(., protectedAreasZ, inv=TRUE) %>% #Omit existing protected areas from analyses
    calc(., fun = log) %>% #density has log normal distbn
    scale(.) %>%
    calc(., fun = rescaleR) %>%
    calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) #for prioritizer bug
nammm1 <- paste0("ecoprofileHerbivore", j)
assign(nammm1, ecoprofileHerbivore)
  
  
## Collect intermediate data (will have 1,3,4 as suffix)-------------------------------------------------------
  
# Combine species level inputs into stacks, by zone and feature type 
Suitability <- stack(evaltext(specieslist, "_habitatSuitability"))  
names(Suitability) <- specieslist
namm1 <- paste0("Suitability", j)
assign(namm1, Suitability)
  
SuitabilityR <- stack(evaltext(specieslist, "_habitatSuitabilityR"))  
names(SuitabilityR) <- specieslist
namm1 <- paste0("SuitabilityR", j)
assign(namm1, SuitabilityR)

Area <- stack(evaltext(specieslist, "_habitatArea"))  
names(Area) <- specieslist
namm2 <- paste0("Area", j)
assign(namm2, Area)

AreaR <- stack(evaltext(specieslist, "_habitatAreaR"))  
names(AreaR) <- specieslist
namm2 <- paste0("AreaR", j)
assign(namm2, AreaR)

Density <- stack(evaltext(specieslist, "_density"))  
names(Density) <- specieslist
namm3 <- paste0("Density", j)
assign(namm3, Density)


# Combine for all features into complete stack
All <- c(evaltext(specieslist, "_habitatSuitability"), 
           evaltext(specieslist, "_habitatArea"),
           evaltext(specieslist, "_density"))
All <- stack(All)
names(All) <- c(paste0(specieslist, "_habitatSuitability"), 
                  paste0(specieslist, "_habitatArea"),
                  paste0(specieslist, "_density"))
namm4 <- paste0("All", j)
assign(namm4, All)

AllR <- c(evaltext(specieslist, "_habitatSuitabilityR"), 
         evaltext(specieslist, "_habitatAreaR"),
         evaltext(specieslist, "_density"))
AllR <- stack(All)
names(All) <- c(paste0(specieslist, "_habitatSuitabilityR"), 
                paste0(specieslist, "_habitatAreaR"),
                paste0(specieslist, "_density"))
namm4 <- paste0("AllR", j)
assign(namm4, AllR)
  
} # Finish ecoregion loop


## Combine transformed rasters into single map per feature type -------------------------------------------

## For output of intermediate products 

# Generic resistance map
genericResAll <- mosaic(genericRes3, genericRes4, fun="max", na.rm=TRUE) %>%
  mosaic(., genericRes1, fun="max", na.rm=TRUE)

# Current density calculated from summarized resistance 
SumResDensity <- mosaic(SumResDensity3, SumResDensity4, fun="max", na.rm=TRUE) %>%
  mosaic(., SumResDensity1, fun="max", na.rm=TRUE)
MeanResDensity <- mosaic(MeanResDensity3, MeanResDensity4, fun="max", na.rm=TRUE) %>%
  mosaic(., MeanResDensity1, fun="max", na.rm=TRUE)
ResDensity <- stack(SumResDensity, MeanResDensity)
names(ResDensity) <- c("SumResDensity", "MeanResDensity")

# Current density calculated for ecoprofile resistances 
  # Scenario 1 taxonomy 
ecoprofileBird <- mosaic(ecoprofileBird3, ecoprofileBird4, fun="max", na.rm=TRUE) %>%
  mosaic(., ecoprofileBird1, fun="max", na.rm=TRUE)
ecoprofileMammal <- mosaic(ecoprofileMammal3, ecoprofileMammal4, fun="max", na.rm=TRUE) %>%
  mosaic(., ecoprofileMammal1, fun="max", na.rm=TRUE)
ecoprofileAmphibian <- mosaic(ecoprofileAmphibian3, ecoprofileAmphibian4, fun="max", na.rm=TRUE) %>%
  mosaic(., ecoprofileAmphibian1, fun="max", na.rm=TRUE)
ResTaxon <- stack(ecoprofileBird, ecoprofileMammal, ecoprofileAmphibian)
names(ResTaxon) <- c("ecoprofileBird", "ecoprofileMammal", "ecoprofileAmphibian")
  # Scenario 2 trophic
ecoprofileOmnivore <- mosaic(ecoprofileOmnivore3, ecoprofileOmnivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., ecoprofileOmnivore1, fun="max", na.rm=TRUE)
ecoprofileInsectivore <- mosaic(ecoprofileInsectivore3, ecoprofileInsectivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., ecoprofileInsectivore1, fun="max", na.rm=TRUE)
ecoprofileCarnivore <- mosaic(ecoprofileCarnivore3, ecoprofileCarnivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., ecoprofileCarnivore1, fun="max", na.rm=TRUE)
ecoprofileHerbivore <- mosaic(ecoprofileHerbivore3, ecoprofileHerbivore4, fun="max", na.rm=TRUE) %>%
  mosaic(., ecoprofileHerbivore1, fun="max", na.rm=TRUE)
ResTrophic <- stack(ecoprofileOmnivore, ecoprofileInsectivore, ecoprofileCarnivore, ecoprofileHerbivore)
names(ResTrophic) <- c("ecoprofileOmnivore", "ecoprofileInsectivore", "ecoprofileCarnivore", "ecoprofileHerbivore")


# Species level features (to be summarized as sum, mean)
Suitability <- mosaic(Suitability3, Suitability4, fun="max", na.rm=TRUE) %>%
  mosaic(., Suitability1, fun="max", na.rm=TRUE)
names(Suitability) <- specieslist

Area <- mosaic(Area3, Area4, fun="max", na.rm=TRUE)%>%
  mosaic(., Area1, fun="max", na.rm=TRUE)
names(Area) <- specieslist

Density <- mosaic(Density3, Density4, fun="max", na.rm=TRUE)%>%
  mosaic(., Density1, fun="max", na.rm=TRUE)
names(Density) <- specieslist

All <- mosaic(All3, All4, fun="max", na.rm=TRUE) %>%
      mosaic(., All1, fun="max", na.rm=TRUE)
names(All) <- c(paste0(specieslist, "_habitatSuitability"), 
                paste0(specieslist, "_habitatArea"),
                paste0(specieslist, "_density"))

# with random variation
SuitabilityR <- mosaic(SuitabilityR3, SuitabilityR4, fun="max", na.rm=TRUE) %>%
  mosaic(., SuitabilityR1, fun="max", na.rm=TRUE)
names(SuitabilityR) <- specieslist

AreaR <- mosaic(AreaR3, AreaR4, fun="max", na.rm=TRUE)%>%
  mosaic(., AreaR1, fun="max", na.rm=TRUE)
names(AreaR) <- specieslist

DensityR <- Density

AllR <- mosaic(AllR3, All4, fun="max", na.rm=TRUE) %>%
  mosaic(., AllR1, fun="max", na.rm=TRUE)
names(AllR) <- c(paste0(specieslist, "_habitatSuitabilityR"), 
                paste0(specieslist, "_habitatAreaR"),
                paste0(specieslist, "_densityR"))


## Output prepared files ----------------------------------------------------------

  # Output values are all scaled and standardized separately per ecoregion then merged into a complete map
  # Output rasterStacks of full Monteregie, w inputs for all scenarios with prioritizR analyses
  # Non-random versions are used for Rscripts 7.1-7.2
writeRaster(genericResAll, 
            file.path(outDir,  "genericResAll.tif"), 
            overwrite=TRUE)  
writeRaster(ResDensity, 
            file.path(outDir,  "ResDensity.tif"), 
            overwrite=TRUE) 
writeRaster(ResTaxon, 
            file.path(outDir,  "ResTaxon.tif"), 
            overwrite=TRUE) 
writeRaster(ResTrophic, 
            file.path(outDir,  "ResTrophic.tif"), 
            overwrite=TRUE) 

writeRaster(Area, 
            file.path(outDir,  "Area.tif"), 
            overwrite=TRUE) 
writeRaster(Density, 
            file.path(outDir,  "Density.tif"), 
            overwrite=TRUE) 
writeRaster(Suitability, 
            file.path(outDir,  "Suitability.tif"), 
            overwrite=TRUE) 
writeRaster(All, 
            file.path(outDir,  "All.tif"), 
            overwrite=TRUE) 

writeRaster(AreaR, 
            file.path(outDir,  "AreaR.tif"), 
            overwrite=TRUE) 
writeRaster(DensityR, 
            file.path(outDir,  "DensityR.tif"), 
            overwrite=TRUE) 
writeRaster(SuitabilityR, 
            file.path(outDir,  "SuitabilityR.tif"), 
            overwrite=TRUE) 
writeRaster(AllR, 
            file.path(outDir,  "AllR.tif"), 
            overwrite=TRUE) 


# Output the ecoregion shapefiles separately for species and input type as well
  # For Zonation inputs
writeRaster(Suitability1, 
              file.path(outDir, "prioritizR", paste(specieslist, "Suitability", "ecoregion1.tif", sep="_")),
              bylayer = TRUE,
              overwrite=TRUE)
writeRaster(Suitability3, 
            file.path(outDir, "prioritizR", paste(specieslist, "Suitability", "ecoregion3.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)            
writeRaster(Suitability4, 
            file.path(outDir, "prioritizR", paste(specieslist, "Suitability", "ecoregion4.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)   
writeRaster(SuitabilityR1, 
            file.path(outDir, "prioritizR", paste(specieslist, "SuitabilityR", "ecoregion1.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)
writeRaster(SuitabilityR3, 
            file.path(outDir, "prioritizR", paste(specieslist, "SuitabilityR", "ecoregion3.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)            
writeRaster(SuitabilityR4, 
            file.path(outDir, "prioritizR", paste(specieslist, "SuitabilityR", "ecoregion4.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)   

#
writeRaster(Density1, 
            file.path(outDir, "prioritizR", paste(specieslist, "Density", "ecoregion1.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)
writeRaster(Density3, 
            file.path(outDir, "prioritizR", paste(specieslist, "Density", "ecoregion3.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)            
writeRaster(Density4, 
            file.path(outDir, "prioritizR", paste(specieslist, "Density", "ecoregion4.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)   
#
writeRaster(Area1, 
            file.path(outDir, "prioritizR", paste(specieslist, "Area", "ecoregion1.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)
writeRaster(Area3, 
            file.path(outDir, "prioritizR", paste(specieslist, "Area", "ecoregion3.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)            
writeRaster(Area4, 
            file.path(outDir, "prioritizR", paste(specieslist, "Area", "ecoregion4.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)   

writeRaster(AreaR1, 
            file.path(outDir, "prioritizR", paste(specieslist, "AreaR", "ecoregion1.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)
writeRaster(AreaR3, 
            file.path(outDir, "prioritizR", paste(specieslist, "AreaR", "ecoregion3.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE)            
writeRaster(AreaR4, 
            file.path(outDir, "prioritizR", paste(specieslist, "AreaR", "ecoregion4.tif", sep="_")),
            bylayer = TRUE,
            overwrite=TRUE) 


## End script