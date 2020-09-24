#####################################################################
# a238 Test Prioritizer library  for 5 test species
# 09-2020                                       					
#                                                                                                        
#  Inputs (per species)
#    - HabitatSuitability, HabitatPatch, and curmap 
#    - Natural areas raster, protected areas raster
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

# Workspace ---------------------------------------------------------

## Packages
  #To initially install
#install.packages("prioritizr", repos = "https://cran.rstudio.com/")

#if(!requireNamespace("BiocManager", quietly=TRUE))
#	install.packages("BiocManager")
#BiocManager::install("lpsymphony")	
#install.packages("Rsymphony")

library(prioritizr)
library(tidyverse)
library(raster)
library(sf)
library("lpsymphony") # to do, install RSymphony
library("grDevices")

#CITE: Hanson JO, Schuster R, Morrell N, Strimas-Mackey M, Watts ME, Arcese P, Bennett J, Possingham HP #(2020). prioritizr: Systematic Conservation Prioritization in R. R package version 5.0.2. Available #at https://CRAN.R-project.org/package=prioritizr.

## Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")


## Test with 5 focal species data ---------------------------------------------------------

## Functions
rescaleR <- function(x, new.min = 0, new.max = 1) {
   x.min = min(x, na.rm=TRUE)
   x.max = max(x, na.rm=TRUE)
   new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}
evaltext <- function(x, y){
			sapply(x, FUN = 
				function(X){
					eval(parse(text = paste0(X, y)))
					})}

## Load files and inputs ------

naturalAreasFocal <- raster(file.path(outDir, "LULCnaturalFocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(outDir, "LULCbinaryFocalArea.tif"))
 
  # Protected areas
protectedAreas <- raster(file.path(outDir, "protectedAreasNaturalFocalArea.tif"))

#values so 0 = protected, 1 = unprotected, reflecting cost of inclusion in network
protectedAreasBinary <-  reclassify(protectedAreas, rcl=matrix(c(0, 1, 1, NA), ncol=2, byrow=T))

## Load files for all species using for loop
species <- read_csv(file.path(paste0(dataDir, "/Focal Species"), "Species.csv"))
species <- species[1:14, 1:3]

specieslist <- species$Code

for(i in specieslist){

species <- i	

  # Load Focal species data 
habitatPatch <- raster(file.path(outDir, paste0(species, "_HabitatPatch_FocalArea.tif")))
				nam4 <- paste0(species, "_habitatPatch")
				assign(nam4, habitatPatch)		

  # Range scale features to have suitabilities 0-1
habitatSuitability <- raster(file.path(outDir, paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
  				calc(., fun = rescaleR)
				nam <- paste0(species, "_habitatSuitability")
				assign(nam, habitatSuitability)

#curmap <- raster(file.path(outDir, paste0(species, "_curmap_FocalArea.tif"))) %>%
#				log(.) %>% #log values
# 				calc(., fun = rescaleR)
#				nam3 <- paste0(species, "_curmap")
#				assign(nam3, curmap)

habitatArea <- raster(file.path(outDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
  				calc(., fun = rescaleR)
				nam4 <- paste0(species, "_habitatArea")
				assign(nam4, habitatArea)

habitatPatchesFocal <- raster(file.path(outDir, paste0(species, "_HabitatID_FocalArea.tif")))

}
#Ignore warnings

  # Combine for all features
allSuitabilities <- c(evaltext(specieslist, "_habitatSuitability"), 
				evaltext(specieslist, "_habitatPatch"),
				evaltext(specieslist, "_habitatArea"))
allSuitabilities <- stack(allSuitabilities)

names(allSuitabilities) <- c(paste0(specieslist, "_habitatSuitability"), paste0(specieslist, "_habitatPatch"), paste0(specieslist, "_habitatArea"))

library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "YlGnBu")
spplot(allSuitabilities, col.regions= pal, cuts = 8, par.settings = list(fontsize = list(text = 8)))


## Run prioritizr

# this includes constraints that 'lock_in' protected areas
testMulti2 <- problem(naturalAreasBinaryFocal, features = allSuitabilities) %>%
			    add_min_set_objective() %>%
      			add_relative_targets(0.17) %>%
      			add_binary_decisions() %>%
      			add_default_solver(gap = 0)  %>%
     			add_locked_in_constraints(protectedAreasBinary)
     			
system.time(PMulti2 <- solve(testMulti2))
#   user  system elapsed 
# 13.472   2.470  16.234
spplot(PMulti2, cuts = 8, col.regions= pal)
plot(PMulti2)
#note that solution map doesn't include protected areas (i.e. you would have to add the two maps)
spplot(stack(PMulti2, protectedAreasBinaryFocal))

feature_representation(testMulti2, PMulti2)



writeRaster(PMulti2, file.path(paste0(projectDir, "/Data/Results"), "example_prioritizRsoln.tif"), overwrite=TRUE)
writeRaster(protectedAreasBinaryFocal, file.path(outDir, "protectedAreasBinaryFocal.tif"), overwrite=TRUE)

