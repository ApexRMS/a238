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

if(!requireNamespace("BiocManager", quietly=TRUE))
	install.packages("BiocManager")
BiocManager::install("lpsymphony")	

library(prioritizr)
library(tidyverse)
library(raster)
library(sf)
library("lpsymphony")

#CITE: Hanson JO, Schuster R, Morrell N, Strimas-Mackey M, Watts ME, Arcese P, Bennett J, Possingham HP #(2020). prioritizr: Systematic Conservation Prioritization in R. R package version 5.0.2. Available #at https://CRAN.R-project.org/package=prioritizr.

## Directories
projectDir <- "~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
outDir <- file.path(projectDir, "Data/Processed")

## Test of prioritzr-----------------------------------------------------------------------
data(sim_pu_polygons)
head(sim_pu_polygons@data)
spplot(sim_pu_polygons, "locked_in", main = "Planning units in protected areas",xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
data(sim_features)

# plot the distribution of suitable habitat for each feature
plot(sim_features, main = paste("Feature", seq_len(nlayers(sim_features))),
     nr = 2)
p1 <- problem(sim_pu_polygons, features = sim_features,
              cost_column = "cost") %>%
      add_min_set_objective() %>%
      add_relative_targets(0.15) %>%
      add_binary_decisions() %>%
      add_default_solver(gap = 0)
      
s1 <- solve(p1)
 spplot(s1)    

## Test with 5 focal species data ---------------------------------------------------------

## Functions
rescaleR <- function(x, new.min = 0, new.max = 1) {
   x.min = min(x, na.rm=TRUE)
   x.max = max(x, na.rm=TRUE)
   new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}


## Load files and inputs ------

naturalAreasFocal <- raster(file.path(outDir, "LULCnaturalFocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(outDir, "LULCbinaryFocalArea.tif"))
 
  # Protected areas
protectedAreas <- raster(file.path(outDir, "protectedAreasFocalArea.tif"))
#values so 0 = protected, 1 = unprotected, reflecting cost of inclusion in network
  # Mask protected areas to naturalAreasFocal
protectedAreasFocal <-  mask(protectedAreas, naturalAreasFocal)  


## Load files for all species using for loop
specieslist <- c("BLBR", "MAAM", "URAM", "RANA", "PLCI")

for(i in specieslist){

species <- i		#BLBR  MAAM  URAM  RANA  PLCI

  # Load Focal species data 
  # Range scale features to have suitabilities 0-1
habitatSuitability <- raster(file.path(outDir, paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
  					  calc(., fun = rescaleR)
				nam <- paste0(species, "_habitatSuitability")
				assign(nam, habitatSuitability)

habitatPatch <- raster(file.path(outDir, paste0(species, "_HabitatPatch_FocalArea.tif"))) %>%
  				calc(., fun = rescaleR)
				nam2 <- paste0(species, "_habitatPatch")
				assign(nam2, habitatPatch)

curmap <- raster(file.path(outDir, paste0(species, "_curmap_FocalArea.tif"))) %>%
  				calc(., fun = rescaleR)
				nam3 <- paste0(species, "_curmap")
				assign(nam3, curmap)

habitatArea <- raster(file.path(outDir, paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
  				calc(., fun = rescaleR)
				nam4 <- paste0(species, "_habitatArea")
				assign(nam4, habitatArea)

}

#BLBR_habitatSuitability
#MAAM_habitatSuitability
#URAM_habitatSuitability
#RANA_habitatSuitability
#PLCI_habitatSuitability

#stack types of suitabilities
allHabitatSuitabilities <- stack(BLBR_habitatSuitability, MAAM_habitatSuitability, URAM_habitatSuitability, RANA_habitatSuitability, PLCI_habitatSuitability)
allHabitatAreas <- stack(BLBR_habitatArea, MAAM_habitatArea, URAM_habitatArea, RANA_habitatArea, PLCI_habitatArea)
allHabitatPatch <- stack(BLBR_habitatPatch, MAAM_habitatPatch, URAM_habitatPatch, RANA_habitatPatch, PLCI_habitatPatch)
allCurMap <- stack(BLBR_curmap, MAAM_curmap, URAM_curmap, RANA_curmap, PLCI_curmap)

#combine for all suitabilities
allSuitabilities <- stack(BLBR_habitatSuitability, MAAM_habitatSuitability, URAM_habitatSuitability, RANA_habitatSuitability, PLCI_habitatSuitability, BLBR_habitatArea, MAAM_habitatArea, URAM_habitatArea, RANA_habitatArea, PLCI_habitatArea, BLBR_curmap, MAAM_curmap, URAM_curmap, RANA_curmap, PLCI_curmap)
  spplot(allSuitabilities)


## Generate input files for prioritizR

  # x here is focal area raster with protected areas coded in locked in/locked out attribute data
 # single species example 
test <- problem(protectedAreasFocal, features = BLBR_habitatSuitability) %>%
			add_min_set_objective() %>%
      		add_relative_targets(0.2) %>%
      		add_binary_decisions() %>%
      		add_default_solver(gap = 0)
     		
P1 <- solve(test)

par(mfrow=c(1,3))
plot(BLBR_habitatSuitability)
plot(protectedAreasFocal)
spplot(P1)

# Multiple species problem
#for all natural areas in Monteregie
# this includes constraints that 'lock_in' protected areas
testMulti <- problem(naturalAreasBinaryFocal, features = allSuitabilities) %>%
			    add_min_set_objective() %>%
      			add_relative_targets(0.2) %>%
      			add_binary_decisions() %>%
      			add_default_solver(gap = 0)  %>%
     			add_locked_in_constraints(protectedAreasFocal)
     			
PMulti<- solve(testMulti)
spplot(PMulti)
