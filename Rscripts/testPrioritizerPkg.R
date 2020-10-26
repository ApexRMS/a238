#####################################################################
# a238 Test Prioritizer library  for 5 test species
#
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
library(Rsymphony)
  #CITE: Hanson JO, Schuster R, Morrell N, Strimas-Mackey M, Watts ME, Arcese P, Bennett J, Possingham HP #(2020). prioritizr: Systematic Conservation Prioritization in R. R package version 5.0.2. Available #at https://CRAN.R-project.org/package=prioritizr.

## Directories
#
projectDir <- "c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity"
dataDir <- file.path(projectDir, "Data/Raw")
outDir <- file.path(projectDir, "Data/Processed")

## Functions
rescaleR <- function(x, new.min = 0, new.max = 1) {
   x.min = min(x, na.rm=TRUE)
   x.max = max(x, na.rm=TRUE)
   new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

evaltext <- function(x, y){
				sapply(x, 
				FUN = function(X){
						eval(parse(text = paste0(X, y)))}
						)}

## Load files 
  # Focal area
naturalAreasFocal <- raster(file.path(outDir, "LULCnatural_FocalArea.tif"))
naturalAreasBinaryFocal <- raster(file.path(outDir, "LULCbinary_FocalArea.tif"))

  # Protected areas
protectedAreas <- raster(file.path(outDir, "protectedAreasNatural_FocalArea.tif"))
protectedAreasBinary <-  reclassify(protectedAreas, rcl=matrix(c(0, 1, 1, 0), ncol=2, byrow=T))
#values so 0 = protected, 1 = unprotected, reflecting cost of inclusion in network

## Load focal species files and rescale ---------------------------------------------------------

speciesID <- read.csv(
					file.path(
					paste0(dataDir, "/Focal Species"), 
					"Species.csv")
					stringsAsFactors = FALSE)

specieslist <- speciesID$Code

## Using for loop
for(i in specieslist){

species <- i	

  # Load Focal species data 
resistance <- raster(
				file.path(outDir, 
				paste0(species, "_Resistance_FocalArea.tif"))) %>%
				calc(., fun = function(x){log(x + 1)}) %>%
				calc(., fun = rescaleR) %>%
				calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) %>%
				calc(., fun = function(x){1 - x}) %>%
				crop(., naturalAreasBinaryFocal)
nam1 <- paste0(species, "_Resistance")
assign(nam1, resistance)		

  # Range scale features to have suitabilities 0-1
habitatSuitability <- raster(
						file.path(outDir, 
						paste0(species, "_HabitatSuitability_FocalArea.tif"))) %>%
  				calc(., fun = function(x){sqrt(x)}) %>%
				calc(., fun = rescaleR) %>%
  				calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) %>%
				crop(., naturalAreasBinaryFocal)
nam2 <- paste0(species, "_habitatSuitability")
assign(nam2, habitatSuitability)

habitatArea <- raster(
				file.path(outDir, 
					paste0(species, "_HabitatArea_FocalArea.tif"))) %>%
				calc(., fun = function(x){log(x+1)}) %>%
  				calc(., fun = rescaleR) %>%
  				calc(., fun = function(x){ifelse(x <= 1e-6, 0, x)}) %>%
				crop(., naturalAreasBinaryFocal)
nam3 <- paste0(species, "_habitatArea")
assign(nam3, habitatArea)

	} #Ignore warnings


## Combine for all features
allSuitabilities <- c(evaltext(specieslist, "_habitatSuitability"), 
					evaltext(specieslist, "_habitatArea"),
					evaltext(specieslist, "_Resistance"))
allSuitabilities <- stack(allSuitabilities)

names(allSuitabilities) <- c(paste0(specieslist, "_habitatSuitability"), 
							paste0(specieslist, "_habitatArea"), 
							paste0(specieslist, "_Resistance"))

  # Prelim plotting
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "YlGnBu")
spplot(allSuitabilities, col.regions= pal, cuts = 8, par.settings = list(fontsize = list(text = 6)))


## Run prioritizr ------------------------------------------------------

propRep <- 0.25

## Min set solutions
  # 1 feature
p2r <- problem(naturalAreasBinaryFocal, allSuitabilities[[1:14]]) %>% #input is the cost surface + features 
      add_min_set_objective() %>% #minimize cost surface
      add_relative_targets(propRep) %>% 
      add_binary_decisions() #inclusion vs no-inclusion	
s2r <- solve(p2r)
feature_representation(p2r, s2r)
freq(s2r)


  # 2 types of features
p2r <- problem(naturalAreasBinaryFocal, allSuitabilities[[1:28]]) %>% #input is the cost surface + features 
      add_min_set_objective() %>% #minimize cost surface
      add_relative_targets(propRep) %>% #Aichi
      add_binary_decisions() #inclusion vs no-inclusion	
s2r <- solve(p2r)
feature_representation(p2r, s2r)
feature_abundances(p2r, na.rm=TRUE)
plot(s2r)

  # 3 types of features
p2r2 <- problem(naturalAreasBinaryFocal, allSuitabilities) %>% #input is the cost surface + features 
      add_min_set_objective() %>% #minimize cost surface
      add_relative_targets(0.04) %>% 
      add_binary_decisions() #inclusion vs no-inclusion	
s2r2 <- solve(p2r2)
feature_representation(p2r2, s2r2)
feature_abundances(p2r2, na.rm=TRUE)
plot(s2r2)
fr <- as.data.frame(freq(s2r2))
fr$count[2]/fr$count[1]

## with protected area constraints

b2 <- problem(naturalAreasBinaryFocal, allSuitabilities) %>% #input is the cost surface + features 
      add_min_set_objective() %>% #minimize cost surface
      add_relative_targets(0.1) %>% #Aichi
      add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
      add_binary_decisions() #inclusion vs no-inclusion	
sb2 <- solve(b2)
feature_representation(b2, sb2)
plot(sb2)
freq(sb2)
# Number of protected sites = 312091/138196 
fr <- as.data.frame(freq(sb2))
fr$count[2]/fr$count[1]

## maximize representation
numCells <- round(450287 * 0.13)

p2rMax <- problem(naturalAreasBinaryFocal, allSuitabilities[[1:14]]) %>% #input is the cost surface + features 
      add_max_features_objective(numCells) %>% #budget to max feature rep to
      add_relative_targets(0.1) %>%
	  add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
	  add_binary_decisions() #inclusion vs no-inclusion	
#presolve_check(p2rMax)
s2rMax <- solve(p2rMax)
plot(s2rMax)
feature_representation(p2rMax, s2rMax)
#does not run successfully
fr <- as.data.frame(freq(s2rMax))
fr$count[2]/fr$count[1]
#numbers don't match


##
numCells <- round(450287 * 0.13)

p2rMax <- problem(naturalAreasBinaryFocal, allSuitabilities) %>% #input is the cost surface + features 
      add_max_features_objective(numCells) %>% #budget to max feature rep to
      add_relative_targets(0.1) %>%
	  add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
	  add_binary_decisions() #inclusion vs no-inclusion	
#presolve_check(p2rMax)
s2rMax <- solve(p2rMax)
plot(s2rMax)
feature_representation(p2rMax, s2rMax)
#does not run successfully
fr <- as.data.frame(freq(s2rMax))
fr$count[2]/fr$count[1]
#numbers don't match




#utility
p2rMax1 <- problem(naturalAreasBinaryFocal, allSuitabilities) %>% #input is the cost surface + features 
      add_max_utility_objective(10^6) %>% #budget to max feature rep to
      add_binary_decisions() #inclusion vs no-inclusion	
presolve_check(p2rMax1)
s2rMax1 <- solve(p2rMax1)
feature_representation(p2rMax1, s2rMax1)
plot(s2rMax1)
print(s2rMax1)
number_of_planning_units(p2rMax1)
number_of_features(p2rMax1)
freq(s2rMax1)




## with boundary constraints

p4 <- problem(naturalAreasBinaryFocal, allSuitabilities[[1:28]]) %>% #input is the cost surface + features 
      add_min_set_objective() %>% #minimize cost surface
      add_relative_targets(propRep) %>% #Aichi
      add_locked_in_constraints(protectedAreasBinary) %>% #force inclusion of protected areas
      add_boundary_penalties(penalty = 0.0005, edge_factor = 1) %>%
      add_binary_decisions() #inclusion vs no-inclusion	
s4 <- solve(p4)
#  Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105







writeRaster(PMulti2, file.path(paste0(projectDir, "/Data/Results"), "example_prioritizRsoln.tif"), overwrite=TRUE)
writeRaster(protectedAreasBinaryFocal, file.path(outDir, "protectedAreasBinaryFocal.tif"), overwrite=TRUE)


## unneeded
  # coverage - inappropriate for our questions
#p2rMax <- problem(naturalAreasBinaryFocal, allSuitabilities) %>% #input is the cost surface + features 
#      add_max_cover_objective(10^10) %>% #budget to max feature rep to
#      add_binary_decisions() #inclusion vs no-inclusion	
#presolve_check(p2rMax)
#s2rMax <- solve(p2rMax)
#print(s2rMax)


