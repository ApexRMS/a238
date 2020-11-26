#####################################################################
# a238 Multispecies prioritization for the Monteregie       
# Plotting and analysis of prioritizing output
# 11-2020                                       					
#                             
#	1.Inputs (for focal species):                                   
#    -solution rasters 
#    
#                                                                   
# Script by C Tucker for ApexRMS 									
#####################################################################

## Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(raster)
library(sp)
library(zonator)
library(corrplot)
library(purrr)
dev.new()

#setwd("~/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity/")
setwd("c:/Users/carol/Dropbox/Documents/ApexRMS/Work/A238 - Multispecies Connectivity")

## Directories
rawDataDir <- "Data/Raw"
procDataDir <- "Data/Processed"
outDir <- "Results/PrioritizationMaps"

# Load species list
speciesID <- read.csv(
  file.path(
    paste0(rawDataDir, "/Focal Species"), 
    "Species.csv"), 
  stringsAsFactors = FALSE)
specieslist <- speciesID$Code


## Load prioritization solution maps
rasters <- list.files(outDir, full.names=T, ".tif$")
output <- list()
  
  for(i in 1:length(rasters)){
    output[[i]] <- raster(rasters[i])
  }
outputAll <- stack(output)

mapNames <- c("SumAll",
              "SumSuit",
              "SumDensity",
              "SumArea",
              "MaxAll",
              "MaxSuit",
              "MaxDensity",
              "MaxArea",
              "MeanAll",
              "MeanSuit",
              "MeanDensity",
              "MeanArea",
              "SumResistDensity",
              "MeanResistDensity",
              "MaxResistDensity",
              "MSAll",
              "MSSuit",
              "MSDensity",
              "MSArea",
              "MUAll",
              "MUSuit",
              "MUDensity",
              "MUArea")
names(outputAll) <- mapNames


## Analysis of solution maps ----------------------------------------------------------------
numModels = 23
## Calculate  repnumModels esentation

# Matrix in which to save representation values, col=models, rows=features
allModelRep <- matrix(NA, nrow=(length(specieslist) * 3 + 1), ncol=numModels)
colnames(allModelRep) <- names(outputAll)[1:numModels]
rownames(allModelRep) <- c("Size", paste(rep(specieslist, each=3), c("Suit", "Density", "Area"), sep="_"))

# Per model output, extract total amount of each feature included in solution

for(m in 1:numModels){
  
  modelSol <- outputAll[[m]]
  allModelRep[1, m] <- freq(modelSol)[2, "count"]
  
  for(l in 1:length(specieslist)){
    sp <- specieslist[l]
    #Suitability
    coverage <- overlay(x=modelSol, y=Suitability[[l]], fun=function(x, y){(x * y)})
    suitcover <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(Suitability[[l]], "sum", na.rm=TRUE)
    # Density
    coverage <- overlay(x=modelSol, y=Density[[l]], fun=function(x, y){(x * y)})
    densecover <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(Density[[l]], "sum", na.rm=TRUE)
    #Area
    coverage <- overlay(x=modelSol, y=Area[[l]], fun=function(x, y){(x * y)})
    areacover <- cellStats(coverage, "sum", na.rm=TRUE)/cellStats(Area[[l]], "sum", na.rm=TRUE)
    
    allModelRep[((l*3)-1):((l*3)+ 1), m] <- c(suitcover, densecover, areacover)
  }
}

#barplot(allModelRep[1, ], cex.names=0.75, angle=30, beside=TRUE) #Confirm solutions are the same size
barplot(allModelRep[-1, ], cex.names=0.75, angle=30, beside=TRUE)

repAll <- allModelRep[-1, c("SumAll", "MaxAll", "MeanAll", "MSAll", "MUAll")]


suitAll <- seq(from=1, to=nrow(repAll), by=3)
densAll <- seq(from=2, to=nrow(repAll), by=3)
areaAll <- seq(from=3, to=nrow(repAll), by=3)

colours <- c("lightsalmon",  "lightsalmon1", "lightsalmon2", "lightsalmon3", "lightsalmon4",
             "indianred1", "indianred2", "indianred",  "indianred3", "indianred4", 
             "mediumpurple ", "mediumpurple1", "mediumpurple2",  "mediumpurple3")

# Plot  models w representation of all features
barplot((repAll[suitAll, ]), 
        beside=TRUE, 
        names.arg= c("Sum All", "Max All", "Mean All", "Min Shortfall", "Max Utility"),
        col = colours,
        ylim=c(0, 0.7), 
        xlab="Model approach", 
        ylab="% Representation", 
        main="Habitat suitability")
legend("topright", 
       legend=specieslist, 
       fill=colours,
       cex=0.75)

barplot((repAll[densAll, ]), 
        beside=TRUE, 
        names.arg= c("Sum All", "Max All", "Mean All", "Min Shortfall", "Max Utility"),
        col = colours,
        ylim=c(0, 0.7), 
        xlab="Model approach", 
        ylab="% Representation", 
        main="Current density")
legend("topright", 
       legend=specieslist, 
       fill=colours, 
       cex=0.75)

barplot((repAll[areaAll, ]), 
        beside=TRUE, 
        names.arg= c("Sum All", "Max All", "Mean All", "Min Shortfall", "Max Utility"),
        col = colours,
        ylim=c(0, 0.7), 
        xlab="Model approach", 
        ylab="% Representation", 
        main="Habitat area")
legend("topright", 
       legend=specieslist, 
       fill=colours,
       cex=0.75)


## Calculate raster correlations

set.seed(10)
#plot all model correlations
subSamp <- sampleRandom(outputAll, size= ncell(outputAll[[1]]) * 0.5)
cors <- cor(subSamp)
corrplot(cors, "ellipse", "upper", tl.cex=0.65, diag=FALSE)

#complete models only
corsAll <- cor(sampleRandom(subset(outputAll, c(1,5,9,16,20)), size= ncell(outputAll[[1]]) * 0.5 ))
corrplot.mixed(corsAll, upper="ellipse", lower="number", tl.srt = 0, lower.col="black", tl.cex=0.65, tl.pos="d", number.cex=0.75)

# density only
corsDensity <- cor(sampleRandom(subset(outputAll, c(3,7,11, 13,14,15,18,22)), size= ncell(outputAll[[1]]) * 0.5 ))

dev.new()
corrplot.mixed(corsDensity, upper="ellipse", lower="number", lower.col="black", tl.pos="d", tl.cex=0.65)



# jaccard
jacs <- cross_jaccard(outputAll, thresholds = 0.05)
jacsAll <- cross_jaccard(subset(outputAll, c(1,5,9,13,17)), thresholds = 0.05)
corrplot(as.matrix(jacsAll$`0.05`), "ellipse", type="upper", diag=FALSE, is.corr=FALSE)



##
mapNames <- c("genericResSol", 1
              "FinalEcoprofileTaxon", 2
              "FinalEcoprofileTrophic", 3
              "FinalSumResDensity", 4
              "FinalMeanResDensity", 5
              "sum_FinalSuitability",6
              "sum_FinalDensity",7
              "sum_FinalArea",8
              "sum_FinalAll",9
              "mean_FinalSuitability",10
              "mean_FinalDensity",11
              "mean_FinalArea",12
              "mean_FinalAll",13
              "minShortSuitSol",14
              "minShortDensitySol",15
              "minShortAreaSol",16
              "minShortAllSol",17
              "maxUtilitySuitSol",18
              "maxUtilityDensitySol",19
              "maxUtilityAreaSol",20
              "maxUtilityAllSol")21