## Workspace ---------------------------------------------------------

options(tibble.width = Inf) 

# Packages
library(tidyverse)
library(raster)

## Directories
rawDataDir <- paste0("Data/Raw")
procDataDir <- paste0("Data/Processed")

##Load files ---------------------------------------------------------

# Natural Areas in Monteregie with LULC codes
naturalAreasMonteregie <- raster(file.path(procDataDir, "LULCnatural_FocalArea.tif"))

# Protected areas in Monteregie 
protectedAreasMonteregie <- raster(file.path(procDataDir, "protectedAreasTerrestrial_FocalArea.tif"))

# Ecoregions in Monteregie
ecoregions <- raster(file.path(procDataDir, "ecoregions_FocalArea.tif")) %>%
  crop(., naturalAreasMonteregie)
ecoregions[ecoregions==-9999]<-NA
ecoregions[ecoregions==0]<-NA


# Ecoregions - zone1 Adirondacks, zone 3 = StL lowlands, Zone 4 = appalachians
ecoregionList <- unique(ecoregions)

# Monteregie extent
# Make a Monteregie study area
studyAreaMonteregie <- Which(ecoregions>0)

# Make a binary version of natural areas for Monteregie
naturalAreasBinaryMonteregie <- Which(naturalAreasFocal)

# Make a protected natural areas map for Monteregie 
protectedNaturalAreasMonteregie <- mask(naturalAreasBinaryMonteregie, protectedAreasMonteregie)



# Divide landscapes into 3 ecoregion zones (1, 3 and 4)
ecoregion1 <- calc(ecoregions, fun=function(x){ifelse(x==1, 1, NA)})
ecoregion3 <- calc(ecoregions, fun=function(x){ifelse(x==3, 1, NA)})
ecoregion4 <- calc(ecoregions, fun=function(x){ifelse(x==4, 1, NA)})
#
naturalAreas1 <-  mask(naturalAreasBinaryMonteregie, ecoregion1)
naturalAreas3 <-  mask(naturalAreasBinaryMonteregie, ecoregion3)
naturalAreas4 <-  mask(naturalAreasBinaryMonteregie, ecoregion4)
#
protectedAreas1 <-  mask(protectedAreasMonteregie, ecoregion1)
protectedAreas3 <-  mask(protectedAreasMonteregie, ecoregion3)
protectedAreas4 <-  mask(protectedAreasMonteregie, ecoregion4)
#
protectedNaturalAreas1 <- mask(naturalAreasBinaryMonteregie, protectedAreas1)
protectedNaturalAreas3 <- mask(naturalAreasBinaryMonteregie, protectedAreas3)
protectedNaturalAreas4 <- mask(naturalAreasBinaryMonteregie, protectedAreas4)




# Tabular summaries for Monteregie and ecoregions
monteregieSummary <- tibble(Name = "Monteregie",
                            ID = NA,
                            TotalPixels = cellStats(studyAreaMonteregie, sum),
                            NaturalAreaPixels = cellStats(naturalAreasBinaryMonteregie, sum),
                            ProtectedAreaPixels = cellStats(protectedAreasMonteregie, sum),
                            ProtectedNaturalAreaPixels = cellStats(protectedNaturalAreasMonteregie, sum))

ecoregion1Summary <- tibble(Name = "Ecoregion1",
               ID = 1,
               TotalPixels = cellStats(ecoregion1, sum),
               NaturalAreaPixels = cellStats(naturalAreas1, sum),
               ProtectedAreaPixels = cellStats(protectedAreas1, sum),
               ProtectedNaturalAreaPixels = cellStats(protectedNaturalAreas1, sum))

ecoregion3Summary <- tibble(Name = "Ecoregion3",
               ID = 3,
               TotalPixels = cellStats(ecoregion3, sum),
               NaturalAreaPixels = cellStats(naturalAreas3, sum),
               ProtectedAreaPixels = cellStats(protectedAreas3, sum),
               ProtectedNaturalAreaPixels = cellStats(protectedNaturalAreas3, sum))
  
ecoregion4Summary <- tibble(Name = "Ecoregion4",
               ID = 4,
               TotalPixels = cellStats(ecoregion4, sum),
               NaturalAreaPixels = cellStats(naturalAreas4, sum),
               ProtectedAreaPixels = cellStats(protectedAreas4, sum),
               ProtectedNaturalAreaPixels = cellStats(protectedNaturalAreas4, sum))
  
overallSummary <- bind_rows(monteregieSummary, ecoregion1Summary, ecoregion3Summary, ecoregion4Summary) %>%
  mutate(NaturalAreaPercent = NaturalAreaPixels/TotalPixels*100) %>%
  mutate(ProtectedAreaPercent = ProtectedAreaPixels/TotalPixels*100) %>%
  mutate(ProtectedNaturalAreaPercent = ProtectedNaturalAreaPixels/TotalPixels*100) %>%
  mutate(TargetPixels0.05 = 0.05*TotalPixels) %>%
  mutate(TargetPixels0.10 = 0.10*TotalPixels) %>%
  mutate(TargetPixels0.17 = 0.17*TotalPixels) %>%
  mutate(PixelsToAdd0.05 = TargetPixels0.05 - ProtectedAreaPixels) %>%
  mutate(PixelsToAdd0.10 = TargetPixels0.10 - ProtectedAreaPixels) %>%
  mutate(PixelsToAdd0.17 = TargetPixels0.17 - ProtectedAreaPixels) %>%
  mutate(TargetNaturalAreaPercent0.05 = TargetPixels0.05/NaturalAreaPixels*100) %>%
  mutate(TargetNaturalAreaPercent0.10 = TargetPixels0.10/NaturalAreaPixels*100) %>%
  mutate(TargetNaturalAreaPercent0.17 = TargetPixels0.17/NaturalAreaPixels*100) %>%
  mutate(TargetAddNaturalAreaPixels0.05 = TargetPixels0.05 - ProtectedNaturalAreaPixels) %>%
  mutate(TargetAddNaturalAreaPixels0.10 = TargetPixels0.10 - ProtectedNaturalAreaPixels) %>%
  mutate(TargetAddNaturalAreaPixels0.17 = TargetPixels0.17 - ProtectedNaturalAreaPixels) %>%
  mutate(TargetAddNaturalAreaPercent0.05 = TargetAddNaturalAreaPixels0.05/NaturalAreaPixels*100) %>%
  mutate(TargetAddNaturalAreaPercent0.10 = TargetAddNaturalAreaPixels0.10/NaturalAreaPixels*100) %>%
  mutate(TargetAddNaturalAreaPercent0.17 = TargetAddNaturalAreaPixels0.17/NaturalAreaPixels*100)
  
write_csv(overallSummary, file.path(procDataDir, "TargetAreaSummary.csv"))




# Calculate zone-specific % protected areas for deliverables (as # of cells)
zone1LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==1, 1, NA)})
zone3LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==3, 1, NA)})
zone4LULC <- calc(ecoregionsLULC, fun=function(x){ifelse(x==4, 1, NA)})
#
zonePA1 <- cellStats(protectedAreasNA1, sum, na.rm=T)#/cellStats(zone1LULC, sum, na.rm=T)
zonePA3 <- cellStats(protectedAreasNA3, sum, na.rm=T)#/cellStats(zone3LULC, sum, na.rm=T)
zonePA4 <- cellStats(protectedAreasNA4, sum, na.rm=T)#/cellStats(zone4LULC, sum, na.rm=T)







## 3) Set prioritization targets  ----------------------------------------------------------------




## Budget - Change budget as desired ##
Budget <- 0.1 

# Identify all natural areas available for prioritization and omit PA areas
costLayer1 <- naturalAreasBinaryFocal1 %>%
  mask(., protectedAreasNA1, inv=TRUE) #omit protected area cells
costLayer3 <- naturalAreasBinaryFocal3 %>%
  mask(., protectedAreasNA3, inv=TRUE) #omit protected area cells
costLayer4 <- naturalAreasBinaryFocal4 %>%
  mask(., protectedAreasNA4, inv=TRUE) #omit protected area cells

# Calculate ecoregion budgets, set by by ecoregion size (number of natural area pixels in zone)
NumSitesGoal1 <- round((Budget * cellStats(zone1LULC, sum, na.rm=T)), 0) - zonePA1
NumSitesGoal1 <- ifelse(NumSitesGoal1 <= 0, 1, NumSitesGoal1)
NumSitesGoal3 <- round((Budget * cellStats(zone3LULC, sum, na.rm=T)), 0) - zonePA3
NumSitesGoal3 <- ifelse(NumSitesGoal3 <= 0, 1, NumSitesGoal3)
NumSitesGoal4 <- round((Budget * cellStats(zone4LULC, sum, na.rm=T)), 0) - zonePA4
NumSitesGoal4 <- ifelse(NumSitesGoal4 <= 0, 1, NumSitesGoal4)