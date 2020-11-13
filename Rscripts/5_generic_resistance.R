
library(raster)
library(tidyverse)
library(rgrass7)

LU <- raster("Data/Processed/LULC_FocalAreaBuffer.tif")

rules <- read_csv("Data/Raw/GenericResistanceCrosswalk.csv")

LuReclassed <- reclassify(LU, rcl = rules[,2:3])
writeRaster(LuReclassed, 
            "Data/Processed/Generic_Resistance_FocalAreaBuffer.tif", 
            overwrite = TRUE)

LUUrban <- LU == 400
LUUrban[LUUrban<1] <- NA
writeRaster(LUUrban, 
            "Data/Processed/LU_Urban.tif",
            overwrite = TRUE)

# Fill gaps
initGRASS(gisBase = "C:/Program Files/GRASS GIS 7.8", gisDbase = "libraries/GRASS/",  
          location = "Mont", mapset = "PERMANENT", override = TRUE)
execGRASS("g.proj", proj4 = projection(LuReclassed), flags="c")
execGRASS("g.mapset", mapset="Monteregie", flags="c")
execGRASS("r.in.gdal", input = "Data/Processed/LU_Urban.tif", 
          output = "urban", flag = c("o", "overwrite"))
execGRASS("g.region", raster = "urban")

# Doesnt work on windows, have to open the console
# execGRASS("r.fill.gaps", input = "LUCreclassed", output = "LUCreclassed_filled ",
#           distance = 3, mode = "mode", 
#           minimum = 1, maximum = 1)

filled <- raster("data/Processed/LU_Urban_filled.tif")
LuReclassed[filled==1] <- 1000
writeRaster(LuReclassed, 
            "Data/Processed/Generic_Resistance_FocalAreaBuffer.tif", 
            overwrite = TRUE)
