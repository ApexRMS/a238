#-------------------------------------------------------------------------------
## 5.1 Gnerate points raster for CS
## 2020
## Inputs: reclassed resistance (names)
## Outputs: INI Files
#-------------------------------------------------------------------------------

library(raster)
library(sf)
library(tidyverse)

set.seed(7777)

FocalAreaWithBuffer <- raster("Data/Processed/LULC_FocalAreaBuffer.tif")
FocalAreaBuffer <- st_read("Data/Processed/buffer.shp")

the_points <- st_as_sf(st_cast(st_sample(st_cast(FocalAreaBuffer, "LINESTRING"), 
                                10000, "regular"), "POINT"))
the_points$ID <- 1:nrow(the_points)

the_distances <- st_distance(the_points)
rownames(the_distances) <- as.character(1:nrow(the_points))
colnames(the_distances) <- as.character(1:nrow(the_points))
                      
the_distances_tbl <- as.data.frame(as.table(the_distances)) %>% 
  rename(P1 = Var1, P2 = Var2, dist = Freq) %>% 
  mutate(P1 = as.numeric(P1), P2 = as.numeric(P2), dist = as.numeric(dist)) %>% 
  filter(dist != 0) %>% 
  mutate(to_keep = dist > 3650 & dist > 0)

to_remove <- the_distances_tbl %>%
  filter(to_keep == FALSE) %>% 
  group_by(dist) %>%
  filter(row_number()==1) %>% 
  pull(P1)

the_points_to_keep <- the_points %>% 
  filter(!(ID %in% to_remove))

dim(the_points_to_keep)

#plot(the_points$x)
#plot(the_points_to_keep$x, col = 2, add = TRUE)
#plot(the_points_to_keep$x, col = 2, add = FALSE)

base_raster <- FocalAreaWithBuffer
values(base_raster) <- 1
the_points_to_keep_rast <- mask(base_raster, the_points_to_keep)

the_points_to_keep_rast_clp <- clump(the_points_to_keep_rast)
the_points_to_keep_rast_clp[!is.na(the_points_to_keep_rast_clp)] <- 
  sample(rep(1:50, each=2))

writeRaster(the_points_to_keep_rast_clp, "Data/Processed/circuitscape_pts.tif", 
            overwrite = TRUE)
