#####################################################################
# a238 ECCC Multi-species connectivity analysis            
# Map overlapping conservation priorities across scenarios        
#                                                                   
# Inputs (per species):                                                           
#                      
# Outputs:                                        
#                                             
#                                                                   
# Script created by Kyle T. Martins (Eco2urb) for ApexRMS                    
#####################################################################

# Packages
library(tidyverse)
library(raster)
library(sf)
library(reshape2)
library(stringr)
library(viridis)
library(magick)

#Set the working directory to the source file location
for(i in 1:2){setwd("..")}

####LOAD THE DATA####

#Specify the paths to the working directories
pathECCC=paste0("Data/")

#File organization data frame
flNms=read.csv(paste0(pathECCC,"Connectivity/FileListGroupings.csv"))

#Hydrological data#

#Load an example file from which to grab the projection and extent
exFile=paste0(pathECCC, "Connectivity/MapsBudget0.17/Ecoprofile-Taxon0.17.tif")
exFile=raster(exFile)

#Hydrological data
l2=st_read(paste0(pathECCC, "Hydrology/MonteregieRivers.shp"))

#Conservation priority data#
#Specify the three folders where the data are stored, the actual files are 
#loaded below
folders=c("MapsBudget0.1", "MapsBudget0.05", "MapsBudget0.17")

#Study region#
studyRegion=st_read(paste0(pathECCC, "StudyRegion/regioMonteregie.shp"))
studyRegion=st_transform(studyRegion, projection(exFile))
#plot(st_geometry(studyRegion))

#Protected areas
protAreas=st_read(paste0(pathECCC, 
               "ProtectedAreas/protectedAreasTerrestrial_FocalArea.shp"))

####COMBINE THE CONSERVATION PRIORITIES####

#Create a list to store the different results
percList=as.list(1:length(folders))
names(percList)=folders

#Run through a loop and add the binary conservation priorities to a reference 
#map, collapsing all the data per conservation budget (e.g. 10, 5, 17%)

#Create three lists for the different sets of files to be analysed
runs=c("Backend-All", "BackFrontend-Density", "BackFrontend-All")
runList=as.list(1:length(runs))
names(runList)=runs
runList[[1]]=as.character(flNms[flNms$run=="Backend-All",]$fileName)
runList[[2]]=as.character(flNms[flNms$run=="BackFrontend-Density",]$fileName)
runList[[3]]=as.character(flNms$fileName)

percList[[1]]=runList; percList[[2]]=runList; percList[[3]]=runList; 

#Specify the three conservation rates at which data were processed
valP=c(0.1, 0.05, 0.17)

for(w in 1:length(folders)){
  for(j in 1:length(runs)){ #loop for run types

    #Specify the file list
    flsECCC=list.files(paste0(pathECCC,"Connectivity/", folders[w]))
    flsECCC=flsECCC[!grepl("desktop", flsECCC)]
    
    #Specify the percentage being analysed
    perc=valP[w]
    
    flsECCCtemp=as.character(runList[[j]])
    flsECCCtemp=gsub("XXXX", perc, flsECCCtemp)
    
    flsECCC=flsECCC[flsECCC %in% flsECCCtemp]
    
    flBase=raster(paste0(pathECCC, "Connectivity/", folders[w],"/", flsECCC[1]))
    
    for(i in 2:length(flsECCC)){
      #plot(flBase)
      flAdd=raster(paste0(pathECCC, "Connectivity/", folders[w],"/", 
                          flsECCC[i]))
      #plot(flAdd)
      flBase=flBase+flAdd
      #plot(flBase)
      print(paste0(folders[w],flsECCC[i]))
    }
    percList[[w]][[j]]=flBase
  }
}

#plot(percList$MapsBudget0.05)
#plot(percList$MapsBudget0.1)
#plot(percList$MapsBudget0.17)

####MAPPING FUNCTIONS####
#Assemble all the pieces and map them together
#Define the map theme
theme_map <- function(...) {
  theme_minimal() +
    theme(
      #text = element_text(family = "Arial", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "white", color = NA), 
      panel.background = element_rect(fill = "white", color = NA), 
      legend.background = element_rect(fill = "white", color = NA),
      #panel.border = element_rect(colour = "grey", fill=NA, size=1),
      legend.position="bottom",
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
}

####MAP THE RESULTS####

for(i in 1:length(folders)){
  for(j in 1:length(runs)){ #loop for run types
    
spDF <- as(percList[[i]][[j]], "SpatialPixelsDataFrame")
spDF <- as.data.frame(spDF)

#which run is it?
runsT=runs[j]

#Specify the breaks in the legend
  if(runsT=="BackFrontend-All"){
  breaksT=c(0,5,10,15,20,25,29)
  } else if(runsT=="BackFrontend-Density"){
    breaksT=c(0,3,6,9,11) 
  } else {
    breaksT=0:6
  }

q=ggplot() + 
  # geom_polygon(data = btp, aes(x = long+2000, #shadow
  #                              y = lat-800,
  #                              group=group),
  #              color = "grey", size = 1, fill="grey")+
  # geom_sf(data = bt, size = 1, fill="white", color="black")+ #btsl
  geom_sf(data=l2, fill="#646464", size=0.2 )+
  geom_sf(data=studyRegion, fill=NA, size=0.5, color="black")+
  geom_sf(data=protAreas, fill="#f8e626ff", color="#c43c39", size=0.2)+
  geom_raster(data=spDF, aes(x=x,y=y, fill=spDF[,1])) + #raster layer
  coord_sf()+
  theme_map()+
  theme(legend.position = c(0.35, 0.8))+
  labs(x = NULL, 
       y = NULL)+
  scale_fill_gradientn( 
    colors = viridis_pal()(9),
    breaks=breaksT,
    name = "No. overlapping layers"
    #limits=c(0,21)
    )+
  guides(fill=
      guide_colorbar(
      direction = "horizontal",
      barheight = unit(2, units = "mm"),
      barwidth = unit(40, units = "mm"),
      draw.ulim = F,
      title.position = 'top',
      # some shifting around
      title.hjust = 0.5,
      label.hjust = 0.5))
  
#q

fileExport=paste0("Results/Maps/",
                  "OverlapMap-",
                  folders[i],"-",runs[j],
                  ".png")

ggsave(fileExport, q)

img=image_read(fileExport)
img=image_trim(img)

image_write(img, path = fileExport, format = "png")


  }
}

#End of script#
