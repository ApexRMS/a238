#####################################################################
# a238 ECCC Multi-species connectivity analysis            
# Conduct overlap analysis for different scenarios           
#                                                                   
# Inputs (per species):                                                           
#                      
# Outputs:                                        
#                                             
#                                                                   
# Script created by Kyle T. Martins (Eco2urb) for ApexRMS                    
#####################################################################

# Packages
library(raster)
library(reshape2)
library(stringr)
library(viridis)
library(ggplot2)
library(tidyr)

#Set the working directory to the source file location
for(i in 1:2){setwd("..")}

####FILE ORGANIZING####

#Specify the directories where the data are loaded
folders=c("MapsBudget0.1", "MapsBudget0.05", "MapsBudget0.17")

#Specify the three conservation rates at which data were processed
valP=c(0.1, 0.05, 0.17)

#Specify the three types of runs for the analysis, create a list to store data
runs=c("Backend-All", "BackFrontend-Density", "BackFrontend-All")
runList=as.list(1:length(runs))
names(runList)=runs

#Create a dataframe to organize file names, run types
dir=c("./Data/Connectivity/")
fls=list.files(paste0(dir, "MapsBudget0.17"))
fls=fls[!grepl("desktop|xml", fls)]

flNms=data.frame(fileName=gsub("0.17", "XXXX", fls))
flNms$fileAbbr=c("ABFAll", "ABFArea", "ABFDens", "ABFSuit",
                 "CAZAll", "CAZArea", "CAZDens", "CAZSuit",
                 "EcoTax", "EcoTroph", "GenRes", 
                 "MaxUtiAll", "MaxUtiArea", "MaxUtiDense","MaxUtiSuit", 
                 "MeanAll", "MeanArea", "MeanDense", "MeanRes", "MeanSuit",
                 "MinShortAll", "MinShortArea", "MinShortDens", "MinShortSuit",
                 "SumAll", "SumArea", "SumDense", "SumRes", "SumSuit")
frontEnd=c("Generic-Resistance","Ecoprofile-Taxon","Ecoprofile-Trophic",
           "Sum-Species-Resistance","Mean-Species-Resistance")
flNms$frontback=ifelse(grepl(paste0(frontEnd, collapse="|"),
                             flNms$fileName), "Frontend", "Backend")
flNms$run=ifelse(grepl("All", flNms$fileName), 
                 "Backend-All",
                 ifelse(grepl("Density", flNms$fileName)| 
                          flNms$frontback=="Frontend",
                          "BackFrontend-Density", 
                          "None"))

#Export a copy of the table to help with file organization
  #write.csv(flNms, 
  #"./Data/Connectivity/FileListGroupings.csv",
  #         row.names=FALSE)

#I updated the table to be able to order the variables better in the 
#correlograms
nmOrder=read.csv("./Data/Connectivity/FileListGroupings.csv")
head(nmOrder)
tail(nmOrder)

#Creating lists#

#Create a list for storing all the different correlation coefficients
JList=as.list(1:length(folders))
names(JList)=folders
tList=as.list(1:length(runs))
names(tList)=runs
JList[[1]]=tList; JList[[2]]=tList; JList[[3]]=tList; 

#Create a list for storing the variable ordering
OrderList=as.list(1:length(runs))
names(OrderList)=runs
orderVarNames=c("BA_Order", "BFD_Order", "BFA_Order")

for(i in 1:length(runs)){
var=orderVarNames[i]
nmOrderT=nmOrder[c("Abbreviation", var)] #need to edit here
nmOrderT=nmOrderT[!is.na(nmOrderT[var]),]
nmOrderT=as.character(nmOrderT[order(nmOrderT[var]),]$Abbreviation)
OrderList[[i]]=nmOrderT
}

####GRAPHING AND INDEX CALCULATION####

#'Comment 20210127 KTM: Would have been more efficient to just calculate the 
#'Jaccard index for all variable combinations once and then subset the resulting
#'table, but for now the following solution works well enough.  

#Processing loop to generate correl graphics summarizing jaccard index
for(w in 1:length(valP)){ #loop for conservation rates
  for(j in 1:length(runs)){ #loop for run types
  #for(j in 2){ 
  
#Take subset of df for organizing file names
if(runs[j]=="BackFrontend-All"){
flNmst=flNms
} else {
flNmst=flNms[flNms$run==runs[j],]
}

#Collate the list of files
perc=valP[w]
flNmst$fileName=gsub("XXXX", perc, flNmst$fileName)

#Reorder the files in the list to match that of the desired output
variableOrder=OrderList[[j]]
flNmst$fileAbbr=as.factor(flNmst$fileAbbr)
flNmst$fileAbbr=factor(flNmst$fileAbbr, levels=rev(variableOrder))
flNmst=flNmst[with(flNmst, order(fileAbbr)),]

#Create list to store rasters
flList=as.list(1:length(flNmst$fileName))
names(flList)=flNmst$fileName

#Load rasters in list
for(i in 1:length(flList)){
flList[[i]]=raster(paste0(dir, folders[w], "/", flNmst$fileName[i]))
}

#Generate index for every pairwise combination of files
nflsVector=1:length(flNmst$fileName)
combDf=data.frame(t(combn(nflsVector, 2)))
nComb=dim(combDf)[1]

#List for storing the jaccard index values
jList=as.list(1:nComb)

#Run through loop and calculate jaccard index for every combo of files
for(i in 1:nComb){
r1=flList[[combDf$X1[i]]]
r2=flList[[combDf$X2[i]]]
nmr1=names(flList)[[combDf$X1[i]]]
nmr2=names(flList)[[combDf$X2[i]]]
names(jList)[i]=paste(nmr1, nmr2, sep="&")
r1r2Intersect=as.numeric(summary(values(r1)==1 & values(r2)==1)[3])
r1Freq=as.numeric(freq(r1)[2,2])
r2Freq=as.numeric(freq(r2)[2,2])
jList[[i]]=round(r1r2Intersect/(r1Freq+r2Freq-r1r2Intersect),4)
print(round(i/nComb,2)*100)
}

#Format the df so can be processed with ggplot (long format), merge in the 
#file name abbreviations
JData=melt(jList)
JData=data.frame(JData, str_split_fixed(JData$L1, "&", 2))
JData=JData[names(JData)!="L1"]
names(JData)=c("jIndex", "v1", "v2")
JData=JData[c("v1", "v2", "jIndex")]
JData=merge(JData, flNmst[c("fileName", "fileAbbr")], 
            by.x="v1", by.y="fileName")
names(JData)[names(JData)=="fileAbbr"]="fileAbbrv1"
JData=merge(JData, flNmst[c("fileName", "fileAbbr")], 
            by.x="v2", by.y="fileName")
names(JData)[names(JData)=="fileAbbr"]="fileAbbrv2"
JData=JData[c(4,5,2,1,3)]

#Store the data in a list
JList[[w]][[j]]=JData

}
} #End of loops

#Loop for plotting the correlation figures
for(w in 1:length(valP)){ #loop for conservation rates
  for(j in 1:length(runs)){ #loop for run types
  #for(j in 3){ #loop for run types
    
JData=JList[[w]][[j]]
perc=valP[w]
    
#Define the text size within the plot, and the legend text size
if(runs[j]=="BackFrontend-All"){
  textSize=1.2
  legendTextSize=7
} else {
  textSize=3
  legendTextSize=11
}

#Order the variables
#Reorder the files in the list to match that of the desired output
variableOrder=OrderList[[j]]
JData$fileAbbrv1=as.factor(JData$fileAbbrv1)
JData$fileAbbrv1=factor(JData$fileAbbrv1, levels=rev(variableOrder))
JData$fileAbbrv2=as.factor(JData$fileAbbrv2)
JData$fileAbbrv2=factor(JData$fileAbbrv2, levels=variableOrder)

#Plotting the results
ggplot(JData, aes(fileAbbrv2, fileAbbrv1, fill = jIndex))+
  geom_tile(color = "white")+
  scale_fill_gradientn(name="Jaccard Index",
                     limits=c(0,1),
                     colours=c("navyblue", "darkmagenta", "darkorange1")) +
  theme_minimal()+ # minimal theme
  scale_y_discrete(position = 'left') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = legendTextSize, hjust = 1),
        axis.text.y = element_text(size = legendTextSize))+
  #scale_x_discrete(labels=lbt$fileAbbr2)+  
  #scale_y_discrete(labels=lbt$fileAbbr2)+
  coord_fixed() + 
  geom_text(aes(fileAbbrv2, fileAbbrv1, label = round(jIndex, 2)), 
            color = "white", size = textSize) +
  #dark_theme_gray()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(1.0, 0.6),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                             title.position = "top", 
                             title.hjust = 0.5))

#Export the data
ggsave(paste0("./Results/Figures/",
              "Correl-",
              runs[j],"-",
              perc,
              ".png"))

}
}

####CREATING COMBINED GRAPHS####

###BackFrontend-Density SUMMARY GRAPH####
#Combine the data into a consolidated df for the run BackFrontend-Density
#JList
JList2=melt(JList)
head(JList2)

JList3=JList2[JList2$L2=="BackFrontend-Density", c("L2", "L1", 
                                                   "fileAbbrv1", 
                                                   "fileAbbrv2", 
                                                   "value")]
JList3=spread(JList3, L1, value )
JList3$annot=paste0("(", round(JList3$MapsBudget0.05,2), ", ", 
                    round(JList3$MapsBudget0.1,2), ")")
head(JList3)

JList3$fileAbbrv1=as.factor(JList3$fileAbbrv1)
JList3$fileAbbrv1=factor(JList3$fileAbbrv1, levels=rev(nmOrderT))
JList3$fileAbbrv2=as.factor(JList3$fileAbbrv2)
JList3$fileAbbrv2=factor(JList3$fileAbbrv2, levels=nmOrderT)

ggplot(JList3, aes(fileAbbrv2, fileAbbrv1, fill = MapsBudget0.17))+
  geom_tile(color = "white")+
  scale_fill_gradientn(name="Jaccard Index",
                       limits=c(0,1),
                       colours=c("navyblue", "darkmagenta", "darkorange1")) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 11, hjust = 1),
        axis.text.y = element_text(size = 11))+
  scale_y_discrete(position = 'left') + 
  #scale_x_discrete(labels=lbt$fileAbbr2)+  
  #scale_y_discrete(labels=lbt$fileAbbr2)+
  coord_fixed() + 
  geom_text(aes(fileAbbrv2, fileAbbrv1, label = annot), 
            color = "white", size = 1.3, nudge_y=0.2) +
  geom_text(aes(fileAbbrv2, fileAbbrv1, label = round(MapsBudget0.17, 2)), 
            color = "white", size = 2.5, nudge_y=-0.1) +
  #dark_theme_gray()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(1, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", 
                               title.hjust = 0.5))

ggsave(paste0("./Results/Figures/",
              "Correl-",
              runs[2],"-",
              "Combined",
              ".png"))

###Backend-All SUMMARY GRAPH####
#Combine the data into a consolidated df for the run Backend-All

JList4=JList2[JList2$L2=="Backend-All", c("L2", "L1", 
                                          "fileAbbrv1", 
                                          "fileAbbrv2", 
                                          "value")]
JList4=spread(JList4, L1, value )
JList4$annot=paste0("(", round(JList4$MapsBudget0.05,2), ", ", 
                    round(JList4$MapsBudget0.1,2), ")")
head(JList4)

fullList=c(unique(JList4$fileAbbrv2), unique(JList4$fileAbbrv1))
fullList[!fullList %in% unique(JList4$fileAbbrv2)]
fullList[!fullList %in% unique(JList4$fileAbbrv1)]

JList4$fileAbbrv1=as.factor(JList4$fileAbbrv1)
JList4$fileAbbrv1=factor(JList4$fileAbbrv1, levels=rev(nmOrderT))
JList4$fileAbbrv2=as.factor(JList4$fileAbbrv2)
JList4$fileAbbrv2=factor(JList4$fileAbbrv2, levels=nmOrderT)

ggplot(JList4, aes(fileAbbrv2, fileAbbrv1, fill = MapsBudget0.17))+
  geom_tile(color = "white")+
  scale_fill_gradientn(name="Jaccard Index",
                       limits=c(0,1),
                       colours=c("navyblue", "darkmagenta", "darkorange1")) +
  theme_minimal()+ # minimal theme
  scale_y_discrete(position = 'left') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 11, hjust = 1),
        axis.text.y = element_text(size = 11))+
  #scale_x_discrete(labels=lbt$fileAbbr2)+  
  #scale_y_discrete(labels=lbt$fileAbbr2)+
  coord_fixed() + 
  geom_text(aes(fileAbbrv2, fileAbbrv1, label = annot), 
            color = "white", size = 2.7, nudge_y=0.22) +
  geom_text(aes(fileAbbrv2, fileAbbrv1, label = round(MapsBudget0.17, 2)), 
            color = "white", size = 3.5, nudge_y=0.0) +
  #dark_theme_gray()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position =  c(1, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", 
                               title.hjust = 0.5))

ggsave(paste0("./Results/Figures/",
              "Correl-",
              runs[1],"-",
              "Combined",
              ".png"))

###BackFrontend-All SUMMARY GRAPH####
#For reference, create a df with all of the data consolidated

JList5=JList2[JList2$L2=="BackFrontend-All", c("L2", "L1", 
                                          "fileAbbrv1", 
                                          "fileAbbrv2", 
                                          "value")]
JList5=spread(JList5, L1, value )
JList5$annot=paste0("(", round(JList5$MapsBudget0.05,2), ", ", 
                    round(JList5$MapsBudget0.1,2), ")")


####TABLE EXPORT####
head(JList3)
JList3Tidy=JList3[names(JList3)!="annot"]
names(JList3Tidy)=c("Comparison", "fileAbbrv1", "fileAbbrv2", "jIndex0.05",
                    "jIndex0.1", "jIndex0.17")
write.csv(JList3Tidy, paste0("./Results/Tables/",
                             "BackFrontend-Density-JaccardIndex-Summary.csv"),
          row.names=FALSE)

head(JList4)
JList4Tidy=JList4[names(JList4)!="annot"]
names(JList4Tidy)=c("Comparison", "fileAbbrv1", "fileAbbrv2", "jIndex0.05",
                    "jIndex0.1", "jIndex0.17")
write.csv(JList4Tidy, paste0("./Results/Tables/",
                 "Backend-All-JaccardIndex-Summary.csv"),
          row.names=FALSE)

head(JList5)
JList5Tidy=JList5[names(JList5)!="annot"]
names(JList5Tidy)=c("Comparison", "fileAbbrv1", "fileAbbrv2", "jIndex0.05",
                    "jIndex0.1", "jIndex0.17")
write.csv(JList5Tidy, paste0("./Results/Tables/",
                             "BackFrontend-All-JaccardIndex-Summary.csv"),
          row.names=FALSE)


#End of script#
