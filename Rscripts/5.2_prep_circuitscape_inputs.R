#-------------------------------------------------------------------------------
## 5.2 Prep Circuitscape inputs
## 2020
## Inputs: reclassed resistance (names)
## Outputs: INI Files
#-------------------------------------------------------------------------------

library(stringr)
library(raster)

## Read in 
reclassed_list <- list.files("Data/Processed/",
                             full.names = TRUE,  
                             pattern = "Resistance")

# Make all the INI files
INI_template <- readLines("config/ini_template.ini")

make_INI_file = function(file) {
  base <- tools::file_path_sans_ext(basename(file))
  INI <- lapply(X=INI_template, FUN=str_replace_all, 
                pattern ="RESISTANCE_FILE",
                replacement = file)
  INI <- lapply(X=INI, FUN=str_replace_all, 
                pattern ="OUTPUT_FILE",
                replacement = paste0("outputs/current_density/",
                                     base, 
                                     "_out")) # No need for .tif here
  fileConn<-file(paste0("config/all/", base, ".ini"))
  writeLines(unlist(INI), fileConn)
  close(fileConn)
}

saved_output <- lapply(reclassed_list, make_INI_file)
