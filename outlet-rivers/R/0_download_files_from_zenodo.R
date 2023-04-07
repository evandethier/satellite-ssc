# Download raw suspended sediment data from Zenodo repository.
# If you run this code from inside an R project, 
# it should download files and set-up the folder structure
# necessary for running the remaining code files (1-10).
# 
# Note: you will need to install the packages loaded in section `i. LIBRARY IMPORTS`
# to run the remaining .R files. For this download, you only need the `zen4R` package.
# 
# Send questions to Evan Dethier
# evan.n.dethier (at) g-mail

#### i. LIBRARY IMPORTS ####
# Tables
library(data.table)
library(tidyverse)
library(lubridate)

# Analysis
library(glmnet)
library(Kendall)

# Plots
library(ggplot2)
library(scales)
library(ggthemes)
library(ggpubr)
library(ggrepel)
library(gstat)
library(ggtext)
library(maps)
library(patchwork)
library(egg)

# Download data from USGS and Zenodo
library(dataRetrieval)
library(zen4R)

#### ii. THEMES ####
# Import world map
world_map <- data.table(map_data(map = 'world'))

## Colors
# Set colors for plotting
# decrease_color <- '#2F9FD2' # More subtle blue
decrease_color <- '#0076AD'
increase_color <- '#FA9E2F'

## Plotting themes
theme_clean <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = 'dashed',color = 'grey70'),
    panel.grid.major.x = element_blank(),
    # panel.grid = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(size = 0.5),
    text = element_text(size=8),
    axis.text = element_text(size = 8), 
    plot.title = element_text(size = 9)
  )

theme_facet <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    # panel.grid = element_blank(),
    # legend.position = 'none',
    panel.border = element_rect(size = 0.5),
    strip.background = element_rect(fill = 'white'),
    text = element_text(size=12),
    axis.text = element_text(size = 12), 
    plot.title = element_text(size = 13)
  )
season_facet <- theme_facet + theme(
  legend.position = 'none', 
  strip.background = element_blank(),
  strip.text = element_text(hjust = 0, margin = margin(0,0,0,0, unit = 'pt'))
)

## Plotting axis formats
# Converts log10 axis values to format 10^x
fancy_scientific_modified <- function(l) { 
  # turn in to character string in scientific notation 
  if(abs(max(log10(l), na.rm = T) - min(log10(l), na.rm = T)) >= 2 | 
     # min(l, na.rm = T) < 0.01 | 
     max(l, na.rm = T) > 1e5){ 
    l <- log10(l)
    label <- parse(text = paste("10^",as.character(l),sep = ""))
  }else{
    label <- parse(text = paste(as.character(l), sep = ""))
  }
  # print(label)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(label)
}

# Format latitude for figures
lat_dd_lab <- function(l){
  label <- c()
  for(i in 1:length(l)){
    label_sel <- ifelse(l[i] < 0, paste0(abs(l[i]), '°S'), 
                        paste0(abs(l[i]), '°N'))
    label <- c(label, label_sel)
  }
  return(label)}

# Format longitude for figures
long_dd_lab <- function(l){
  label <- c()
  for(i in 1:length(l)){
    label_sel <- ifelse(l[i] < 0, paste0(abs(l[i]), '°W'), 
                        paste0(abs(l[i]), '°E'))
    label <- c(label, label_sel)
  }
  return(label)}

# Format years with abbreviation (1990 becomes '90)
year_abbrev_lab <- function(l){
  label <- c()
  for(i in 1:length(l)){
    label_sel <- paste0("'", substr(l[i], 3,4))
    label <- c(label, label_sel)
  }
  return(label)}

# Set up log-10 breaks for plot axes
breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#### iii. SET DIRECTORIES ####
# Set root directory
wd_root <- getwd()

# Imports folder (store all import files here)
wd_imports <- paste0(wd_root,"/outlet-rivers/imports/")
# Exports folder (save tables and reports here. May remove if replaced by wd_figures)
wd_exports <- paste0(wd_root,"/outlet-rivers/exports/")

wd_figures <- paste0(wd_root, "/outlet-rivers/figures/")

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_imports, wd_exports, wd_figures)
for(i in 1:length(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}
# Set root directory
wd_root <- getwd()

# Imports folder (store all import files here)
wd_imports <- paste0(wd_root,"/outlet-rivers/imports/")
# Exports folder (save tables and reports here. May remove if replaced by wd_figures)
wd_exports <- paste0(wd_root,"/outlet-rivers/exports/")

wd_figures <- paste0(wd_root, "/outlet-rivers/figures/")

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_imports, wd_exports, wd_figures)
for(i in 1:length(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}

#### 1. SET-UP AND EXECUTE DOWNLOAD(S) ####
#### 1A. SET-UP DOWNLOAD(S) ####
# Specify Zenodo DOI
zenodo_doi <- '10.5281/zenodo.7772047'

# Select certain Zenodo files

# All files needed to create project, split into sediment files and discharge files
# (RECOMMENDED APPROACH)
zenodo_files_sediment <- 'evandethier_2022_global_sediment_flux_required_sediment_files.zip'
zenodo_files_discharge <- c('africa.zip','asia.zip','europe.zip','namerica.zip','pacific.zip','samerica.zip')
zenodo_files <- c(zenodo_files_sediment,zenodo_files_discharge)

## OTHER OPTIONS: DOWNLOAD ALL FILES USED IN PROJECT. ##
# CAREFUL: IF YOU USE THIS, YOU'LL NEED TO MODIFY THE DOWNLOAD CODE
# If you don't want to run the code from scratch, I suggest downloading manually from Zenodo.
# Then you can skip straight to `8_figures_and_final_calculations.R` OR
# Analyze the files as you wish.

# # All files in project (even those created during project)
# zenodo_files <- 'evandethier_2022_global_sediment_flux_all_imported_files.zip'

#### 1B. DOWNLOAD FROM ZENODO ####
# Download the files from Zenodo
# Stored as a zip file in folder specified by `path`

# First set timeout to longer so files will download
getOption('timeout')
options(timeout=400)

# Then download files
download_zenodo(
  doi = zenodo_doi,
  path = wd_imports,
  files = zenodo_files
)

#### 2. UNZIP DOWNLOADED FILES AND ARRANGE IN FILE STRUCTURE ####
# Unzip all downloaded files
# Arrange in folder structure for next .R files
for(i in 1:length(zenodo_files)){
  # Select file
  zipfile_sel <- zenodo_files[i]
  # Downloaded file name doesn't have the 'evandethier' prefix, so remove it
  zenodo_file_download <- paste0(wd_imports, gsub('evandethier/', '', zipfile_sel))
  
  # Determine whether to add files to subfolder
  # (Discharge data go in subfolders, everything else stays out)
  exdir_sel <- ifelse(zipfile_sel %chin% zenodo_files_discharge, 
                      paste0(wd_imports, gsub('.zip', '', zipfile_sel)), 
                      wd_imports)
  
  # Unzip and store in main folder or subfolders
  unzip(zenodo_file_download, # pathname of the zip file 
        exdir = exdir_sel,    # pathname to put extracted files (subfolder if discharge)
        junkpaths = T)        # Just adds the files without making any sub-folders
}

