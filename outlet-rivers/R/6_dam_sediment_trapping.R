# Calculate timeseries of dam building & reservoir capacity for global outlet rivers.
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

# Download data from USGS
library(dataRetrieval)

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
#### 1. IMPORT DATA ####
#### 1A. IMPORT DAM DATA, CALCULATE CUMULATIVE ANNUAL STORAGE ####
# Import dams data GRAND DAMS dataset
GRAND_dams_all <- fread(paste0(wd_imports,'GRAND_dams.csv'))
# Dams in each outlet river watershed
watershed_dams_import <- fread(paste0(wd_imports,'global-outlet-rivers-dams-08012021.csv'))[
  ,':='(site_no = paste0('st_', site_no))
][site_no %chin% ssc_rivers$site_no]


#### 1B. IMPORT METADATA ####
# Simple metadata
ssc_rivers <- fread(paste0(wd_imports, 'ssc_outlet_rivers_names.csv'),
                    colClasses = c('ID' = 'character'))

# Fuller metadata
ssc_rivers_stats <- fread(paste0(wd_imports, 'ssc_rivers_stats.csv'),
                          colClasses = c('ID' = 'character'))

# Full by-river statistics and metadata
ssc_rivers_stats_full <- fread(paste0(wd_imports, 'ssc_rivers_stats_full.csv'),
                               colClasses = c('ID' = 'character'))

#### TESTING ####
watershed_dams_import[MAIN_BASIN == 'Nile'][,.(N = .N), by = site_no]
ggplot(watershed_dams_import[site_no == 'st_00000000000000001525'], 
       aes(x = dam_long, y = dam_lat, size = CAP_MCM, color = YEAR)) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point() +
  season_facet

#### 2. CALCULATE RESERVOIR STORAGE TIMESERIES ####
#### 2A. RESERVOIR STORAGE TIMESERIES BY RIVER ####
# Calculate cumulative reservoir storage and area for each river, by year
watershed_dams <- watershed_dams_import[
  site_no %chin% ssc_rivers$site_no &
    YEAR > 1800 & YEAR < 2020][
      ssc_rivers, on = c('site_no')][ # Join to add large river name
        ,':='(decade = ifelse(YEAR < 1985, 'Baseline', YEAR - YEAR%%5))
      ][order(YEAR),':='(cumulative_dam_storage_mcm = cumsum(CAP_MCM),
                         cumulative_dam_area_km2 = cumsum(AREA_SKM)),
        by = .(site_no, RiverName)]

# For each river, calculate cumulative dam area for every year (not just every year with a new dam)
# For each river
for(i in 1:nrow(ssc_rivers)){
  site_no_sel <- ssc_rivers$site_no[i]
  riverName_sel <- ssc_rivers$RiverName[i]
  # Select all dams
  river_sel <- watershed_dams[site_no == site_no_sel][
    ,':='(year = ifelse(YEAR < 1900, 1900, YEAR),
          year_known = ifelse(YEAR < 1900, 'unknown', 'known'))
  ]
  # Find year of first dam
  first_yr <- min(river_sel$YEAR, na.rm = T)
  if(!is.infinite(first_yr)){
    
    annual_cap_mcm <- river_sel[,
                                # Set year of all dams without a date to first year
                                ':='(year = ifelse(year < first_yr, first_yr, year))][
                                  # Sum dam area for each year
                                  ,.(cap_mcm = sum(CAP_MCM, na.rm = T),
                                     area_km2 = sum(AREA_SKM, na.rm = T)),
                                  by = .(RiverName, site_no, year)][
                                    ,.(year, cap_mcm, area_km2)
                                  ]
    
    # Make an empty datatable with a row for every year from first year to 2020
    annual_every_cap_mcm <- merge(data.table(year = c(1900:2020)),
                                  annual_cap_mcm,
                                  by = 'year', 
                                  all.x = T)[
                                    ,':='(cap_mcm = ifelse(is.na(cap_mcm), 0, cap_mcm),
                                          area_km2 = ifelse(is.na(area_km2), 0, area_km2))][
                                            ,':='(cap_mcm_cml = cumsum(cap_mcm),
                                                  area_km2_cml = cumsum(area_km2))
                                          ][,':='(site_no = site_no_sel,
                                                  RiverName = riverName_sel)]
    
    
  }else{
    annual_every_cap_mcm <- data.table(year = c(1900:2020))[
      ,':='(cap_mcm = 0, cap_mcm_cml = 0, area_km2 = 0, area_km2_cml = 0, site_no = site_no_sel, RiverName = riverName_sel)
    ]
  }
  
  
  if(i == 1){
    annual_cap_mcm_comb <- annual_every_cap_mcm
  }else{
    annual_cap_mcm_comb <- rbind(annual_every_cap_mcm, annual_cap_mcm_comb, use.names = T, fill = T)
  }
}

# Timeseries of dam building at each river
dam_timeseries <- annual_cap_mcm_comb[year > 1919][
  ssc_rivers_stats_full[,.(site_no, RiverName, Continent_Region)], on = c('site_no','RiverName')][
    ,':='(Continent_Region = ifelse(Continent_Region %in% c('Europe','Eurasia'),'Europe/Eurasia', Continent_Region))]



#### 3. EXPORT DAM DATA TO DRIVE ####
# Timeseries of dam building, by river
fwrite(dam_timeseries, file = paste0(wd_imports, 'dam_timeseries.csv'))




