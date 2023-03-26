# Find nearest reference in situ sediment flux station.
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
#### 1A. IMPORT REFERENCE DATA ####
# Baseline data from Milliman and Farnsworth, 2012
milliman_river_database <- fread(paste0(
  wd_imports,'milliman-2012-sediment-database.csv'))[,':='(
    Latitude_mil = LATITUDE,
    Longitude_mil = LONGITUDE,
    LATITUDE = NULL,
    LONGITUDE = NULL)
  ]

#### 1B. IMPORT METADATA FOR LANDSAT OBSERVATIONS ####
# River metadata, including site ID, Lat/Long, drainage area, width, number of images, number of pixels
outlet_rivers_ls_metadata <- fread(paste0(wd_imports, 'outlet_rivers_ls_metadata.csv'))


# Get duplicated Milliman and Farnsworth rivers (based on RiverName field)
# Not sure this does anything
# milliman_duplicate_names <- milliman_river_database[RiverName %chin% milliman_river_database[which(duplicated(milliman_river_database$RiverName))]$RiverName]

#### 2. DISTANCE TO REFERENCE STATION ####
# For all Landsat SSC sites, find nearest station in the Milliman and Farnsworth database
# !! REQUIRED !! data.table of stations with columns: Latitude, Longitude, drainage_area_km2
getNearestStn <- function(site_no_sel){
  dt <- outlet_rivers_ls_metadata[site_no == site_no_sel]
  dt_lat <- dt$Latitude
  dt_long <- dt$Longitude
  dt_drainage_area <- dt$drainage_area_km2

    # River drainage area must be within ~30% of Milliman and Farnsworth drainage area
  closest_stns <- milliman_river_database[
    (Area * 1000) > (dt_drainage_area * 0.2) &
      (Area * 1000) < (dt_drainage_area * 5)]
  # print(paste0(round(dt_lat, 3), ', ', round(dt_long,3))) # QA
  # print(closest_stns[1:3]) # QA
  # Calculate distance
  # As written, keeps all reference sites that meet criteria.
  closest_stns <- closest_stns[
        ,':='(distance_deg = sqrt((Latitude_mil - dt_lat)^2 + (Longitude_mil - dt_long)^2))][
          order(distance_deg)
        ][
          # 1 # Should this just be the closest station ([1], or more than one station [1:N])
        ]
  
  min_dist <- min(closest_stns$distance_deg, na.rm = T)
  # print(min_dist) # QA
  closest_stn <- closest_stns[,':='(relative_dist = (distance_deg - min_dist)/min_dist)]
  closest_stn <- closest_stns[distance_deg < 6 & relative_dist < 0.25]
  closest_stn <- cbind(dt, closest_stn)
  return(closest_stn)
}

# Find nearest station: apply function to list of sites
# Takes about 5 seconds for 500 sites
ls_site_nos <- outlet_rivers_ls_metadata$site_no
nearby_stns <- rbindlist(lapply(ls_site_nos, getNearestStn),
                         use.names = T, fill = T)

# Remove Milliman data if longer than threshold distance
# Without this step, reference stations from far away would get applied to rivers
threshold_distance_dd <- 6
milliman_cols <- colnames(milliman_river_database)
nearby_stns_clean <- nearby_stns[distance_deg > threshold_distance_dd, (milliman_cols) := NA]

fwrite(nearby_stns_clean[,.(site_no, Latitude, Longitude, drainage_area_km2, ID, RiverName, n_imgs)],
       file = paste0(wd_imports, 'river_name_id_table.csv'))

# Import sites with corrected IDs
# I corrected this table manually, identifying duplicates and special characters in names (Evan Dethier, 2021)
nearby_stns_ID_correction <- fread(paste0(wd_imports,'river_name_id_table_corrected.csv'))
# Select only sites that need to be corrected (River Name is wrong, or no coordinates, or tributary/no river sites)
nearby_stns_ID_correction_filtered <- nearby_stns_ID_correction[
              !(ID_corrected == 'None' | ID_corrected == '') &
              !is.na(Latitude) & !(Main_Trib %chin% c('Trib', 'NoRiver'))]

# Convert ID column to numeric for merging, select columns
nearby_stns_ID_correction_filtered <- nearby_stns_ID_correction_filtered[,':='(
              ID = as.numeric(ID_corrected))][
                ,.(site_no, ID, drainage_area_km2, Latitude, Longitude, n_imgs)]

# Merge with river database table using corrected ID
nearby_stns_ID_replace <- milliman_river_database[nearby_stns_ID_correction_filtered, on = 'ID']

# Select sites that do not need to be corrected, also excluding sites that have duplicates due to site correction from previous step
nearby_stns_unaltered <- nearby_stns_ID_correction[ID_corrected == '' & !is.na(Latitude) &
                                                     !(ID %in% nearby_stns_ID_replace$ID)]

# Identify duplicate sites in full reference database
# Make table for subsequent binding
nearby_stns_all <- rbind(nearby_stns_ID_replace, 
                         milliman_river_database[
                           nearby_stns_ID_correction[ID_corrected == '' & !is.na(Latitude)][
                             ,.(site_no, ID, drainage_area_km2, Latitude, Longitude, n_imgs)],
                           on = 'ID'], 
                         use.names = T, fill = T)

# Find duplicate name/id combinations
true_duplicates <- duplicated(nearby_stns_all[,.(site_no, ID)])

# Remove duplicates from metadata table
nearby_stns_all_clean <- nearby_stns_all[!true_duplicates]

# Find sites without a match in the Milliman dataset
unmatched_sites <- outlet_rivers_ls_metadata[!(site_no %chin% nearby_stns_all_clean$site_no)]

# Save unmatched sites to file
fwrite(unmatched_sites, file = paste0(wd_imports, 'unmatched_sites.csv'))

# Import table with corrected names
# I corrected this table manually, identifying duplicates and special characters in names (Evan Dethier, 2021)
unmatched_sites_corrected_raw <- fread(paste0(wd_imports,'unmatched_sites_corrected.csv'))

# Update unmatched sites in database with names from corrected names table
unmatched_sites_corrected <- milliman_river_database[
  unmatched_sites_corrected_raw[!(ID %chin% c('None', 'New')),
                                .(ID, site_no, Latitude, Longitude, drainage_area_km2, n_imgs)][
                                  ,':='(ID = as.numeric(ID))], 
  on = 'ID']

unmatched_sites_new <- unmatched_sites_corrected_raw[ID == 'New',
                                                     .(ID, site_no, Latitude, Longitude, drainage_area_km2, n_imgs, RiverName)]

# Add all sites together into one complete reference data table
nearby_stns_all_clean <- rbind(nearby_stns_all_clean, 
                               unmatched_sites_corrected,
                               unmatched_sites_new,
                               use.names = T, fill = T)

# Add distance in decimal degrees to each site
nearby_stns_all_clean <- nearby_stns_all_clean[,':='(distance_deg = sqrt((Latitude_mil - Latitude)^2 + (Longitude_mil - Longitude)^2))]

#### 3. EXPORT DATA FOR NEXT STEPS ####
fwrite(nearby_stns_all_clean, file = paste0(wd_imports, 'nearby_stns_all_clean.csv'))

#### 4. EXPORT VISUAL OF RIVER MATCH-UPS ####
# Plot distance to nearest station
site_distance_map <- ggplot(nearby_stns_all_clean[!is.na(Continent_Region)], 
                            aes(x = Longitude, y = Latitude, color = distance_deg)) +
  # geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
  #              fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(size = 2) + 
  geom_point(aes(x = Longitude_mil, y = Latitude_mil), color = 'black', size = 1) +
  geom_segment(aes(x = Longitude_mil, y = Latitude_mil, xend = Longitude , yend = Latitude),
               color = 'black', size = 0.1) +
  season_facet +
  theme(legend.position = c(0.1, 0.3)) +
  scale_color_gradientn(colors = c('#CABFD9','#565902','#BF9004')) +
  labs(
    x = 'Longitude',
    y = 'Latitude',
    color = 'Distance\n(deg)'
  )

site_distance_map_labeled <- site_distance_map + 
  geom_label_repel(aes(label = RiverName), size = 2) +
  facet_wrap(.~Continent_Region, scales = 'free')

ggsave(site_distance_map_labeled, filename = paste0(wd_figures, 'Site_distance_map.png'),
       width = 20, height = 20)
