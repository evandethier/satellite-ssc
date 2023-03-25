#### CACLULATE MONTHLY DISCHARGE TIMESERIES
# Calculate monthly discharge for each river with satellite-derived suspended sediment data
# Takes about 20 minutes

# Only run from scratch if `discharge_month_yr_data_all.csv` AND/OR
# Only run from scratch if `discharge_metadata_all.csv`
# does not yet exist or needs to be updated.

# Right now, the code that creates these files is commented out
# **To run it, comment in the lines in section 7**

# `discharge_month_yr_data_all.csv` and `discharge_metadata_all.csv` 
# can be downloaded from:


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
#### 1A. IMPORT REFERENCE SITES ####
# Reference data from Milliman and Farnsworth, 2012
# Reference data has been merged with Landsat-derived metadata, from this study (Dethier, 2022)
# Discharge data from the Global Runoff Data Center (GRDC), first established in Lehner et al., 2002
# Landsat data identified with unique site id, stored in the `site_no` column
nearby_stns_all_clean <- fread(paste0(wd_imports, 'nearby_stns_all_clean.csv'))

#### 1B. IMPORT SEASONAL DISCHARGE DATA #### 
# This takes a while, so best to do just once
# Make a vector of possible metadata column names
# Possible column names
col_possible <- c('GRDC-No.', 'River','Station','Country','Latitude (DD)','Longitude (DD)',
                  'Catchment area ()', 'Altitude (m ASL)')
# Replacement column names
col_replace <- c('site_no_GRDC', 'RiverName_GRDC','Station_GRDC','Country','Latitude_GRDC','Longitude_GRDC',
                 'drainage_area_km2', 'Elevation_m')


#### 1C. SET UP MONTHLY DISCHARGE IMPORT FUNCTION ####
getGRDC_monthly <- function(file_sel, folder_sel){
  # Import select file without metadata
  q_import_step1 <- fread(paste0(folder_sel,'/', file_sel),
                          fill = T)
  # Import select file to preserve metadata columns
  q_import_metadata <- fread(paste0(folder_sel,'/', file_sel),
                             fill = T, sep = ':', header = T)
  
  if(colnames(q_import_metadata)[2] == 'GRDC STATION DATA FILE'){
    # Number of columns in data file
    import_ncol <- ncol(q_import_step1)
    # Vector of column headers from import (default VX naming)
    col_headers <- paste0('V',3:import_ncol)
    
    # Transpose first 15 rows of metadata table to columns
    # Each column is a single metadata value
    import_metadata <- data.table(t(q_import_metadata[1:20, `GRDC STATION DATA FILE`]))
    # Rename column names with column headers transposed from metadata file "Title" column
    colnames(import_metadata) <- gsub('# |km\xb2', '', t(q_import_metadata[1:20, `# Title`]))
    
    # Clean up monthly discharge data
    # Monthly discharge data has a semicolon in the row, so select only those rows from the full table
    # The header for the monthly discharge is MM;, so don't select that column
    # Convert default column headers V1 (month) and V2 (discharge) to Month and Q_cms headers
    monthly_q_cms <- data.table(q_import_step1)[grepl(';',V1) & V3 == '' & V1 != 'MM;'][
      ,.(V1, V2)][
        ,':='(month = as.numeric(gsub(';','',V1)),
              Q_cms = as.numeric(gsub(';','',V2)),
              V1 = NULL,
              V2 = NULL)]
    
    # ID column names from the master list of possible columns
    col_sel <- colnames(import_metadata)[
      colnames(import_metadata) %chin% col_possible]
    
    # Only select columns from master list
    metadata <- import_metadata[,..col_sel]
    # Rename columns using replacement list
    colnames(metadata) <- col_replace[which(col_sel %in% col_possible)]
    
    # Bind together tables, include latitude and longitude from GRDC station
    monthly_q_md <- cbind(monthly_q_cms, metadata)[
      ,':='(Latitude_GRDC = as.numeric(Latitude_GRDC),
            Longitude_GRDC = as.numeric(Longitude_GRDC),
            drainage_area_km2 = as.numeric(drainage_area_km2))
    ]
    
    return(monthly_q_md)
  }}

# Folders with runoff data
all_folders <- list.dirs(wd_imports)
runoff_folders <- all_folders[all_folders != wd_imports]

#### 1D. APPLY FUNCTION TO IMPORT FROM DAILY DATABASE ####
# Loop through folders with daily-derived data in them
# Import and process monthly data files
for(i in 1:length(runoff_folders)){
  
  runoff_folder_sel <- runoff_folders[i] 
  # List of all files
  continent_files <- list.files(runoff_folder_sel)
  # List of LTVD files
  continent_files_ltvd <- continent_files[which(grepl('LTVD', continent_files))]
  
  for(j in 1:length(continent_files_ltvd)){
    file_sel <- continent_files_ltvd[j]
    import_monthly <- getGRDC_monthly(file_sel, runoff_folder_sel)
    if(i == 1 & j == 1){
      
      import_monthly_all <- import_monthly
    }else{
      import_monthly_all <- rbind(import_monthly_all, import_monthly)
    }
  }
}

#### 1E. APPLY FUNCTION TO IMPORT FROM MONTHLY DATABASE ####
# Loop through folders with monthly-derived data in them
# Import and process monthly data files
for(i in 1:length(runoff_folders)){
  
  runoff_folder_sel <- runoff_folders[i] 
  # List of all files
  continent_files <- list.files(runoff_folder_sel)
  # List of LTVM files
  continent_files_ltvm <- continent_files[which(grepl('LTVM', continent_files))]
  
  for(j in 1:length(continent_files_ltvm)){
    file_sel <- continent_files_ltvm[j]
    import_monthly_mean <- getGRDC_monthly(file_sel, runoff_folder_sel)
    if(i == 1 & j == 1){
      
      import_monthly_mean_all <- import_monthly_mean
    }else{
      import_monthly_mean_all <- rbind(import_monthly_mean_all, import_monthly_mean)
    }
  }
}

ggplot(monthly_q_all_stns[,.(lat = mean(as.numeric(Latitude_GRDC), na.rm = T), 
                             long = mean(as.numeric(Longitude_GRDC), na.rm = T)),
                          by = .(site_no_GRDC)]) + 
  geom_point(aes(x = long, y = lat))

#### 2. COMBINE DAILY- AND MONTHLY-DERIVED DATASETS ####
# Preference is for data derived from daily measurements

# Merge daily-derived and monthly-derived discharge observations
monthly_q_all_stns <- import_monthly_all[
  import_monthly_mean_all[,':='(Q_cms_mnth = Q_cms)][
    ,.(month, Q_cms_mnth, site_no_GRDC)],
  on = c('month', 'site_no_GRDC')]

# Replace values missing from daily-calculated discharge (Step B)
# With values calculated from monthly-calculated mean discharge (Step C)
monthly_q_all_stns <- monthly_q_all_stns[,':='(
  Q_cms_dly = Q_cms,
  Q_cms = ifelse(Q_cms > 0, Q_cms, Q_cms_mnth)
)]

#### 3. IMPORT ANNUAL DISCHARGE ####

#### 3A. SET UP FUNCTION FOR ANNUAL IMPORT ####
getGRDC_annual <- function(file_sel, folder_sel){
  file_sel_yvd <- gsub('YVM','YVD',file_sel)
  # Import select YVD file without metadata
  q_import_step1_yvm <- fread(paste0(folder_sel,'/', file_sel),
                              fill = T)
  # Import select YVM file without metadata
  q_import_step1_yvd <- fread(paste0(folder_sel,'/', file_sel_yvd),
                              fill = T)
  # Import select file to preserve metadata columns
  q_import_metadata <- fread(paste0(folder_sel,'/', file_sel),
                             fill = T, sep = ':', header = T)
  
  
  
  if(colnames(q_import_metadata)[2] == 'GRDC STATION DATA FILE'){
    
    # Transpose first 15 rows of metadata table to columns
    # Each column is a single metadata value
    import_metadata <- data.table(t(q_import_metadata[1:20, `GRDC STATION DATA FILE`]))
    # Rename column names with column headers transposed from metadata file "Title" column
    colnames(import_metadata) <- gsub('# |km\xb2', '', t(q_import_metadata[1:20, `# Title`]))
    
    # ID column names from the master list of possible columns
    col_sel <- colnames(import_metadata)[
      colnames(import_metadata) %chin% col_possible]
    
    # Only select columns from master list
    metadata <- import_metadata[,..col_sel]
    # Rename columns using replacement list
    colnames(metadata) <- col_replace[which(col_sel %in% col_possible)]
    
    # Clean up annual discharge data
    # Annual discharge data has a "YYYY;" in the row, so select only those rows from the full table
    # The header for the monthly discharge is "YYYY;", so don't select that column
    # Convert default column headers V1 (month) and V2 (discharge) to Month and Q_cms headers
    # Find row with mean value header
    
    mean_row_num_yvd <- which(q_import_step1_yvd$V1 == 'YYYY;')[3]
    mean_row_num_yvm <- which(q_import_step1_yvm$V1 == 'YYYY;')[2]
    
    # Find number of years with data (two lines before header row)
    n_years_yvd <- as.numeric(gsub(',','',q_import_step1_yvd$V4[(mean_row_num_yvd-2)]))
    n_years_yvm <- as.numeric(q_import_step1_yvm$V4[(mean_row_num_yvm-2)])
    
    row_range_yvd <- c((mean_row_num_yvd + 1):(mean_row_num_yvd + n_years_yvd))
    row_range_yvm <- c((mean_row_num_yvm + 1):(mean_row_num_yvm + n_years_yvm))
    # Parse annual data into a clean, labeled datatable
    annual_q_cms_yvd <- data.table(q_import_step1_yvd)[row_range_yvd][,
                                                                      .(V1, V2)][,':='(year = as.numeric(gsub(';','',V1)),
                                                                                       Q_cms_yvd = as.numeric(gsub(';','',V2)),
                                                                                       meas_interval = 'Daily',
                                                                                       V1 = NULL,
                                                                                       V2 = NULL)]
    
    annual_q_cms_yvm <- data.table(q_import_step1_yvm)[row_range_yvm][,
                                                                      .(V1, V2)][,':='(year = as.numeric(gsub(';','',V1)),
                                                                                       Q_cms_yvm = as.numeric(gsub(';','',V2)),
                                                                                       meas_interval = 'Monthly',
                                                                                       V1 = NULL,
                                                                                       V2 = NULL)]
    
    # Combine monthly and daily discharge data, selecting daily when available
    annual_q_cms <- merge(merge(data.table(year = c(1900:2020),
                                           Q_cms = NA),
                                annual_q_cms_yvd[Q_cms_yvd >= 0][,.(year, Q_cms_yvd)], all.x = T, by = 'year'),
                          annual_q_cms_yvm[Q_cms_yvm >= 0][,.(year, Q_cms_yvm)], all.x = T, by = 'year')[
                            ,':='(Q_cms = ifelse(!is.na(Q_cms_yvd), Q_cms_yvd, Q_cms_yvm),
                                  Q_cms_yvd = NULL,
                                  Q_cms_yvm = NULL
                            )
                          ]
    
    # Fill missing record with mean values
    annual_q_cms_avg <- mean(annual_q_cms$Q_cms, na.rm = T)
    
    annual_q_cms <- annual_q_cms[,':='(Q_cms_avg = annual_q_cms_avg)]
    for(n in c(3, 5, 10, 30)){
      col_sel <- c(paste0('Q_cms_',n,'yr'))
      annual_q_cms_avg_roll <- data.table(frollmean(annual_q_cms$Q_cms, n = n, 
                                                    align = 'center', na.rm = T))
      colnames(annual_q_cms_avg_roll) <- col_sel
      annual_q_cms <- cbind(annual_q_cms, annual_q_cms_avg_roll)
    }
    
    annual_q_cms <- annual_q_cms[
      ,':='(Q_cms_sel = ifelse(!is.na(Q_cms), Q_cms,
                               ifelse(!is.na(Q_cms_3yr), Q_cms_3yr,
                                      ifelse(!is.na(Q_cms_5yr), Q_cms_5yr,
                                             ifelse(!is.na(Q_cms_10yr), Q_cms_10yr, Q_cms_avg)))),
            Q_cms_period = ifelse(!is.na(Q_cms), 'Annual',
                                  ifelse(!is.na(Q_cms_3yr), 'Avg_3yr',
                                         ifelse(!is.na(Q_cms_5yr), 'Avg_5yr',
                                                ifelse(!is.na(Q_cms_10yr), 'Avg_10yr', 'Avg')))))]
    
    # Get final 10-year running avg. estimate to fill all unfilled observations at end of record
    if(sum(annual_q_cms$Q_cms_10yr,na.rm = T) > 0){
      final_10yr_yr <- max(annual_q_cms[!is.na(Q_cms_10yr)]$year)-5
      final_10yr_obs <- annual_q_cms[year == final_10yr_yr]$Q_cms_10yr
      annual_q_cms <- annual_q_cms[year > final_10yr_yr, 
                                   ':='(Q_cms_sel = final_10yr_obs,
                                        Q_cms_period = 'Avg_10yr_end')]
    }
    
    # Select columns to write
    annual_q_cms_final <- annual_q_cms[,.(year, Q_cms_sel, Q_cms_period)][
      ,':='(Q_cms = Q_cms_sel,
            Q_cms_sel = NULL)
    ]
    
    annual_q_md <- cbind(annual_q_cms_final, metadata)[
      ,':='(Latitude_GRDC = as.numeric(Latitude_GRDC),
            Longitude_GRDC = as.numeric(Longitude_GRDC),
            drainage_area_km2 = as.numeric(drainage_area_km2))
    ]
    
    return(annual_q_md)
  }}

#### 3B. APPLY FUNCTION FOR ANNUAL IMPORT ####
# IMPORT FROM ANNUAL DATABASE
# First from daily datasets
for(i in 1:length(runoff_folders)){
  
  runoff_folder_sel <- runoff_folders[i]
  # List of all files
  continent_files <- list.files(runoff_folder_sel)
  # List of YVD files
  continent_files_yvm <- continent_files[which(grepl('YVM', continent_files))]
  
  for(j in 1:length(continent_files_yvm)){
    file_sel <- continent_files_yvm[j]
    import_annual_mean <- getGRDC_annual(file_sel, runoff_folder_sel)
    if(i == 1 & j == 1){
      
      import_annual_mean_all <- import_annual_mean
    }else{
      import_annual_mean_all <- rbind(import_annual_mean_all, import_annual_mean)
    }
  }
}


#### 4. JOIN DISCHARGE AND SSC STATIONS ####
# Every SSC station needs a Q station
# For a given SSC station, find nearby Q stations
nearby_stns_final <- nearby_stns_all_clean
site_nos_all <- nearby_stns_final[RiverName != '']$site_no

# Name changes:
# Hawkesbury should be Nepean R. (Australia)
# East Fitzroy should be Fitzroy (Australia)
# Murray-Darling should be Murray

site_no_sel <- 'st_0000000000000000059b'

#### 4A. SET UP FUNCTION TO FIND NEAREST DISCHARGE STATION ####
getNearest_Q_stn <- function(site_no_sel){
  # For each SSC Landsat site, get nearby sites with comparable drainage area
  dt <- nearby_stns_final[RiverName != ''][site_no == site_no_sel]
  dt_name <- dt$RiverName
  dt_drainage_area <- dt$drainage_area_km2
  
  # Get distance to each site with similar name
  dt_lat <- dt$Latitude
  dt_long <- dt$Longitude
  closest_stn <- monthly_q_all_stns[
    # drainage_area_km2 > 0.2 * dt_drainage_area &
    # drainage_area_km2 < 5 * dt_drainage_area
  ][
    ,':='(distance_deg = sqrt((Latitude_GRDC - dt_lat)^2 + 
                                (Longitude_GRDC - dt_long)^2))][
                                  order(distance_deg)
                                ][distance_deg < 6][
                                  ,':='(
                                    name_match = agrepl(dt_name, RiverName_GRDC, ignore.case = T),
                                    name_dist = adist(dt_name, RiverName_GRDC, ignore.case = T))
                                ]
  
  min_name_dist <- min(closest_stn[name_match == TRUE]$name_dist)
  
  # If no name + proximity match can be found
  # Use proximity and drainage area correspondence to find nearby station
  if(min_name_dist == Inf){
    # Filter to only nearby stations with similar discharge
    nearby_q_stns <- closest_stn[
      drainage_area_km2 > 0.2 * dt_drainage_area &
        drainage_area_km2 < 5 * dt_drainage_area
    ]
    # Get minimum distance of nearby stations
    nearby_stns_min_dist <- min(nearby_q_stns$distance_deg, na.rm = T)
    # Select only rows that are the minimum distance (from one station)
    closest_stn_sel <- nearby_q_stns[distance_deg == nearby_stns_min_dist][
      ,':='(match_type = 'Nearby',
            drainage_area_km2_GRDC = drainage_area_km2)
    ]
  }else{
    # Standard approach with near-exact name match + proximity
    closest_stn_sel <- closest_stn[name_match == TRUE & name_dist <= min_name_dist*1.5][
      ,':='(match_type = 'Exact',
            drainage_area_km2_GRDC = drainage_area_km2)
    ]
  }
  # Bind with selected columns from SSC landsat metadata
  closest_stn_final <- cbind(dt[,.(site_no, drainage_area_km2, RiverName, Latitude, Longitude)], closest_stn_sel[,-c('drainage_area_km2')])
  return(closest_stn_final)
}

#### 4B. APPLY FUNCTION TO FIND NEAREST DISCHARGE STATION ####
# Loop through all sites in SSC Landsat database to find nearest Runoff Station
for(i in 1:length(site_nos_all)){
  site_no_sel <- site_nos_all[i]
  ID_sel <- nearby_stns_final[site_no == site_no_sel]$ID
  # Apply getNearest_Q_stn function to each site
  nearest_q_stn <- getNearest_Q_stn(site_no_sel)[,':='(ID = as.character(ID_sel))]
  
  # Write to a combined table for future analysis
  if(i == 1){
    ssc_q_stns <- nearest_q_stn
  }else{
    ssc_q_stns <- rbind(ssc_q_stns, nearest_q_stn, fill = T)
  }}

# Number of stations with exact name match vs. nearby match
n_stations_exact_vs_nearby <- ssc_q_stns[
  , .(`Runoff Match Type` = uniqueN(RiverName)),
  by = match_type]

#### 5. ASSESS SUCCESS OF DISCHARGE/SSC SITE MATCHING, MANUALLY CLEAN UP ####
# To limit the size of the data table, only use outlet river discharge stations
ssc_q_stns <- ssc_q_stns[,':='(unique_id = paste0(RiverName, '_', site_no_GRDC))]
# Make vectors of:
# a) GRDC station number (this is never used)
# b) Outlet river stations 
outlet_river_grdc <- unique(ssc_q_stns$unique_id)
outlet_river_names <- unique(ssc_q_stns[!is.na(Q_cms)]$RiverName)

# Combine GRDC station name and number into single vector
# This is never used
unique_river_names_site_nos <- unique(ssc_q_stns[!is.na(Q_cms)][,.(paste0(RiverName, site_no))])

# Make a table of unique discharge stations and metadata
ssc_q_stns_unique <- ssc_q_stns[,.(distance_deg = mean(distance_deg, na.rm = T),
                                   drainage_area_km2 = mean(drainage_area_km2, na.rm = T)),
                                by = .(site_no, ID, RiverName, RiverName_GRDC, match_type, Latitude, Longitude)]

# Stations without an exact match (wrong name or too far from GRDC station)
ssc_q_stns_nearby <- ssc_q_stns_unique[match_type != 'Exact', 
                                       .(ID, RiverName, RiverName_GRDC, match_type, Latitude, Longitude,distance_deg, drainage_area_km2)]

# Make a table of GRDC stations (this is never used)
GRDC_stations <- monthly_q_all_stns[,.(Q_cms_avg = mean(Q_cms, na.rm = T)),
                                    by = .(site_no_GRDC, RiverName_GRDC, Station_GRDC, Country,
                                           Latitude_GRDC, Longitude_GRDC, drainage_area_km2, Elevation_m)]

# Write metadata files to drive
fwrite(ssc_q_stns_nearby, file = paste0(wd_imports, 'SSC_stations_missing_Q_station.csv'))
fwrite(GRDC_stations, file = paste0(wd_imports, 'GRDC_stations.csv'))

# I made manual name corrections to this file -Evan Dethier, Aug. 2021
# This file needs to be loaded into directory or created manually
ssc_q_stns_nearby_corrected <- fread(paste0(wd_imports, 'SSC_stations_missing_Q_station_corrected.csv'))

# Make a vector of site names for which correcting name resulted in an exact match
ssc_q_stns_exact_name_corrected <- ssc_q_stns_nearby_corrected[
  match_type_corrected == 'Exact']$ID

# Update number of exact matches in SSC/Discharge metadata table
ssc_q_stns_unique <- ssc_q_stns_unique[ID %chin% ssc_q_stns_exact_name_corrected,
                                       ':='(match_type = 'Exact')]

# Summarize exact vs. nearby matches (and no matches)
n_stations_exact_vs_nearby_corrected <- ssc_q_stns_unique[
  , .(`Runoff Match Type` = uniqueN(ID)),
  by = match_type][
    ,':='(`% Total` = round((`Runoff Match Type`/nrow(ssc_q_stns_unique))*100))
  ]


#### 6. JOIN MONTHLY AND ANNUAL DISCHARGE TIMESERIES, INTERPOLATE MISSING DATA ####
i <- 499
i <- 463
i <- 403
i <- 357
i <- 367

# For each river, get monthly discharge
# Replace missing months with interpolated data from:
# a) computing likely monthly discharge based on annual total OR
# b) taking long-term running average
# (if there is data within 3-yrs, use that, otherwise use 5-yr or 10-year avg.)
for(i in 1:length(outlet_river_names)){
  # Selected river name
  river_sel <- outlet_river_names[i]
  # print(i)
  # Selected station
  # station_sel <- outlet_river_grdc[grepl('Yukon', outlet_river_grdc)]
  station_sel <- unique(ssc_q_stns[RiverName == river_sel]$unique_id)
  # Monthly data (may include multiple stations)
  monthly_sel <- ssc_q_stns[unique_id %in% station_sel]
  
  # ggplot(monthly_sel, aes(x = month, y = Q_cms)) +
  #   geom_line(aes(group = site_no_GRDC)) +
  #   season_facet
  
  # Compute avg. annual Q from monthly Q as a reference
  Q_cms_avg_monthly <- monthly_sel[,.(Q_cms = mean(Q_cms, na.rm = T)),
                                   by = site_no_GRDC]
  
  # Only select sites with > 80% of largest discharge
  # Ensures discharge sites are approximately the right scale for SSC estimates
  largest_sites <- Q_cms_avg_monthly[
    Q_cms > max(Q_cms_avg_monthly$Q_cms, na.rm = T)*0.8]
  # Monthly data stations
  unique_grdc_sel <- largest_sites$site_no_GRDC
  # Annual data
  annual_sel <- import_annual_mean_all[site_no_GRDC %chin% unique_grdc_sel]
  
  # ggplot(annual_sel, aes(x = year, y = Q_cms)) +
  #   geom_line(aes(group = site_no_GRDC)) +
  #   geom_point(aes(color = Q_cms_period)) +
  #   season_facet +
  #   theme(legend.position = c(0.2, 0.8)) +
  #   labs(title = annual_sel$RiverName_GRDC[1])

  # Length of annual timeseries
  annual_ts_length <- annual_sel[Q_cms_period == 'Annual' & Q_cms > 0][
    ,.(n_years = .N), by = site_no_GRDC]
  
  # Get the nearby site with the longest record
  longest_large_site <- annual_ts_length[order(-n_years)]$site_no_GRDC[1]
  Q_cms_avg_monthly_sel <- Q_cms_avg_monthly[site_no_GRDC == longest_large_site]$Q_cms
  
  # Modify annual timeseries with monthly data
  annual_sel <- annual_sel[site_no_GRDC == longest_large_site][
    ,':='(Q_ann_rel = Q_cms/Q_cms_avg_monthly_sel)
  ]
  # Create empty data.table
  # Table with a row for unique SSC site, year (1900-2020), month (1-12)
  month_yr_data <- data.table(year = sort(rep(c(1900:2020), 12)), 
                              month = rep(c(1:12), length(1900:2020)))
  # Join with month data
  month_yr_data <- month_yr_data[
    monthly_sel[site_no_GRDC == longest_large_site][
      order(paste0(match_type, distance_deg))][ # Added to choose closest exact match site, if there are multiple options
        1:12][
          ,':='(Q_cms_month = Q_cms)], on = 'month'][
            ,-c('Q_cms')][
              annual_sel[,.(year, Q_ann_rel, Q_cms_period)], on = 'year']
  
  # Compute monthly Q (scaled by annual total relative to reported monthly avg.)
  month_yr_data <- month_yr_data[,':='(Q_cms = Q_cms_month * Q_ann_rel)]
  
  month_yr_metadata <- month_yr_data[1][,.(site_no, drainage_area_km2, RiverName, Latitude, Longitude, 
                                        site_no_GRDC, RiverName_GRDC, Station_GRDC, Country, Latitude_GRDC, Longitude_GRDC,
                                        Elevation_m, distance_deg, name_match, name_dist, match_type,
                                        drainage_area_km2_GRDC, unique_id)]
  
  month_yr_data <- month_yr_data[,.(site_no, year, month, # Metadata/match data
                                    # Q_cms_mnth, Q_cms_dly, Q_cms_month, # Data from original files
                                    Q_ann_rel, Q_cms_period, Q_cms)] # Aggregated data
  
  # ggplot(month_yr_data[year > 1950 & year < 2020], aes(x = ymd(paste(year, month, 1, sep = '-')), y = Q_cms)) +
  #   geom_line() +
  #   geom_point(aes(color = Q_cms_period)) +
  #   facet_wrap(.~(year - year%%10), scales = 'free_x', ncol = 1) +
  #   season_facet +
  #   theme(legend.position = 'right') +
  #   labs(title = annual_sel$RiverName_GRDC[1])
  
  
  if(i == 1){
    discharge_month_yr_data_all <- month_yr_data
    discharge_metadata_all <- month_yr_metadata
  }else{
    discharge_month_yr_data_all <- rbind(discharge_month_yr_data_all, month_yr_data, fill = T, use.names = T)
    discharge_metadata_all <- rbind(discharge_metadata_all, month_yr_metadata, fill = T, use.names = T)
  }
  
}

#### 7. WRITE DATA TO FILE (IF IT DOESN'T EXIST) ####
# Split file into monthly data and metadata
# Monthly data
# Careful!! This will overwrite an existing file
# fwrite(discharge_month_yr_data_all, paste0(wd_imports, 'discharge_month_yr_data_all.csv'))

# Metadata
# Careful!! This will overwrite an existing file
# fwrite(discharge_metadata_all, paste0(wd_imports, 'discharge_metadata_all.csv'))

#### TEST ####
# continent_files <- list.files(runoff_folder_sel)
# # List of YVD files
# continent_files_yvd <- continent_files[which(grepl('YVD', continent_files))]
# # List of YVM files
# continent_files_yvm <- continent_files[which(grepl('YVM', continent_files))]
# 
# 
# runoff_folder_sel <- paste0(wd_imports,'asia')
# file_sel <- continent_files_yvm[145]
# file_sel <- '2335950_Q_YVM.csv'
# 
# test_runoff <- getGRDC_annual(file_sel, runoff_folder_sel)[]
