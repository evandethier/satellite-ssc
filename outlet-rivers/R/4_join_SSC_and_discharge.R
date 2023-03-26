# Join satellite-derived SSC dataset with GRDC discharge data.
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
#### 1. IMPORT DATA: SSC, DISCHARGE, AND REFERENCE DATA #### 
#### 1A. IMPORT LANDSAT-DERIVED SSC DATA AND METADATA ####
# Landsat-derived daily SSC estimates
river_landsat_pred <- fread(paste0(wd_imports, 'river_landsat_pred.csv'))

# Cluster assignments for each river (by month)
river_cluster_assignment_monthly <- fread(paste0(wd_imports, 'river_cluster_assignment_monthly.csv'))

# River metadata, including site ID, Lat/Long, drainage area, width, number of images, number of pixels
outlet_rivers_ls_metadata <- fread(paste0(wd_imports, 'outlet_rivers_ls_metadata.csv'))

# Join with river metadata to have river name and ID associated
# river_landsat_pred <- river_landsat_pred[ssc_q_stns_unique[,.(site_no, ID, RiverName, match_type)],
#                                          on = c('site_no')]  

#### 1B. IMPORT DISCHARGE DATA TO CALCULATE FLUX ####
# Discharge data from nearest GRDC station
# Monthly data for entire period of record (including interpolated data)
discharge_month_yr_data_all <- fread(paste0(wd_imports, 'discharge_month_yr_data_all.csv'))
# Discharge metadata
discharge_metadata_all <- fread(paste0(wd_imports, 'discharge_metadata_all.csv'))
# Join discharge and metadata together
discharge_month_yr_data_all <- discharge_month_yr_data_all[discharge_metadata_all, on = 'site_no']

discharge_month_yr_data_all <- discharge_month_yr_data_all[
  ,':='(
    # decade = ifelse(year < 1990, 1990, 
    #           ifelse(year >= 2020, 2010, year - year%%10)))
    decade = ifelse(year < 1990, 1990,
                    ifelse(year > 2019, 2015, 
                           year - year%%5)))
]

#### 1C. IMPORT REFERENCE RIVER SEDIMENT DATA ####
# Baseline data from Milliman and Farnsworth, 2012
milliman_river_database <- fread(paste0(
  wd_imports,'milliman-2012-sediment-database.csv'))[,':='(
    Latitude_mil = LATITUDE,
    Longitude_mil = LONGITUDE,
    LATITUDE = NULL,
    LONGITUDE = NULL)
  ]

# Get duplicated Milliman and Farnsworth rivers (based on RiverName field)
# Not sure this does anything
# milliman_duplicate_names <- milliman_river_database[RiverName %chin% milliman_river_database[which(duplicated(milliman_river_database$RiverName))]$RiverName]

#### 1D. IMPORT NEAREST REFERENCE STATIONS METADATA TABLE ####
nearby_stns_final <- nearby_stns_final[,':='(ID = as.character(ID))]

#### 2. JOIN SSC WITH MONTHLY DISCHARGE DATA ####
#### 2A. SELECT DISCHARGE COLUMNS AND ADD METADATA ####
# Join discharge data with reference station metadata
# Column names shared by nearby_stns_final, discharge_month_yr_data_all
unique_metadata_col_names <- c('site_no', 
  colnames(nearby_stns_final)[
      !(colnames(nearby_stns_final) %in% colnames(discharge_month_yr_data_all))])
# Select reference station metadata columns for join with discharge data
nearby_stns_simple <- nearby_stns_final[,..unique_metadata_col_names]

# Add METADATA columns to discharge data
discharge_month_yr_data_all <- discharge_month_yr_data_all[nearby_stns_simple, on = 'site_no']

# Summarize discharge data by river name, site number, and average discharge
decadal_monthly_q_rivers <- discharge_month_yr_data_all[year >= 1980][
  ,.(Q_cms_avg = mean(Q_cms)), 
  by = .(RiverName, ID, site_no)]

# Add river metadata to river name, site number, avg. discharge
decadal_monthly_q_wMeta <- merge(decadal_monthly_q_rivers, nearby_stns_final, by = c('RiverName', 'ID','site_no'))


#### 2B. CALCULATE MONTHLY SSC FOR EACH YEAR ####
river_landsat_pred <- merge(river_landsat_pred[,':='(year = year(sample_dt))],
                            decadal_monthly_q_rivers,
                            by = c('site_no'))

# Calculate monthly average SSC for each station, half-decade
river_landsat_annual <- river_landsat_pred[,.(
  SSC_mgL = mean(SSC_mgL, na.rm = T),
  SSC_mgL_sd = sd(SSC_mgL, na.rm = T),
  N_imgs = .N), # Std. error of SSC
  by = .(site_no, RiverName, drainage_area_km2, width, Latitude, Longitude, decade, year, month)][
    ,':='(SSC_mgL_se_rel = ifelse(is.na(SSC_mgL_sd), 0.73, 
                                  sqrt((SSC_mgL_sd/SSC_mgL)^2 + 0.73^2))/sqrt(N_imgs)) # Relative Std. error (includes 37% meas. error)
  ]


#### 2C. JOIN SSC WITH ANNUAL DISCHARGE DATA ####
# Join monthly SSC with monthly discharge, by river, decade, and month
# Keep all months from discharge record, which is complete 
# (every month has a value during satellite record).
# Some months have interpolated data, but all have data.
ssc_ann_monthly <- merge(river_landsat_annual, discharge_month_yr_data_all[year >= 1984], 
by = c('site_no','RiverName','month','decade','year'), all.y = T)[
  ,':='(tons_month = 0.0850354 * Q_cms * SSC_mgL * 30.5)
][,':='(tons_month_se = tons_month*SSC_mgL_se_rel)][
  ,.(site_no, RiverName, decade, month, year,
     SSC_mgL, SSC_mgL_sd, SSC_mgL_se_rel, Q_cms,Q_cms_period,
     tons_month, tons_month_se)
]

ssc_ann_monthly_n_obs <- ssc_ann_monthly[!is.na(SSC_mgL)][,.(N = .N),
                                                          by = .(site_no, RiverName)][
                                                            N > 75
                                                          ]
#### 3. REPLACE MISSING MONTHS USING DISCHARGE-QSS RATING CURVE ####
# For each site, replace missing values with rating-curve derived Qss
for(i in 1:nrow(ssc_ann_monthly_n_obs)){
  site_no_sel <- ssc_ann_monthly_n_obs[i]$site_no
  
  ssc_ann_monthly_sel <- ssc_ann_monthly[site_no == site_no_sel][
    Q_cms > 0
  ]
  
  # Plot timeseries of sediment flux
  # river_qss_annual_test_ts_plot <- ggplot(ssc_ann_monthly_sel,
  #                                         aes(x = year + month/12, y = tons_month/1e6)) +
  #   geom_line() +
  #   geom_point(size = 3, pch = 21, fill = NA) +
  #   geom_smooth(color = 'orange',se = F) +
  #   season_facet
  # 
  # print(river_qss_annual_test_ts_plot)
  
  # Number of observations per month
  months_N <- ssc_ann_monthly_sel[!is.na(tons_month)][
    ,.(n = .N), by = month
  ]
  # Select all the months with reliable data
  months_with_data <- months_N[n >= 10]$month
  
  
  ssc_ann_monthly_sel <- ssc_ann_monthly_sel[month %in% months_with_data]
  
  # Make a monthly Qss prediction function for each river, month, year
  # (?could make year a categorical variable?)
  qss_predict <- lm(
    log10(tons_month) ~ log10(Q_cms) + factor(month) + year
    # + I(year^2)
    ,
    data = ssc_ann_monthly_sel)
  # summary(qss_predict)
  
  ssc_ann_monthly_sel <- ssc_ann_monthly_sel[
    ,':='(tons_month_pred = 10^predict(qss_predict, newdata = ssc_ann_monthly_sel))][
      ,':='(tons_month = ifelse(!is.na(tons_month), tons_month, tons_month_pred))
    ]
  
  # Plot actual vs. prediction sediment flux
  # river_qss_annual_sel_regression_plot <- ggplot(ssc_ann_monthly_sel,
  #                                         aes(x = tons_month/1e6, y = tons_month_pred/1e6)) +
  #   geom_point(size = 3, pch = 21, fill = NA) +
  #   geom_abline(slope = 1, intercept = 0) +
  #   geom_smooth(method = 'lm', color = 'orange',se = F) +
  #   season_facet
  # 
  # print(river_qss_annual_sel_regression_plot)
  # 
  # Plot timeseries of actual and predicted sediment flux
  # river_qss_annual_sel_ts_plot <- ggplot(ssc_ann_monthly_sel,
  #                                         aes(x = year + month/12, y = tons_month/1e6)) +
  #   geom_line() +
  #   geom_point(size = 3, pch = 21, fill = NA) +
  #   geom_line(aes(y = tons_month_pred/1e6), color = 'blue') +
  #   geom_point(aes(y = tons_month_pred/1e6), color = 'blue') +
  #   geom_smooth(color = 'orange',se = F) +
  #   season_facet
  # 
  # print(river_qss_annual_sel_ts_plot)
  
  if(i == 1){
    ssc_ann_filled <- ssc_ann_monthly_sel
  }else{
    ssc_ann_filled <- rbind(ssc_ann_filled, ssc_ann_monthly_sel,
                            use.names = T, fill = T)
  }
}

# Formula for SSC based on tons/month
# SSC_mgL = tons_month/(Q_cms * 0.0850354 *  30.5)

ssc_ann_filled <- ssc_ann_filled[,':='(
  SSC_mgL = ifelse(!is.na(SSC_mgL), SSC_mgL, tons_month/(Q_cms * 0.0850354 *  30.5)),
  SSC_source = ifelse(!is.na(SSC_mgL), 'Landsat observation','Rating curve estimate')
  )]


#### 4. WRITE DATA TO DRIVE ####
# These files will be used in subsequent analysis steps
# Save monthly suspended sediment transport data
fwrite(ssc_ann_filled, file = paste0(wd_imports, 'sediment_transport_month_yr_data.csv'))
# Save metadata
fwrite(decadal_monthly_q_wMeta, file = paste0(wd_imports, 'discharge_and_reference_site_metadata.csv'))

# Validation data
validation_data_monthly <- ssc_ann_filled[RiverName %in% c('Mississippi','Changjiang','Huaihe','Zhujiang')]
fwrite(validation_data_monthly, file = paste0(wd_imports, 'satellite_validation_data_monthly.csv'))
#### TESTING ####
# ggplot(ssc_ann_filled[RiverName == 'Zhujiang'], aes(x = year + month/12, y = SSC_mgL, color = SSC_source)) +
#   geom_point() +
#   geom_line(aes(group = RiverName))
# 
# ggplot(ssc_ann_filled[RiverName == 'Zhujiang'], aes(x = Q_cms, y = tons_month, color = SSC_source)) +
#   geom_point() +
#   scale_x_log10(labels = fancy_scientific_modified) +
#   scale_y_log10(labels = fancy_scientific_modified) +
#   geom_smooth(method = 'lm', aes(group = decade, color = factor(decade)))


