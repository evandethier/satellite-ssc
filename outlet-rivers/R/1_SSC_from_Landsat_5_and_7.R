# Import Landsat data (raw output from Google Earth Engine)
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
#### 1. IMPORT LANDSAT DATA AND DEFINE COLUMNS ####
# Import Landsat data from river outlets
river_landsat_import <- fread(paste0(wd_imports, 'global_outlet_rivers_120km_river_training_ls57_rawBands_b7lt500_15km2.csv'))

# PREPARE DATA COLUMNS, MAKE LANDSAT MATCHUP, CALCULATE SSC ESTIMATE
# Landsat data do have station information
# They also have latitude and longitude
river_landsat_import <- na.omit(river_landsat_import[,
                                  ':='(
                                      site_no = paste0('st_', site_no),
                                       # Rename columns for simplicity
                                       B1 = B1_median,
                                       B2 = B2_median,
                                       B3 = B3_median,
                                       B4 = B4_median,
                                       B5 = B5_median,
                                       B6 = B6_median,
                                       B7 = B7_median,
                                       num_pix = B2_count,
                                       sample_dt = ymd(date),
                                       landsat_dt = ymd(date)
                                  )]
                             , cols = c('B1','B2','B3','B4','B5','B7'))[
                               B1 > 0 & B2 > 0 & B3 > 0 & B4 > 0 & B5 > 0 & B7 > 0 &
                               B1 < 5000 & B2 < 5000 & B3 < 5000 & B4 < 5000 & B6 < 4000][
                                 ,':='( 
                                   # add new columns with band ratios
                                   B1.2 = B1^2,
                                   B2.2 = B2^2,
                                   B3.2 = B3^2,
                                   B4.2 = B4^2,
                                   B5.2 = B5^2,
                                   B7.2 = B7^2,
                                   B2.B1 = B2/B1,
                                   B3.B1 = B3/B1,
                                   B4.B1 = B4/B1,
                                   B5.B1 = B5/B1,
                                   B7.B1 = B7/B1,
                                   B3.B2 = B3/B2,
                                   B4.B2 = B4/B2,
                                   B5.B2 = B5/B2,
                                   B7.B2 = B7/B2,
                                   B4.B3 = B4/B3,
                                   B5.B3 = B5/B3,
                                   B7.B3 = B7/B3,
                                   B5.B4 = B5/B4,
                                   B7.B4 = B7/B4,
                                   B7.B5 = B7/B5,
                                   Latitude = gage_lat,
                                   Longitude = gage_long
                                   )][ 
                                     # select only columns of interest
                                     ,.(site_no, station_nm, 
                                        width, drainage_area_km2,
                                        Latitude,Longitude,sample_dt, num_pix, snow_ice_qa_count, landsat_dt,
                                        B1,B2,B3,B4,B5,B6,B7,B2.B1,B3.B1,B4.B1,B5.B1,B7.B1,B3.B2,B4.B2,B5.B2,
                                        B7.B2,B4.B3,B5.B3,B7.B3,B5.B4,B7.B4,B7.B5, B1.2,B2.2,B3.2,B4.2,B5.2,B7.2
                                     )][site_no != "0"][
# Filter images with too much snow/ice interference through various decision trees
                                       !((B6 < 2800 & B1 > 900 & B2 > 900 & B3 > 900 & B5 > 300 & B7 > 200 & B1 > B3 & B1 < B4) | # Elimate snowy & cold images
    (B1 > 700 & B1/B2 > 1.2 & B5 > 200)|
    ((B1 + B2 + B3 + B4) > 3200 & B3 < B1 & B3/B1 < 1.5 & B6 < 2750 & B5 > 300) |
    (B4 > 1500 & B4/B3 > 1.5 & B6 < 2800)| # This eliminates many cloudy/snowy images at high latitudes
  ((B1 + B2 + B3 + B4) > 4000 & B6 < 2750 & B5 > 300 & abs(Latitude) > 40)|
    (B1 > 700 &
       snow_ice_qa_count > (num_pix * 10) & 
       snow_ice_qa_count > 500 &
       B3/B1 < 1.5))
]

# Generate metadata for each river
# Lat/Long, drainage area, number of images per river, avg. number of pixels per river
outlet_rivers_ls_metadata <- river_landsat_import[,.(n_imgs = uniqueN(sample_dt),
                                     n_pixels = mean(num_pix, na.rm = T)),
                                  by = .(site_no, Latitude, Longitude, drainage_area_km2, width)][
                                    width > 30
                                  ]
                                  
#### 2. ASSIGN CLUSTER TO ALL DATA ####
# Import cluster center object (for n clusters = 1 to 10)
# This object was generated in the initial SSC calibration from Dethier, 2020 JGR-ES
clusters_calculated_list <- readRDS(paste0(wd_imports, 'cluster_centers.rds'))

# This approach assigns a cluster to each river
# Other approaches allow for more dynamic cluster assignment
# e.g., by month/season or decade
getCluster <- function(df,clustering_vars,n_centers, kmeans_object){
  # Compute band median at each site for clustering variables
  site_band_quantiles_all <- df[
    ,':='(month = month(sample_dt))][
      # n_insitu_samples_bySite][N_insitu_samples > 12
      ,.(N_samples = .N,
         B1 = median(B1),
         B2 = median(B2),
         B3 = median(B3),
         B4 = median(B4),
         # B5 = median(B5),
         # B7 = median(B7),
         B2.B1 = median(B2.B1),
         B3.B1 = median(B3.B1),
         B4.B1 = median(B4.B1),
         B3.B2 = median(B3.B2),
         B4.B2 = median(B4.B2),
         B4.B3 = median(B4.B3),
         B4.B3.B1 = median(B4.B3/B1)), 
      by = .(station_nm,site_no)]
  
  site_band_quantile_scaled <- scale(site_band_quantiles_all[,..clustering_vars], 
                                     center = attributes(site_band_scaling)$`scaled:center`, 
                                     scale = attributes(site_band_scaling)$`scaled:scale`)
  
  closest.cluster <- function(x) {
    cluster.dist <- apply(kmeans_object$centers, 1, function(y) sqrt(sum((x-y)^2)))
    return(which.min(cluster.dist)[1])
  }
  site_band_quantiles_all$cluster <- apply(site_band_quantile_scaled, 1, closest.cluster)
  
  df_cluster <- merge(df,site_band_quantiles_all[,c('site_no','cluster')], by = 'site_no')
  df_cluster$cluster_sel <- df_cluster$cluster
  return(df_cluster)
  
}

# (NOT USED IN THIS APPLICATION. Uncomment in below line to choose monthly assignment)
# This approach assigns a cluster to each river, !!BY MONTH!!
# Allows for a dynamic cluster assignment
# e.g., by month/season or decade
# Function to create river clustering function: break out by month
getCluster_monthly <- function(df,clustering_vars,n_centers, kmeans_object){
  # Compute band median at each site for clustering variables
  site_band_quantiles_all <- df[
    ,':='(month = month(sample_dt))][
      # n_insitu_samples_bySite][N_insitu_samples > 12
      ,.(N_samples = .N,
         B1 = median(B1),
         B2 = median(B2),
         B3 = median(B3),
         B4 = median(B4),
         # B5 = median(B5),
         # B7 = median(B7),
         B2.B1 = median(B2.B1),
         B3.B1 = median(B3.B1),
         B4.B1 = median(B4.B1),
         B3.B2 = median(B3.B2),
         B4.B2 = median(B4.B2),
         B4.B3 = median(B4.B3),
         B4.B3.B1 = median(B4.B3/B1)), 
      by = .(station_nm,site_no, month)]
  
  site_band_quantile_scaled <- scale(site_band_quantiles_all[,..clustering_vars], 
                                     center = attributes(site_band_scaling)$`scaled:center`, 
                                     scale = attributes(site_band_scaling)$`scaled:scale`)
  
  closest.cluster <- function(x) {
    cluster.dist <- apply(kmeans_object$centers, 1, function(y) sqrt(sum((x-y)^2)))
    return(which.min(cluster.dist)[1])
  }
  site_band_quantiles_all$cluster <- apply(site_band_quantile_scaled, 1, closest.cluster)
  
  df_cluster <- merge(df,site_band_quantiles_all[,c('site_no','cluster', 'month')], by = c('site_no', 'month'))
  df_cluster$cluster_sel <- df_cluster$cluster
  return(df_cluster)
  
}
# Get cluster for each site based on typical spectral profile
# Apply cluster to every Landsat daily observation
# For this application, use single cluster for river
# Other applications may benefit from a monthly breakdown
# river_landsat_cl <- getCluster_monthly(
river_landsat_cl <- getCluster(
  river_landsat_import, 
  clustering_vars,cluster_n_best, 
  clusters_calculated_list[[cluster_n_best]])

# Calculate monthly median color for each river
site_monthly_band_median_color <- river_landsat_cl[
    ,.(N_samples = .N,
       B1 = median(B1),
       B2 = median(B2),
       B3 = median(B3),
       B4 = median(B4),
       B5 = median(B5),
       B7 = median(B7)
    ),
  by = .(site_no, month)]
  
# Save cluster centers (CAREFUL! SEE BELOW)
# Commented out to avoid overwriting. 
# Do not comment back in unless you are absolutely sure you want to overwrite!
# saveRDS(clusters_calculated_list, paste0(wd_imports, 'cluster_centers.rds'))


#### 3. CALCULATE SSC ####
# Import SSC prediction functions (for n clusters = 1 to 10)
ssc_cluster_funs <- readRDS(paste0(wd_imports,'SSC_cluster_function.rds'))

# Save SSC estimation algorithm. (CAREFUL! SEE BELOW)
# Commented out to avoid overwriting. 
# Do not comment back in unless you are absolutely sure you want to overwrite!
# saveRDS(ssc_cluster_funs, file = paste0(wd_imports,'SSC_cluster_function.rds'))

# Function to apply SSC calibration models to make predictions based on new surface reflectance inputs (cluster needed)
# For base, cluster, and site predictions
getSSC_pred <- function(lm_data, regressors, cluster_funs){ # Version that includes site specification
  lm_data$pred_st <- NA
  lm_data[,ssc_subset:=cluster_sel] # clusters
  subsets <- unique(lm_data$ssc_subset)
  for(i in subsets){ # for individual cluster models
    # print(i)
    regressors_sel <- regressors[-which(regressors == 'site_no')]
    lm_data_lm <- lm_data[ssc_subset == i] # only chooses sites within cluster
    
    ssc_lm <- cluster_funs[[i]]
    glm_pred <- predict(ssc_lm, newx = as.matrix(lm_data_lm[,..regressors_sel]), s = "lambda.1se")
    lm_data[ssc_subset == i, pred_cl:= glm_pred]
    lm_data_lm <- NA
    # lm_data$res_cl[which(lm_data$ssc_subset == i)] <- resid(ssc_lm)
  }
  return(lm_data)
}

# Run SSC prediction algorithm to get clustered prediction for SSC
river_landsat_pred <- getSSC_pred(na.omit(river_landsat_cl, cols = c(regressors_all, 'cluster_sel')), 
                                  regressors_all, ssc_cluster_funs)[,':='(
                                    SSC_mgL = 10^pred_cl,
                                    month = month(sample_dt),
                                    # decade = ifelse(year(sample_dt) < 1990, 1990,
                                    #                 ifelse(year(sample_dt) > 2019, 2010, 
                                    #                        year(sample_dt) - year(sample_dt)%%10)))]
                                    decade = ifelse(year(sample_dt) < 1990, 1990,
                                                    ifelse(year(sample_dt) > 2019, 2015, 
                                                           year(sample_dt) - year(sample_dt)%%5)))][
                                                             SSC_mgL < 20000
                                                           ]


# Select only necessary columns
river_landsat_pred <- river_landsat_pred[,.(site_no, sample_dt, month, decade, SSC_mgL, 
                                            cluster, Latitude, Longitude, 
                                            drainage_area_km2, num_pix, width)]

# Address known issue with Congo River
# Based on HYBAM in situ measurements
river_landsat_pred <- river_landsat_pred[site_no == 'st_000000000000000013b9', ':='(SSC_mgL = SSC_mgL/8)]

#### CALCULATE BAND AVERAGES FOR COLOR VISUALIZATION ####                                                 
# Get median monthly cluster assignment for visualization
river_cluster_assignment_monthly <- river_landsat_pred[
  ,.(cluster = median(cluster, na.rm = T)),
  by = .(site_no, Latitude, Longitude, month)]

# Get river color for each cluster and SSC category
cluster_ssc_cat_colors <- merge(
      river_landsat_pred,
      river_landsat_cl[,.(site_no, sample_dt, B1, B2, B3)],
      by = c('site_no','sample_dt'))[,':='(
        SSC_cat = cut(SSC_mgL,
                      breaks = c(-Inf, 20, 50, 100, 300, 500, Inf),
                      labels = c('0-20', '21-50', '51-100', '100-300', '300-500', '>500')))][
                        ,.(blue = mean(B1, na.rm = T),
                           green = mean(B2, na.rm = T),
                           red = mean(B3, na.rm = T)),
                        by = .(SSC_cat, cluster)
                        ]


#### 4. EXPORT DATA FOR SUBSEQUENT STEPS ####
# Landsat-derived daily SSC estimates
fwrite(river_landsat_pred, file = paste0(wd_imports, 'river_landsat_pred.csv'))

# Cluster assignments for each river (by month)
fwrite(river_cluster_assignment_monthly, file = paste0(wd_imports, 'river_cluster_assignment_monthly.csv'))

# River metadata, including site ID, Lat/Long, drainage area, width, number of images, number of pixels
fwrite(outlet_rivers_ls_metadata, file = paste0(wd_imports, 'outlet_rivers_ls_metadata.csv'))

# River color, by month
fwrite(site_monthly_band_median_color, file = paste0(wd_imports, 'site_monthly_band_median_color.csv'))

# River color for each cluster and SSC category
fwrite(cluster_ssc_cat_colors, file = paste0(wd_imports, 'cluster_ssc_cat_colors.csv'))



#### 5. LANDSAT SSC TESTING AND INSPECTION ####
## Test and visualize clusters
# Plot band 1 (blue) and band 3 (red) to visualize cluster groups
blue_vs_red_plot_by_cluster <- ggplot(river_landsat_cl[year(sample_dt) == 2011 & 
                                                         B1 > 0 & B1 < 2500 &
                                                         B3 > 0 & B3 < 2500
], 
aes(x = B1, y = B3, color = factor(cluster))) +
  geom_point(alpha = 0.1) +
  # geom_density_2d(bins = 10) +
  season_facet +
  # scale_x_continuous(limits = c(0, 2500)) +
  # scale_y_continuous(limits = c(0, 2500)) +
  # facet_wrap(.~factor(cluster)) +
  labs(
    x = 'Band 1 (blue)',
    y = 'Band 3 (red)'
  )


#### 6. VISUALIZE LANDSAT IMAGE STATISTICS ####
# How many images are there, images per river, pixels per image, etc.

## CALCULATION FOR MANUSCRIPT ##
# N Images per site
summary(outlet_rivers_ls_metadata)
## Fig. S1A
# Average discharge-weighted SSC at each station
ls_stations_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = outlet_rivers_ls_metadata, pch = 21, stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = n_imgs,
                 size = drainage_area_km2/1000)) +
  scale_fill_gradientn(trans = 'log10', colors = c('black','blue','grey','red')) +
  scale_size_continuous(range = c(2,5), trans = 'sqrt', breaks = c(50, 100, 500, 1000),
                        guide = guide_legend(override.aes = list(fill = 'black'))) +
  season_facet +
  scale_x_continuous(limits = c(-183, 195), expand = expansion(mult = c(0, 0)),
                     breaks = seq(-180,180,60),
                     labels = long_dd_lab,
                     sec.axis = dup_axis()) +
  scale_y_continuous(limits = c(-61, 86), expand = expansion(mult = c(0, 0)),
                     breaks = seq(-60,80,20),
                     labels = lat_dd_lab,
                     sec.axis = dup_axis()) + 
  labs(
    x = '',
    y = '',
    fill = '**N images**',
    size = '**Drainage Area**<br>**(1000s km<sup>2</sup>)**'
  ) +
  theme(legend.position = c(0.01, 0.39),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(ls_stations_map, filename = paste0(wd_figures, 'ls_stations_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(ls_stations_map, filename = paste0(wd_figures, 'ls_stations_map.png'),
       width = 11, height = 6)

## Fig. S1B
# Plot images per site
images_per_site_plot <- ggplot(outlet_rivers_ls_metadata, aes(x = n_imgs)) +
  # c('#CABFD9','#565902','#BF9004')
  geom_histogram(binwidth = 50, fill = '#368096', color = 'black', lwd = 0.2) +
  season_facet + 
  scale_x_continuous(expand = expansion(mult = c(0,0))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  theme(panel.background = element_rect(fill = 'grey90'),
        panel.grid.major.y = element_line(color = 'white', size = 0.25),
        panel.grid.minor.y = element_line(color = 'white', size = 0.25)) +
  labs(x = 'N images per river',
       y = 'Count')

## Fig. S1C
# Plot pixels per sample
pixels_per_sample_plot <- ggplot(outlet_rivers_ls_metadata, aes(x = n_pixels)) +
  # c('#CABFD9','#565902','#BF9004')
  geom_histogram(bins = 30, fill = '#D67258', color = 'black', lwd = 0.2) +
  season_facet + 
  # scale_x_continuous(expand = expansion(mult = c(0,0))) +
  scale_x_log10(labels = fancy_scientific_modified) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  theme(panel.background = element_rect(fill = 'grey90'),
        panel.grid.major.y = element_line(color = 'white', size = 0.25),
        panel.grid.minor.y = element_line(color = 'white', size = 0.25)) +
  labs(x = 'N pixels per river sample',
       y = 'Count')

## Fig. S1D
# Plot number of sample pixels vs. width
pixels_vs_width_plot <- ggplot(outlet_rivers_ls_metadata, aes(x = width, y = n_pixels)) +
  # c('#CABFD9','#565902','#BF9004')
  geom_point() +
  season_facet + 
  scale_x_log10(labels = fancy_scientific_modified) +
  scale_y_log10(labels = fancy_scientific_modified) +
  theme(panel.background = element_rect(fill = 'grey90'),
        panel.grid.major.y = element_line(color = 'white', size = 0.25),
        panel.grid.minor.y = element_line(color = 'white', size = 0.25)) +
  labs(y = 'N pixels per river sample',
       x = 'Approximate river width (m)')

## Fig. S1B-D
# Combine image statistics plots
image_statistics_panel_plot <- 
  ggpubr::ggarrange(images_per_site_plot, pixels_per_sample_plot, pixels_vs_width_plot,
            ncol = 3, align = 'hv', labels = c('B','C','D'))

# Save combined image statistics plots
# (Just a part of what is used in the paper)
ggsave(image_statistics_panel_plot, filename = paste0(wd_figures, 'image_statistics_panel_plot.pdf'),
       width = 8, height = 4, useDingbats = F)
ggsave(image_statistics_panel_plot, filename = paste0(wd_figures, 'image_statistics_panel_plot.png'),
       width = 8, height = 4)

## Fig. S1
# Combine map of N sites with image statistics
image_statistics_panel_plot_wmap <- ggpubr::ggarrange(
  ls_stations_map, NULL, image_statistics_panel_plot,
  heights = c(1, 0, 0.7), ncol = 1, labels = c('A','',''), align = 'v')

ggsave(image_statistics_panel_plot_wmap, filename = paste0(wd_figures, 'FigS1_image_statistics_panel_plot_wmap.pdf'),
       width = 11, height = 9, useDingbats = F)
ggsave(image_statistics_panel_plot_wmap, filename = paste0(wd_figures, 'FigS1_image_statistics_panel_plot_wmap.png'),
       width = 11, height = 9)

