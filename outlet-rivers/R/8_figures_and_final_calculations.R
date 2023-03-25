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
############## I. IMPORT DATA ############## 
#### 1A. SUSPENDED SEDIMENT TRANSPORT DATA AND METADATA ####
# Import monthly suspended sediment transport data
sed_transport_month_yr <- fread(paste0(wd_imports, 'sediment_transport_month_yr_data.csv'))
# Import discharge and reference site metadata
Q_stn_metadata <- fread(paste0(wd_imports, 'discharge_and_reference_site_metadata.csv'),
                        colClasses = c('ID' = 'character'))

# Merge SSC data with metadata
sed_transport_month_yr <- merge(
  sed_transport_month_yr,
  Q_stn_metadata, 
  by = c('RiverName','site_no'), 
  all.y = T)

# Cluster assignments for each river (by month)
river_cluster_assignment_monthly <- fread(file = paste0(wd_imports, 'river_cluster_assignment_monthly.csv'))
river_cluster_assignment_single <- river_cluster_assignment_monthly[,.(cluster = get_mode(cluster)), 
                                                                    by = .(site_no, Latitude, Longitude)]
# River color, by month
site_monthly_band_median_color <- fread(paste0(wd_imports, 'site_monthly_band_median_color.csv'))

# River color for each cluster and SSC category
cluster_ssc_cat_colors <- fread(file = paste0(wd_imports, 'cluster_ssc_cat_colors.csv'))
cluster_ssc_cat_colors <- cluster_ssc_cat_colors[
  ,':='(SSC_cat = factor(SSC_cat, 
                         levels = c("0-20", "21-50", "51-100", "100-300", "300-500",">500")))]

#### 1B. SUMMARY DATA OF SUSPENDED SEDIMENT TRANSPORT, BY-RIVER ####
# Just name, site_no, and id
ssc_rivers <- fread(file = paste0(wd_imports, 'ssc_outlet_rivers_names.csv'))

# Add discharge and area metadata
ssc_rivers_stats <- fread(file = paste0(wd_imports, 'ssc_rivers_stats.csv'))

# Full by-river statistics and metadata
ssc_rivers_stats_full <- fread(file = paste0(wd_imports, 'ssc_rivers_stats_full.csv'))

# Pre-dam estimates from Milliman and Farnsworth, 2012
Qss_pre_byRiver <- fread(file = paste0(wd_imports, 'Qss_pre_byRiver.csv'))

#### 1C. TIMESERIES OF SUSPENDED SEDIMENT TRANSPORT, BY-RIVER ####
# Annual by-river flux with metadata
Qss_river_annual <- fread(file = paste0(wd_imports, 'Qss_river_annual.csv'))

# Annual by-river flux with metadata
Qss_river_decadal <- fread(file = paste0(wd_imports, 'Qss_river_decadal.csv'))

# Annual by-river SSC with metadata. 
# Only direct observations and no discharge, so weighted towards low SSC
SSC_annual_avg <- fread(file = paste0(wd_imports, 'SSC_annual_avg_direct_observations_no_discharge.csv'))

#### 1D. MONTHLY SUSPENDED SEDIMENT TRANSPORT SUMMARIES, BY-RIVER ####
# Monthly averages of discharge, sediment flux
avg_monthly_Q_Qss <- fread(file = paste0(wd_imports, 'avg_monthly_Q_Qss.csv'))

# Timing of peak suspended sediment concentration
Qss_peak_month <- fread(file = paste0(wd_imports, 'Qss_peak_month.csv'))

#### 1E. SUSPENDED SEDIMENT TRANSPORT TIMESERIES AND SUMMARIES, BY-CONTINENT ####
# Drainage area and total discharge, by continent
continent_area_Q <- fread(file = paste0(wd_imports, 'continent_area_Q.csv'))

# Annual sediment flux, by continent
Qss_global_annual_ref <- fread(file = paste0(wd_imports, 'Qss_global_annual_ref.csv'))

#### 1F. DAM-BUILDING TIMESERIES DATA ####
# Timeseries of dam building, by river
dam_timeseries <- fread(file = paste0(wd_imports, 'dam_timeseries.csv'))

#### 1G. DEFORESTATION DATA ####
watershed_deforestation <- fread(paste0(wd_imports, 'watershed_deforestation.csv'))
#### 1H. WATERSHED CLIMATE ####
watershed_climate <- fread(paste0(wd_imports, 'outlet-rivers-watershed-climate-stats.csv'))[
  ,':='(
    precip_m_yr = total_precipitation*12,
    site_no = paste0('st_',site_no))
][
  # Join with general watershed statistics
  # Could do annual statistics for more detail
  ssc_rivers_stats[,.(site_no, RiverName, Latitude, Longitude, Q_cms_avg)], on = 'site_no']


############## II. ANALYZE CHANGES ############## 
#### 2. MAP SSC AVG, SEASONALITY FOR ALL RIVERS ####
## Fig. 1B
# Average discharge-weighted SSC at each station
avg_ssc_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = ssc_rivers_stats_full, pch = 21, stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = SSC_mgL,
                 size = Q_cms_avg/1000)) +
  scale_fill_gradientn(trans = 'log10', colors = c('black','#286098','#529AE3','grey80','#966D27','#966D27'),
                       guide = guide_colorbar(order = 1)) +
  scale_size_continuous(range = c(2,10), breaks = c(1, 10, 100)) +
  guides(size = guide_legend(override.aes = list(fill = "black"))) +
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
    fill = '**SSC (mg/L)**',
    size = '**Avg. Discharge**<br>**(1,000s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.01, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(avg_ssc_map, filename = paste0(wd_figures, 'avg_ssc_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(avg_ssc_map, filename = paste0(wd_figures, 'avg_ssc_map.png'),
       width = 11, height = 6)


## Fig. S4
# Plot month of peak SSC
peak_ssc_month_map <- ggplot(data = Qss_peak_month) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_polygon(aes(x = 0, y = 0, fill = factor(month))) +
  geom_point(pch = 21, stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = factor(month),
                 # fill = month,
                 size = Q_cms/1000), show.legend = F) +
  scale_fill_discrete(
    breaks = c(1:12),
    # colors = rev(c('#FE2712','#FC600A','#FB9902','#FCCC1A','#FEFE33','#B2D732',
    #            '#66B032','#347C98','#0247FE','#4424D6','#8601AF','#C21460')),
    guide = guide_legend(
      keyheight = unit(1.25, units = "mm"), 
      keywidth=unit(7, units = "mm"), 
      override.aes = list(alpha = 1,
                          size = 2),
      label.position = "bottom", 
      title.position = 'top', nrow=1)) +
  # scale_fill_gradientn(trans = 'log10', colors = c('black','#286098','#529AE3','grey80','#966D27','#966D27')) +
  # scale_fill_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  #                              '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')) +
  scale_size_continuous(range = c(2,10), breaks = c(1, 10, 100)) +
  guides(size = guide_legend(override.aes = list(fill = "black"))) +
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
    fill = 'Month of Peak Sediment Flux',
    size = '**Avg. Discharge**<br>**(100s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.42, 0.08),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(peak_ssc_month_map, filename = paste0(wd_figures, 'FigS4_peak_ssc_month_map.png'),
       width = 11, height = 6)
ggsave(peak_ssc_month_map, filename = paste0(wd_figures, 'FigS4_peak_ssc_month_map.pdf'),
       width = 11, height = 6, useDingbats = F)

 # Plot average ssc as a function of latitude
# Color by climate
# (Not used in paper)
avg_ssc_by_latitude_climate_plot <- ggplot(ssc_rivers_stats_full,
       aes(x = Latitude, y = SSC_mgL)) +
  geom_point(color = 'black', pch = 21, stroke = 0.2, 
             aes(fill = Climate_T, size = Q_cms_avg/100)) +
  geom_smooth(color = 'black', se = F) +
  season_facet +
  scale_y_log10(limits = c(1, 10000), labels = fancy_scientific_modified) +
  scale_size_continuous(range = c(2,10), breaks = c(100, 500, 1000)) +
  rotate()

#### 3. LARGEST 20 RIVERS -- ANALYSIS ####
largest_20_rivers <- ssc_rivers_stats_full[order(Q_cms_avg, decreasing = T)][1:21]$RiverName

# Remove Congo because of known issue. This issue has been corrected, but still avoiding it.
sed_transport_month_yr_largest_20_rivers <- sed_transport_month_yr[RiverName %chin% largest_20_rivers & RiverName != 'Congo']

## Fig. S3
# Rating curve between monthly discharge and monthly sediment flux
# for the 20 largest rivers
# Observed and interpolated data points are shown
large_river_Q_Qss_rating_curves_plot <- ggplot(sed_transport_month_yr_largest_20_rivers, 
                                               aes(x = Q_cms/1000, y = tons_month_pred/1e6)) +
  geom_point(aes(y = tons_month/1e6, color = 'Satellite')) +
  geom_point(aes(color = 'Regression')) + 
  scale_color_manual(values = c('Regression' = 'red', 'Satellite' = 'blue')) +
  scale_x_log10(guide = guide_axis(check.overlap = TRUE)) +
  scale_y_log10() +
  season_facet +
  facet_wrap(.~ifelse(RiverName == 'Changjiang','Yangtze', RiverName), scales = 'free', ncol = 4) +
  theme(axis.title.x = element_markdown()) +
  labs(
    x = 'Avg. Monthly Discharge (1,000s m<sup>3</sup>/s)',
    y = 'Monthly sediment flux (Mt)'
  )

ggsave(large_river_Q_Qss_rating_curves_plot, filename = paste0(wd_figures, 'FigS3_large_river_Q_Qss_rating_curves_plot.pdf'),
       width = 9, height = 14, useDingbats = F)
ggsave(large_river_Q_Qss_rating_curves_plot, filename = paste0(wd_figures, 'FigS3_large_river_Q_Qss_rating_curves_plot.png'),
       width = 9, height = 14)



#### 4. DAM BUILDING ####
#### 4A. DAM RESERVOIR STORAGE TIMESERIES BY RIVER, CONTINENT ####
# Calculate continent area and annual discharge
# (TO DO: calculate annual precipitation here too)

# Dam area and storage, cumulative by continent
dam_timeseries_continent <- dam_timeseries[Continent_Region != ''][
  ,.(area_km2_cml = sum(area_km2_cml, na.rm = T),
     storage_mcm_cml = sum(cap_mcm_cml, na.rm = T)),
  by = .(year, Continent_Region)][
    continent_area_Q, on = 'Continent_Region']

dam_building_landsat_era <- dcast.data.table(Continent_Region ~ year,
                                             value.var = c('area_km2_cml','storage_mcm_cml'),
                                             data = dam_timeseries_continent[year %in% c(1984, 2020)])

# Combine dam building and suspended sediment transport timeseries
dam_ssc_comb <- sed_transport_month_yr[
  dam_timeseries[,
                 .(RiverName, site_no, year, cap_mcm_cml, area_km2_cml)],
  on = c('RiverName', 'site_no', 'year')][
    ,':='(cap_mcm_cml = ifelse(is.na(cap_mcm_cml), 0, cap_mcm_cml))
  ]

#### 4B. DAM IMPACTS ON A PER-RIVER BASIS ####
# Calculate dams per river in 1984
dams_per_river_pre1984 <- na.omit(dam_ssc_comb[year < 1984,
                                               .(cap_mcm_cml_pre1984 = max(cap_mcm_cml, na.rm = T)),
                                               by = .(RiverName, site_no)],
                                  cols = c('cap_mcm_cml_pre1984'))

# Dams per river and sediment impacts
dams_per_river <- na.omit(dam_ssc_comb[,.(tons_month_decadal = mean(tons_month, na.rm = T),
                                          tons_month_decadal_sd = sd(tons_month, na.rm = T),
                                          cap_mcm_cml = max(cap_mcm_cml, na.rm = T)),
                                       by = .(RiverName, site_no, Continent_Region, Country, Q_cms_avg, drainage_area_km2, Latitude, Longitude)],
                          cols = c('Q_cms_avg','tons_month_decadal', 'cap_mcm_cml'))[
                            ,':='(RCI = round(cap_mcm_cml/(Q_cms_avg * (3600 * 24 * 365)/(1000^3)) * 0.001 * 100,2))
                          ][dams_per_river_pre1984, on = c('RiverName','site_no')]

# Add categories for dam timing based on 1984 Landsat period cutoff
dams_per_river <- dams_per_river[,':='(
  when_dams = ifelse(cap_mcm_cml == 0, 'No Dams',
                     ifelse((cap_mcm_cml - cap_mcm_cml_pre1984) > 0 & cap_mcm_cml_pre1984 == 0, 'Only Post-1984 Dams',
                            ifelse((cap_mcm_cml - cap_mcm_cml_pre1984) > 0 & cap_mcm_cml_pre1984 > 0, 'Pre- and Post-1984 Dams',
                                   'Only Pre-1984 Dams'))))][!is.na(when_dams)]


fwrite(dams_per_river, file = paste0(wd_exports, 'dams_per_river.csv'))

## CALCULATION FOR PAPER ##
# When were dams built for each river, pre- and post-1984 analysis
dam_timing_global <- dams_per_river[,.(.N), by = when_dams]
rivers_with_major_dams_post1984 <- dams_per_river[(cap_mcm_cml - cap_mcm_cml_pre1984) > 0]
rivers_with_major_dams <- dams_per_river[cap_mcm_cml > 0]

## CALCULATION FOR PAPER ##
# (Not used)
# What is the RCI for China?
china_rci <- dams_per_river[Country == 'China', .(RCI = round(
  sum(cap_mcm_cml, na.rm = T)/
    (sum(Q_cms_avg, na.rm = T) * 
       (3600 * 24 * 365)/(1000^3)) * 0.001 * 100,2)
)]
russia_rci <- dams_per_river[Country == 'Russia' & Continent_Region == 'Asia', .(RCI = round(
  sum(cap_mcm_cml, na.rm = T)/
    (sum(Q_cms_avg, na.rm = T) * 
       (3600 * 24 * 365)/(1000^3)) * 0.001 * 100,2)
)]
india_rci <- dams_per_river[Country == 'India' & Continent_Region == 'Asia', .(RCI = round(
  sum(cap_mcm_cml, na.rm = T)/
    (sum(Q_cms_avg, na.rm = T) * 
       (3600 * 24 * 365)/(1000^3)) * 0.001 * 100,2)
)]

#### 5. CALCULATE CHANGE IN SEDIMENT TRANSPORT ####
#### 5A. CHANGE IN FLUX 1980S vs. 2020 ####
# Change in Qss, beginning of Landsat record to present decade 
# (Calculated using half-decade annual avg. Qss), from avg. SSC * monthly Q)
Qss_river_decadal_change <- dcast.data.table(
  Qss_river_decadal,
  RiverName + site_no + ID + Latitude + Longitude ~ decade,
  value.var = 'tons_month_avg')[
    ,':='(tons_month_2015_1990 = ifelse(!is.na(`1990`), (`2015` - `1990`)/`1990`,
                                        (`2015` - `2000`)/`2000`)*100)
  ][
    # If the difference between epochs is significant, get the sign of the change
    ,':='(Qss_sign = ifelse(abs(tons_month_2015_1990) < 25, 'No change', 
                            ifelse(tons_month_2015_1990 > 0, 'Increase','Decrease')))]


#### 5B. KENDALL TREND IN FLUX AND SSC FOR EACH RIVER ####
# Only make plots if specified: significantly speeds up process if no plots 
# (< 15 seconds vs. 22 minutes)
# make_ts_plots <- 'Yes'
make_ts_plots <- 'No'

for(i in 1:length(unique(Qss_river_annual$site_no))){
  
  site_no_sel <- unique(Qss_river_annual$site_no)[i]
  dt_sel <- Qss_river_annual[site_no == site_no_sel]
  location_sel <- dt_sel[1, .(Latitude, Longitude)]
  n_years <- nrow(dt_sel[tons_month_avg > 0])
  if(n_years > 10){
    
    trend_test <- Kendall(dt_sel$year, dt_sel$tons_month_avg)
    trend_pval <- trend_test$sl[1]
    trend_tau <- trend_test$tau[1]
    
    
    
    dt_ssc <- SSC_annual_avg[site_no == site_no_sel]
    max_months <- max(dt_ssc$N, na.rm = T)
    first_reliable_yr <- min(dt_ssc[N >= (max_months * 0.6)]$year, na.rm = T)
    trend_test_ssc <- Kendall(dt_ssc[year >= first_reliable_yr]$year, dt_ssc[year >= first_reliable_yr]$SSC_mgL)
    trend_pval_ssc <- trend_test_ssc$sl[1]
    trend_tau_ssc <- trend_test_ssc$tau[1]
    
    dam_impact_sel <- dam_ssc_comb[year == 2020 & site_no == site_no_sel][
      ,.(Q_cms = mean(Q_cms, na.rm = T)),
      by = .(RiverName, ID, site_no, cap_mcm_cml)
    ][,':='(dam_storage_percent = cap_mcm_cml/(Q_cms * 3600 * 24 * 365.25/1e6) * 100)]
    deforest_impact_sel <- watershed_deforestation[site_no == site_no_sel & year == 2020][
      ,.(site_no, cumulative_deforestation_rel)
    ]
    
    trend_test_final <- cbind(dam_impact_sel[deforest_impact_sel, on = 'site_no'][
      ,':='(trend_pval = trend_pval,
            trend_tau = trend_tau,
            trend_pval_ssc = trend_pval_ssc,
            trend_tau_ssc = trend_tau_ssc,
            first_reliable_yr = first_reliable_yr,
            N_years = n_years)
    ], location_sel)
   
    # Only make plots if specified: significantly speeds up process if no plots 
    # (< 15 seconds vs. 22 minutes)
    if(make_ts_plots == 'Yes'){ 
    ts_plot <- ggplot(dt_sel, aes(x = year, y = tons_month_avg*12/1e6)) +
      geom_point(size = 3, fill = 'grey70', color = 'black', pch = 21) +
      geom_smooth(method = 'loess',formula = y~x, lty = 'dashed', color = 'black') +
      season_facet +
      scale_x_continuous(limits = c(1984, 2021)) +
      labs(
        x = '',
        y = 'Mt/yr',
        title = dt_sel$RiverName[1],
        subtitle = paste0(ifelse(trend_pval < 0.05, 'Significant, p = ','Not significant, p = '),
                          round(trend_pval, 3))
      )
    
    ts_plot_ssc <- ggplot(dt_ssc,
                          aes(x = year, y = SSC_mgL)) +
      geom_vline(xintercept = first_reliable_yr) +
      geom_point(size = 3, fill = 'light blue', color = 'black', pch = 21) +
      geom_smooth(data = dt_ssc[year >= first_reliable_yr], method = 'loess',formula = y~x, lty = 'dashed', color = 'black') +
      geom_text_repel(aes(label = N), size = 3, direction = 'y', point.padding = 1) +
      season_facet +
      scale_x_continuous(limits = c(1984, 2021)) +
      labs(
        x = '',
        y = 'SSC (mg/L)',
        subtitle = paste0(ifelse(trend_pval_ssc < 0.05, 'Significant, p = ','Not significant, p = '),
                          round(trend_pval_ssc, 3))
      )
    
    ts_plot_comb <- ggpubr::ggarrange(ts_plot, ts_plot_ssc, ncol = 1, labels = c('A','B'))
    ggsave(ts_plot_comb, filename = paste0(wd_exports, 'annual_Qss_ts_', gsub('/', '', dt_sel$RiverName[1]), '.png'),
           width = 5, height = 6.5) 
    }
    if(i == 1){
      trend_test_final_comb <- trend_test_final
    }else{
      trend_test_final_comb <- rbind(trend_test_final_comb, trend_test_final,
                                     use.names = T, fill = T)
    }
  }}

# Add a text column for trend sign
trend_test_final_comb <- trend_test_final_comb[
  ,':='(trend_sign = ifelse(trend_pval_ssc > 0.05, 'No change',
                            ifelse(trend_tau_ssc < 0, 'Decrease','Increase')),
        significance = ifelse(trend_pval_ssc < 0.05, 'Siginficant','Not significant'))
]

#### 5C. KENDALL NON-PARAMETRIC TREND SUMMARIES AND MAPS ####
# SSC trend summary
# Number and percent of rivers with increasing and decreasing trends
trend_test_summary <- trend_test_final_comb[!is.na(trend_sign)][,.(N = .N),
                                                                by = trend_sign][
                                                                  ,':='(`% of Total` = round(N/sum(N)*100,1))
                                                                ]
# Rivers in northern latitudes that are increasing
northern_increasing_rivers <- trend_test_final_comb[trend_sign == 'Increase' & Latitude > 20]

# Write final trend data to drive
fwrite(trend_test_final_comb, file = paste0(wd_exports, 'outlet_rivers_trend_test.csv'))

## Fig. 3A
# Map Kendall trend test results
ssc_trend_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = trend_test_final_comb[trend_pval_ssc > 0.05], stroke = 0.2,
             aes(x = Longitude, y = Latitude,
                 # fill = trend_tau_ssc # decided to leave non-significant trends white
                 size = Q_cms/1000,
                 shape = trend_sign,
                 ), fill = 'white') +
  geom_point(data = trend_test_final_comb[trend_pval_ssc < 0.05], stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = trend_tau_ssc,
                 size = Q_cms/1000,
                 shape = trend_sign)) +
  scale_fill_gradientn(limits = c(-0.75, 0.75), colors = rev(c(increase_color,'grey80',decrease_color)), oob = squish) +
  scale_size_continuous(range = c(2,10), breaks = c(1, 10, 100), guide = guide_legend(order = 1)) +
  scale_shape_manual(values = c('No change' = 21, 'Increase' = 24, 'Decrease' = 25),
                     guide = 'none') +
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
    fill = '**Kendall Tau statistic,<br>Sediment Flux**<br>**(1984-2020)**',
    size = '**Avg. Discharge**<br>**(1,000s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.005, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())


ggsave(ssc_trend_map, filename = paste0(wd_figures, 'ssc_trend_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(ssc_trend_map, filename = paste0(wd_figures, 'ssc_trend_map.png'),
       width = 11, height = 6)

#### 5D. KENDALL NON-PARAMETRIC TREND TEST, BY LATITUDE ####
# Trend summary by latitude
# Average trend, by latitude bin
trends_by_latitude <- trend_test_final_comb[,':='(Latitude_binned = Latitude - Latitude%%3)][
  ,.(avg_trend = mean(trend_tau_ssc, na.rm = T),
     N = .N),
  by = .(Latitude_binned)
][order(Latitude_binned)]


## Fig. 3B
# Plot average trend by latitude bins
trends_by_latitude_mean_kendall_tau_plot <- 
  ggplot(trends_by_latitude, aes(x = Latitude_binned, y = avg_trend)) +
  geom_bar(stat = 'identity', aes(fill = ifelse(avg_trend > 0, 'Increasing','Decreasing'))) +
  geom_line() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c(decrease_color,increase_color)) +
  # scale_fill_manual(values = c('#3399FF','#FABD5A','#DB592A')) +
  scale_x_continuous(limits = c(-61, 86), expand = expansion(mult = c(0, 0)),
                     breaks = seq(-60,80,20),
                     labels = lat_dd_lab,
                     sec.axis = dup_axis()) +
  scale_y_continuous(sec.axis = dup_axis()) +
  season_facet +
  theme(axis.text.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.x.top = element_blank(),
        axis.ticks.y.left = element_blank()) +
  labs(
    x = 'Latitude',
    y = 'Kendall Tau Statistic'
  ) +
  rotate()

ggsave(trends_by_latitude_mean_kendall_tau_plot, 
       filename = paste0(wd_figures, 'trends_by_latitude_mean_kendall_tau_plot.pdf'),
       width = 6, height = 6, useDingbats = F)
ggsave(trends_by_latitude_mean_kendall_tau_plot, 
       filename = paste0(wd_figures, 'trends_by_latitude_mean_kendall_tau_plot.png'),
       width = 6, height = 6)

## Fig. 3A, Fig. 3B
# Combine Kendall trends map (Fig. 3A) and by-latitude bar plots (Fig. 3B)
latitude_trends_map_and_panel <- ggpubr::ggarrange(
  ssc_trend_map + theme(axis.text.y.right = element_blank(),
                        axis.title.y.right = element_blank(),
                        axis.ticks.y.right = element_blank()),
  NULL, 
  trends_by_latitude_mean_kendall_tau_plot + 
    theme(plot.background = element_blank(),
          axis.title = element_text(size = 14),
          panel.background = element_blank()),
  align = 'hv', widths = c(1, -0.1, 0.3), nrow = 1,
  labels = c('A','','B'), label.x = c(0, 0, 0.1))

ggsave(latitude_trends_map_and_panel, filename = paste0(wd_figures, 'latitude_trends_map_and_panel.pdf'),
       width = 13, height = 6, useDingbats = F)
ggsave(latitude_trends_map_and_panel, filename = paste0(wd_figures, 'latitude_trends_map_and_panel.png'),
       width = 13, height = 6)


#### 6. CALCULATE TRENDS BY LATITUDE AND PLOT ####
#### 6A. CALCULATE ANOMALY TIMESERIES FOR SSC, PRECIP, TEMP ####
# Can take avg. based on all-time (N.B.: I'm not doing it this way: see next code block)
SSC_alltime_avg <- SSC_annual_avg[,.(SSC_mgL_avg_alltime = mean(SSC_mgL, na.rm = T),
                                     SSC_mgL_sd_alltime = sd(SSC_mgL, na.rm = T)),
                                  by = .(site_no, RiverName)]
# Can take avg. based on reference period: 1984-1990, inclusive (N.B.: I am doing it this way)
SSC_alltime_avg <- SSC_annual_avg[year %in% c(1985:1990)][
  ,.(SSC_mgL_avg_alltime = mean(SSC_mgL, na.rm = T),
     SSC_mgL_sd_alltime = sd(SSC_mgL, na.rm = T)),
  by = .(site_no, RiverName)]

SSC_annual_avg <- SSC_annual_avg[SSC_alltime_avg, on = c('site_no','RiverName')]

# Normalize annual SSC by all-time average for each river
SSC_annual_avg <- SSC_annual_avg[,':='(SSC_norm = SSC_mgL/SSC_mgL_avg_alltime,
                                       SSC_stdized = (SSC_mgL - SSC_mgL_avg_alltime)/
                                         SSC_mgL_sd_alltime)]

SSC_annual_avg_global <- SSC_annual_avg[,.(SSC_norm = mean(SSC_norm, na.rm = T)),
                                        by = year]

# Get annual precipitation anomaly
watershed_climate_site_avg <- watershed_climate[year %in% c(1985:1990)][,.(precip_m_yr_avg = mean(precip_m_yr, na.rm = T),
                                                                           precip_m_yr_sd = sd(precip_m_yr, na.rm = T),
                                                                           temp_c_avg = mean((mean_2m_air_temperature - 273.15), na.rm = T),
                                                                           temp_c_sd = sd((mean_2m_air_temperature - 273.15), na.rm = T)),
                                                                        by = .(site_no, RiverName)]

# Standardize anomaly in temperature and precipitation
watershed_climate_annual_avg <- watershed_climate[watershed_climate_site_avg,
                                                  on = c('site_no','RiverName')][
                                                    ,':='(precip_norm = precip_m_yr/precip_m_yr_avg,
                                                          precip_stdized = (precip_m_yr - precip_m_yr_avg)/
                                                            precip_m_yr_sd,
                                                          temp_norm = (mean_2m_air_temperature - 273.15)/temp_c_avg,
                                                          temp_stdized = ((mean_2m_air_temperature - 273.15) - temp_c_avg)/
                                                            temp_c_sd,
                                                          temp_anomaly = (mean_2m_air_temperature - 273.15) - temp_c_avg)]

# Combine annual SSC data with precipitation and temperature anomaly table
SSC_annual_avg <- SSC_annual_avg[watershed_climate_annual_avg[
  ,.(site_no, RiverName, year, precip_norm, precip_stdized, temp_norm, temp_stdized, temp_anomaly)
], on = c('site_no','RiverName','year')]

# Split data by north and south of 20° N Latitude
# Compute anomalies in SSC, temperature, and precipitation
SSC_annual_avg_north_south <- na.omit(SSC_annual_avg[
  ,':='(north_south = ifelse(Latitude > 20, 
                              'North of 20° N Latitude', 
                              'South of 20° N Latitude'))][
      ,.(SSC_norm = mean(SSC_norm, na.rm = T),
         SSC_stdized = mean(SSC_stdized, na.rm = T),
         precip_norm = mean(precip_norm, na.rm = T),
         precip_stdized = mean(precip_stdized, na.rm = T),
         temp_norm = mean(temp_norm, na.rm = T),
         temp_stdized = mean(temp_stdized, na.rm = T),
         temp_anomaly = mean(temp_anomaly, na.rm = T)),
      by = .(year, north_south)],
  cols = 'north_south')

## CALCULATION FOR PAPER ##
# Calculate recent anomaly for SSC
SSC_anomaly_post_2010 <- na.omit(SSC_annual_avg_north_south, cols = 'north_south',)[
  year > 2014,.(SSC_norm = mean(SSC_norm, na.rm = T),
                SSC_norm_se = sd(SSC_norm, na.rm = T)),
  by = .(north_south)]

# SSC chnge vs. temp. and precip. anomaly ts analysis
# Trend is not significant between precipitation and SSC for Southern rivers
Kendall(SSC_annual_avg_north_south[grepl('South', north_south)]$precip_norm, 
        SSC_annual_avg_north_south[grepl('South', north_south)]$SSC_norm)

#### 6B. PLOT SSC, PRECIP, TEMP ANOMALY TIMESERIES ####
## Fig. 3C
# Plot avg. ssc anomaly per year, by North vs. South
# Include Pinatubo excursion
ssc_standardized_anomaly_plot <- ggplot(SSC_annual_avg_north_south[year > 1984], aes(x = year)) + 
  # geom_smooth(span = 0.5, color = 'black', lty = dashed) +
  geom_ribbon(data = data.table(year = c(1991:1997), ymin = -Inf, ymax = Inf),
              aes(x = year, ymin = ymin, ymax = ymax), fill = 'grey80', color = NA) +
  # geom_smooth(aes(y = SSC_norm, group = north_south, color = north_south, fill = north_south),
  #             span = 0.3, lty = 'dashed', se = F, lwd = 0.5) +
  geom_smooth(data = SSC_annual_avg_north_south[year > 1984][
    ((year %in% c(1991:1997)) & north_south == 'North')
    ,':='(SSC_stdized = NA,
          SSC_norm = NA)],
    aes(y = SSC_norm, group = north_south, color = north_south),
    span = 0.3, lty = 'dashed', se = T, lwd = 0.5, fill = 'grey80') +
  geom_step(data = SSC_annual_avg_north_south[year > 1984],
            aes(y = SSC_norm, group = north_south, color = north_south)) +
  scale_color_manual(values = c('North of 20° N Latitude' = decrease_color, 'South of 20° N Latitude' = increase_color)) +
  scale_fill_manual(values = c('North of 20° N Latitude' = decrease_color, 'South of 20° N Latitude' = increase_color)) +
  season_facet +
  facet_wrap(.~north_south, ncol = 1) +
  labs(
    x = '',
    y = 'SSC*'
  )

## Fig. 3D
# Plot avg. precip anomaly per year, by North vs. South
# Include Pinatubo excursion
precipitation_anomaly_north_south_plot <- ggplot(SSC_annual_avg_north_south[year > 1984 & year < 2020],
                                                 aes(x = year)) + 
  # geom_smooth(span = 0.5, color = 'black', lty = dashed) +
  geom_ribbon(data = data.table(year = c(1991:1997), ymin = -Inf, ymax = Inf),
              aes(x = year, ymin = ymin, ymax = ymax), fill = 'grey80', color = NA) +
  # geom_smooth(aes(y = precip_norm, group = north_south, color = north_south, fill = north_south),
  #             span = 0.3, lty = 'dashed', se = F, lwd = 0.5) +
  geom_smooth(data = SSC_annual_avg_north_south[year > 1984 & year < 2020],
              aes(y = precip_norm, group = north_south, color = north_south, fill = north_south),
              span = 0.3, lty = 'dashed', se = F, lwd = 0.5) +
  geom_step(data = na.omit(SSC_annual_avg_north_south[year > 1984 & year < 2020], cols = 'north_south'),
            aes(y = precip_norm, group = north_south, color = north_south)) +
  scale_color_manual(values = c('North of 20° N Latitude' = decrease_color, 'South of 20° N Latitude' = increase_color)) +
  scale_fill_manual(values = c('North of 20° N Latitude' = decrease_color, 'South of 20° N Latitude' = increase_color)) +
  scale_x_continuous(labels = year_abbrev_lab) +
  season_facet +
  facet_wrap(.~north_south, ncol = 1) +
  labs(
    x = '',
    y = 'Precip*'
  )

## Fig. 3E
# Plot avg. temperature anomaly per year, by North vs. South
# Include Pinatubo excursion
temperature_anomaly_north_south_plot <- ggplot(SSC_annual_avg_north_south[year > 1984 & year < 2020],
                                               aes(x = year)) + 
  # geom_smooth(span = 0.5, color = 'black', lty = dashed) +
  geom_ribbon(data = data.table(year = c(1991:1997), ymin = -Inf, ymax = Inf),
              aes(x = year, ymin = ymin, ymax = ymax), fill = 'grey80', color = NA) +
  # geom_smooth(aes(y = temp_anomaly, group = north_south, color = north_south, fill = north_south),
  #             span = 0.3, lty = 'dashed', se = F, lwd = 0.5) +
  geom_smooth(data = SSC_annual_avg_north_south[year > 1984 & year < 2020],
              aes(y = temp_anomaly, group = north_south, color = north_south, fill = north_south),
              span = 0.3, lty = 'dashed', se = F, lwd = 0.5) +
  geom_step(data = SSC_annual_avg_north_south[year > 1984 & year < 2020],
            aes(y = temp_anomaly, group = north_south, color = north_south)) +
  scale_color_manual(values = c('North of 20° N Latitude' = decrease_color, 'South of 20° N Latitude' = increase_color)) +
  scale_fill_manual(values = c('North of 20° N Latitude' = decrease_color, 'South of 20° N Latitude' = increase_color)) +
  scale_x_continuous(labels = year_abbrev_lab) +
  season_facet +
  facet_wrap(.~north_south, ncol = 1) +
  labs(
    x = '',
    y = 'Temperature anomaly (°C)'
  )

# Blank plot with year labels
year_label_plot <- ggplot() + 
  geom_blank() +
  season_facet +
  theme_transparent() +
  theme(axis.title.x = element_text(size = 14)) +
  labs(x = 'Year')

# Combine pieces of Fig. 3
# Lower panel (Fig. 3C-E)
global_ssc_precip_t_trends_by_latitude <- ggpubr::ggarrange(ggpubr::ggarrange(ssc_standardized_anomaly_plot +
                                                                theme(strip.text = element_text(size = 14),
                                                                      axis.title.y = element_text(size = 14),
                                                                      axis.text.x = element_text(size = 14)),
                                                              precipitation_anomaly_north_south_plot + 
                                                                theme(strip.text = element_blank(),
                                                                      axis.title.y = element_text(size = 14),
                                                                      axis.text.x = element_text(size = 14)), 
                                                              temperature_anomaly_north_south_plot + 
                                                                theme(strip.text = element_blank(),
                                                                      axis.title.y = element_text(size = 14),
                                                                      axis.text.x = element_text(size = 14)),
                                                              nrow = 1, align = 'hv', widths = c(1, 0.5, 0.5), 
                                                              labels = c('C','D','E')),
                                                    NULL, year_label_plot, ncol = 1, heights = c(1, -0.04, 0.05))

## Fig. 3
# Combine upper (map and latitude trends) 
# and lower (SSC and climate anomalies)
# for final plot

latitude_trends_figure <- ggpubr::ggarrange(
  latitude_trends_map_and_panel, 
  NULL,
  global_ssc_precip_t_trends_by_latitude, 
  ncol = 1,
  heights = c(1, -0.01, 0.6))

ggsave(latitude_trends_figure, filename = paste0(wd_figures, 'Fig3_latitude_trends_figure.pdf'),
       width = 13, height = 10, useDingbats = F)
ggsave(latitude_trends_figure, filename = paste0(wd_figures, 'Fig3_latitude_trends_figure.png'),
       width = 13, height = 10)


#### 7. LAND USE CHANGE AND DEFORESTATION ANALYSIS AND MAP ####
# Calculate increase in deforestation since 2000
pre_post_deforestation <- watershed_deforestation[year == 2020][
  ssc_rivers_stats_full[,.(RiverName, Q_cms_avg, Latitude, Longitude)], on = 'RiverName']

# Watershed deforestation for each river in 2020
watershed_deforestation_2020 <- watershed_deforestation[year == 2020][
  ssc_rivers_stats_full[Continent_Region != ''], on = c('RiverName', 'site_no','ID')
][
  ,':='(Continent_Region = ifelse(Continent_Region %in% c('Europe','Eurasia'),
                                  'Europe/Eurasia', 
                                  Continent_Region))]

# Watershed deforestation by continent
watershed_deforestation_continent <- na.omit(watershed_deforestation_2020[
  ,.(deforestation_km2 = sum(cumulative_deforestation_km2, na.rm = T)),
  by = .(Continent_Region)])

## Fig. S11
# Map deforestation change for each watershed
watershed_deforestation_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = pre_post_deforestation, 
             stroke = 0.2, pch = 21,
             aes(x = Longitude, y = Latitude, 
                 fill = cumulative_deforestation_rel,
                 size = Q_cms_avg/100)) +
  scale_fill_gradientn(limits = c(0, 30), oob = squish,
                       colors = c('#222601','#565902','#BF9004','#D9B504'),
                       guide = guide_colorbar(order = 1)) +
  scale_size_continuous(range = c(2,10), breaks = c(100, 500, 1000)) +
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
    fill = '**% Deforested**<br>2000-2020',
    size = '**Avg. Discharge**<br>**(100s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(watershed_deforestation_map, filename = paste0(wd_figures, 'FigS11_watershed_deforestation_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(watershed_deforestation_map, filename = paste0(wd_figures, 'FigS11_watershed_deforestation_map.png'),
       width = 11, height = 6)

## Land use change from my (Evan Dethier) survey in 2021
land_use_survey_data <- data.table(`Land Use` = c(
  'Land Use Change',
  'Deforestation',
  'Alluvial Gold Mining',
  'Sand Mining',
  'Palm oil',
  'Agriculture (not palm oil)',
  'No Land Use Change',
  'Discharge'),
  `N Rivers` = c(58, 32, 24, 13, 12, 4, 12, 2),
  color = c(rep('Land Use Change', 6), 'None','None'),
  `Land Use Cause` = factor(c('Overview',rep('Observed Land Use Change', 5), 'Overview','None'), 
                            levels = c('Overview','Observed Land Use Change', 'None')))[
                              ,':='(`Land Use` = fct_reorder(`Land Use`, `N Rivers`))
                            ]


land_use_change_bar_chart <- ggplot(land_use_survey_data, aes(x = `Land Use`, y = `N Rivers`, fill = color)) +
  geom_col(alpha = 0.8, width = 0.85) +
  facet_grid(rows = vars(`Land Use Cause`), scales = 'free_y', switch = "y", space = "free_y") +
  coord_flip() +
  scale_fill_manual(values = rev(c('#222601','#D9B504'))) +
  season_facet +
  theme_minimal() +
  labs(
    title = "Land Use Change in Watersheds with SSC increase",
    subtitle = "",
    caption = "",
    y = "Number of Rivers"
  ) +
  # theme_minimal(base_family = "Roboto Condensed") +
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
    plot.title = element_text(face = "bold", hjust = 0),
    strip.text.y = element_text(angle = 270, face = "bold"),
    strip.placement = "outside",
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
  )

ggsave(land_use_change_bar_chart, filename = paste0(wd_figures, 'land_use_change_bar_chart.pdf'),
       width = 5, height = 5, useDingbats = F)

#### 8. PRE- AND POST- SEDIMENT FLUX ANALYSIS ####
#### 8A. SEDIMENT FLUX CHANGE CALCULATED BY CONTINENT ####
# Pre- and post-dam sediment flux, by continent
Qss_pre_post_dam_continent <- dcast.data.table(Continent_Region ~ year, value.var = 'tons_yr',
                                               data = Qss_global_annual_ref[year %in% c(1970, 2020)])

# Pre- and post-dam change in sediment flux UNCERTAINTY, by continent
Qss_post_dam_continental_rel_error <- Qss_global_annual_ref[year %in% c(2015:2020)][
  ,.(tons_yr_se_rel = sqrt(sqrt(sum(tons_yr_se_rel^2))/6/sqrt(6) + 0.5^2)),
  by = .(Continent_Region)]

# Half-decadal change in suspended flux, by continent
# Includes UNCERTAINTY
Qss_global_annual_ref_decadal <- Qss_global_annual_ref[
  ,':='(decade = ifelse(year == 1984, 1985, year - year%%5))
][
  ,.(tons_yr = mean(tons_yr, na.rm = T),
     tons_yr_se = sd(tons_yr, na.rm = T)/sqrt(.N)),
  by = .(decade, Continent_Region)
][
  ,':='(tons_yr_se = ifelse(is.na(tons_yr_se), 0, tons_yr_se))
]

# Change in suspended sediment flux, by continent
# (since deforestation quantified by Hansen et al., 2012)
Qss_pre_post_deforestation <- dcast.data.table(Continent_Region ~ decade, value.var = c('tons_yr', 'tons_yr_se'),
                                               data = Qss_global_annual_ref_decadal[decade %in% c(2015, 1995)])
# (since beginning of Landsat record)
Qss_pre_post_landsat <- dcast.data.table(Continent_Region ~ decade, value.var = c('tons_yr', 'tons_yr_se'),
                                         data = Qss_global_annual_ref_decadal[decade %in% c(2015, 1990)])


# A. Combined change in sediment flux, dam building, and watershed deforestation
# (Since pre-dam)
dam_deforest_impacts <- dam_timeseries_continent[year == 2020][
  Qss_pre_post_dam_continent, on = 'Continent_Region'
][watershed_deforestation_continent, on = 'Continent_Region']

# Add Reservoir capacity index (RCI)
# (Fraction of total annual discharge storable in reservoirs)
dam_deforest_impacts <- dam_deforest_impacts[,':='(RCI = storage_mcm_cml/Q_km3_yr * 0.001 * 100,
                                                   delSSL_2020_1970 = (`2020` - `1970`)/`1970` * 100)]



# B. Combined change in sediment flux, dam building, and watershed deforestation
# (Since 2000)
dam_deforest_impacts_2000s <- dam_timeseries_continent[year == 2020][
  Qss_pre_post_deforestation, on = 'Continent_Region'
][watershed_deforestation_continent, on = 'Continent_Region']

# (Since start of Landsat 5 record, 1990)
dam_deforest_impacts_landsat <- dam_building_landsat_era[
  Qss_pre_post_landsat, on = 'Continent_Region'
][watershed_deforestation_continent, on = 'Continent_Region'][
  continent_area_Q, on = 'Continent_Region'
]

# Add UNCERTAINTY columns to Landsat record beginning and end comparison
dam_deforest_impacts_landsat <- dam_deforest_impacts_landsat[
  ,':='(Qss_change_percent_2020_1984 = (tons_yr_2015 - tons_yr_1990)/tons_yr_1990 * 100,
        Qss_change_percent_2020_1984_se = sqrt((tons_yr_se_2015/tons_yr_2015)^2 + 
                                                 (tons_yr_se_1990/tons_yr_1990)^2 + 
                                                 (tons_yr_se_1990/tons_yr_1990)^2) * 100)
]

## CALCULATION FOR PAPER ##
contintental_dam_impacts_lm <- lm(delSSL_2020_1970 ~ RCI, data = dam_deforest_impacts)
summary(contintental_dam_impacts_lm)

#### 8B. CHANGE CALCULATED GLOBALLY ####
## CALCULATION FOR PAPER ##
# Figures and tables are not used except here
# Pre- and post-dam change in sediment flux, global
Qss_pre_post_dam_global <- (sum(Qss_pre_post_dam_continent$`2020`, na.rm = T) - 
                              sum(Qss_pre_post_dam_continent$`1970`, na.rm = T))/
  sum(Qss_pre_post_dam_continent$`1970`, na.rm = T) * 100

# # Pre- and post-dam change in sediment flux UNCERTAINTY, global
Qss_post_dam_global_rel_error <- sqrt(sum(Qss_global_annual_ref[year %in% c(2015:2020)]$tons_yr_se_rel^2))/6/sqrt(6)
Qss_pre_post_dam_global_rel_error <- sqrt(Qss_post_dam_global_rel_error^2 + 0.50^2)*49


## CALCULATION FOR PAPER ##
# Different way of calculating the same as above
# (Gets the same answer)
Qss_global_pre_post_dam <- Qss_global_annual_ref[decade %in% c(1970, 2015)][
  ,':='(decadal_sum = sum(tons_yr, na.rm = T)/5),
  by = .(decade)
][
  ,':='(continent_percent = ifelse(decade == 1970, tons_yr/(decadal_sum*5),
                                   tons_yr/decadal_sum))
]

Qss_global_pre_post_dam[,.(Gt_yr = sum(tons_yr, na.rm = T)/1e9),
                        by = year]

## CALCULATION FOR PAPER ##
Qss_global_pre_post_dam_avg <- Qss_global_pre_post_dam[
  ,.(tons_yr = mean(tons_yr, na.rm = T),
     tons_yr_sd = sd(tons_yr, na.rm = T),
     continent_percent = mean(continent_percent, na.rm = T)*100,
     continent_percent_se = sd(continent_percent)*100),
  by = .(Continent_Region, decade)
]

Qss_global_pre_post_dam_total <- Qss_global_pre_post_dam_avg[
  ,.(Mt_yr = sum(tons_yr, na.rm = T)/1e6,
     Mt_yr_sd = sqrt(sum(tons_yr_sd^2, na.rm = T))/1e6),
  by = .(decade)
]

Qss_global_post_dam_percent_reduction <- round(Qss_global_pre_post_dam_total[decade == 2015]$Mt_yr/
                                                 Qss_global_pre_post_dam_total[decade == 1970]$Mt_yr * 100)


#### 8C. PLOTS OF CHANGES IN SED. FLUX, RES. CAPACITY, DEFOREST. ####
# Fig. 4A
# Timeseries plot of cumulative reservoir storage (RCI)
# 1920 to present
# (By continent)
continent_dam_annual_ts_plot <- ggplot(
  dam_timeseries_continent, 
  aes(x = year)) +
  geom_step(aes(y = storage_mcm_cml/Q_km3_yr * 0.001 * 100, 
                group = Continent_Region), color = 'black') +
  # geom_point(data = dam_timeseries_continent[year == 2020],
  #            aes(x = year, y = storage_mcm_cml/Q_km3_yr * 0.001)) +  
  geom_text_repel(nudge_y = 0.015, data = dam_timeseries_continent[year == 2020],
                  aes(x = year, y = storage_mcm_cml/Q_km3_yr * 0.001 * 100, 
                      label = Continent_Region),
                  box.padding = unit(0.5, "lines")) +
  # geom_line(aes(y = tons_yr/1e6, group = Continent_Region)) +
  # geom_point(pch = 21, color = 'black', aes(fill = Continent_Region),
  #            size = 3) +
  # geom_smooth(data = Qss_global_annual_ref[year > 1970 & Continent_Region == continent_sel],
  #             aes(y = tons_yr/1e6, group = Continent_Region), method = 'loess', formula = y~x,
  #             color = 'black') +
  season_facet +
  # theme(legend.position = 'top') +
  scale_x_continuous(breaks = c(1920, 1940,1960,1980, 2000,2020),
                     expand = expansion(mult = c(0,0.01))) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
  #                    limits = c(0,4000)) + # Keep all y-axis scales equal
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))
  ) +
  guides(color = 'none') +
  # scale_y_log10() +
  # facet_wrap(.~Continent_Region) +
  theme(
    strip.text = element_text(size = 13) 
  ) +
  labs(
    x = 'Year',
    y = 'Reservoir capacity as % of annual discharge',
    fill = 'Continent'
  )

# Fig. 4B
# Dam capacity vs. change in sediment flux
# Pre-dam vs. present
# (by continent)
dam_capacity_vs_qss_reduction_continental <- ggplot(dam_deforest_impacts,
                                                    # [Qss_post_dam_continental_rel_error, on = 'Continent_Region'], 
                                                    aes(x = RCI, y = delSSL_2020_1970)) +
  # geom_linerange(aes(ymin = delSSL_2020_1970 - delSSL_2020_1970*tons_yr_se_rel,
  #                    ymax = delSSL_2020_1970 + delSSL_2020_1970*tons_yr_se_rel)) +
  geom_point(size = 4) +
  geom_smooth(lty = 'dashed', color = 'black', method = 'lm', formula = y ~ x, se = F) +
  geom_text_repel(aes(label = Continent_Region),
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines")) +
  season_facet +
  labs(
    x = 'Reservoir capacity as % of annual discharge',
    y = '% Change in sediment flux, 2020 vs. pre-dam estimate'
  )

# NEW Dam capacity vs. change in sediment flux
# 1984 vs. present
# (by continent)
dam_capacity_recent_vs_qss_reduction_continental <- ggplot(dam_deforest_impacts_landsat, 
                                                           aes(x = (storage_mcm_cml_2020 - storage_mcm_cml_1984)/Q_km3_yr * 0.001  * 100, 
                                                               y = Qss_change_percent_2020_1984)) +
  geom_point(size = 4) +
  geom_smooth(lty = 'dashed', color = 'black', method = 'lm', formula = y ~ x, se = F) +
  geom_text_repel(aes(label = Continent_Region),
                  box.padding = unit(0.5, "lines")) +
  season_facet +
  labs(
    x = 'New reservoir capacity as % of annual discharge',
    y = '% Change in sediment flux, 2020 vs. 1984'
  )

ggsave(dam_capacity_vs_qss_reduction_continental, filename = paste0(wd_figures, 'dam_capacity_vs_qss_reduction_continental.pdf'),
       width = 5, height = 5, useDingbats = F)



# Deforestation vs. changes in sediment flux
# (Not used in paper)
deforestation_vs_qss_reduction_continental <- ggplot(
  dam_deforest_impacts_2000s, aes(
      x = deforestation_km2/drainage_area_km2, 
      y = (tons_yr_2015 - tons_yr_1995)/tons_yr_1995 * 100)) +
  geom_point(size = 4) +
  geom_smooth(lty = 'dashed', color = 'black', method = 'lm', formula = y ~ x, se = F) +
  geom_text_repel(aes(label = Continent_Region),
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines")) +
  season_facet +
  labs(
    x = 'Deforestation as % of watershed area',
    y = '% Change in sediment flux, 2020 vs. 2000'
  )



#### 8D. MAP OF RESERVOIR CAPACITY INDEX (RCI) ####
# Dam reservoir storage in 2020
pre_post_dam <- dam_timeseries[year == 2020][
  ssc_rivers_stats_full[,.(RiverName, Q_cms_avg, Latitude, Longitude, drainage_area_km2)], on = 'RiverName']

# Dam reservoir area in 1984
dam_area_1984 <- dam_timeseries[year == 1984][,
                                              ':='(cap_mcm_cml_1984 = cap_mcm_cml,
                                                   area_km2_cml_1984 = area_km2_cml)][
                                                     ,.(cap_mcm_cml_1984, area_km2_cml_1984, site_no, RiverName)
                                                   ]

# Dam reservoir 1984 and present (changes in Landsat record)
pre_post_dam <- pre_post_dam[dam_area_1984, on = c('RiverName','site_no')]

# Fig. 4C
# Map reservoir capacity (RCI) in 2020
watershed_dam_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = pre_post_dam, 
             stroke = 0.2, pch = 21,
             aes(x = Longitude, y = Latitude, 
                 fill = cap_mcm_cml/(Q_cms_avg * 60 * 60 * 24 * 365.25 / 1e6) * 100,
                 # size = cut(Q_cms_avg/1000, breaks = c(0, 0.1, 1, 10, 100, Inf)))) +
                 size = Q_cms_avg/1000)) +
  scale_fill_gradientn(limits = c(0, 50), oob = squish,
                       colors = c('#7475B0','#ACCAE3','#B08874','#E6AC76')) +
  scale_size_continuous(breaks = c(1, 10, 100), range = c(2,10), 
                        guide = guide_legend(order = 1, override.aes = list(fill = 'black'))) +
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
    fill = '**Reservoir capacity**<br>as % of annual discharge',
    size = '**Avg. Discharge**<br>(1,000s m<sup>3</sup>/s)'
  ) +
  theme(legend.position = c(0.005, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(watershed_dam_map, filename = paste0(wd_figures, 'watershed_dam_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(watershed_dam_map, filename = paste0(wd_figures, 'watershed_dam_map.png'),
       width = 11, height = 6)


dam_panel_plot <- ggpubr::ggarrange(continent_dam_annual_ts_plot, NULL, 
                            dam_capacity_vs_qss_reduction_continental +
                              labs(y = '% Change sediment flux, 2020 vs. pre-dam'),
                            widths = c(1,0.05,1), nrow = 1, align = 'v', labels = c('A','','B'))

ggsave(dam_panel_plot, filename = paste0(wd_figures, 'dam_panel_plot.pdf'),
       width = 10, height = 5, useDingbats = F)
ggsave(dam_panel_plot, filename = paste0(wd_figures, 'dam_panel_plot.png'),
       width = 10, height = 5)

## Fig. 4
# Combine reservoir capacity (RCI) plots to make Fig. 4
# Fig. 4A: Reservoir capacity index over time
# Fig. 4B: Relation between RCI and % change in sediment flux
# Fig. 4C: Map of RCI for each river
dam_panel_plot_w_map <- ggpubr::ggarrange(dam_panel_plot, NULL, watershed_dam_map,
                                  nrow = 3, align = 'v', labels = c('','', 'C'),
                                  heights = c(0.8,-0.01, 1))

ggsave(dam_panel_plot_w_map, filename = paste0(wd_figures, 'Fig4_dam_panel_w_map_plot.pdf'),
       width = 11, height = 10, useDingbats = F)
ggsave(dam_panel_plot_w_map, filename = paste0(wd_figures, 'Fig4_dam_panel_w_map_plot.png'),
       width = 11, height = 10)


#### 9. PLOT ANNUAL, DECADAL CHANGES IN SEDIMENT FLUX ####
#### 9A. PLOT CHANGE IN SSC: 2020 vs. 1990 ####
# Add discharge to pre and post SSC table
pre_post_ssc <- Qss_river_decadal_change[ssc_rivers_stats_full[,.(RiverName, Q_cms_avg)], on = 'RiverName']

## Fig. 2B
# Percent change in discharge-weighted ssc
# 2015-2021 vs. 1984-1991
ssc_change_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = pre_post_ssc, stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = tons_month_2015_1990,
                 size = Q_cms_avg/1000,
                 shape = Qss_sign)) +
  scale_fill_gradientn(values = c(0, 0.4, 0.5, 0.6, 1), limits = c(-75, 75), colors = rev(c(increase_color,'grey80','grey80','grey80',decrease_color)), oob = squish) +
  # scale_size_continuous(range = c(2,10), breaks = c(100, 500, 1000)) +
  scale_size_continuous(range = c(2,10), breaks = c(1, 10, 100), guide = guide_legend(order = 1)) +
  scale_shape_manual(values = c('No change' = 21, 'Increase' = 24, 'Decrease' = 25),
                     guide = 'none') +
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
    fill = '**% Change,\nSediment Flux**<br>**(2015-2021 vs. 1984-1991)**',
    size = '**Avg. Discharge**<br>**(1,000s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())


# Combine change with SSC average
# (Not used in paper)
ssc_avg_change_map_comb <- ggpubr::ggarrange(avg_ssc_map, ssc_change_map,
                                     ncol = 1, align = 'hv',
                                     labels = c('A','B'))

ggsave(ssc_avg_change_map_comb, filename = paste0(wd_figures, 'ssc_avg_change_map_comb.pdf'),
       width = 12.25, height = 12, useDingbats = F)



#### 9B. BAR PLOT BY-CONTINENT ANNUAL FLUX ####
continent_list <- unique(Qss_global_annual_ref$Continent_Region)
continent_qss_annual_ts_plotlist <- vector('list', length(continent_list))
names(continent_qss_annual_ts_plotlist) <- continent_list
for(i in 1:length(continent_list)){
  continent_sel <- continent_list[i]
  continent_qss_annual_ts_plot <- ggplot(
    Qss_global_annual_ref[Continent_Region == continent_sel & year > 1983], 
    aes(x = year)) +
    # geom_errorbar(aes(ymin = tons_yr/1e6 - tons_yr*tons_yr_se_rel/1e6,
    #                 ymax = tons_yr/1e6 + tons_yr*tons_yr_se_rel/1e6,
    #                 fill = Continent_Region, group = Continent_Region),
    #             alpha = 0.25) +
    geom_point(data = Qss_global_annual_ref[Continent_Region == continent_sel & year < 1983],
               aes(x = 1970, y = tons_yr/1e6), size = 7) +
    geom_line(data = Qss_global_annual_ref[Continent_Region == continent_sel & year < 1985],
              aes(x = year, y = tons_yr/1e6), lty = 'dashed') +
    geom_bar(stat = 'identity', aes(y = tons_yr/1e6, group = Continent_Region)) +
    # geom_line(aes(y = tons_yr/1e6, group = Continent_Region)) +
    # geom_point(pch = 21, color = 'black', aes(fill = Continent_Region),
    #            size = 3) +
    geom_linerange(aes(x = year, ymin = tons_yr/1e6 - tons_yr/1e6*tons_yr_se_rel,
                       ymax = tons_yr/1e6 + tons_yr/1e6*tons_yr_se_rel)) +
    geom_vline(xintercept = 1975) +
    geom_smooth(data = Qss_global_annual_ref[year > 1970 & Continent_Region == continent_sel],
                aes(y = tons_yr/1e6, group = Continent_Region), method = 'loess', formula = y~x,
                color = 'black') +
    season_facet +
    # theme(legend.position = 'top') +
    scale_x_continuous(breaks = c(1970, 1980,1990,2000,2010,2020),
                       labels = c('Pre-dam',1980,1990,2000,2010,2020)
    ) +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
    #                    limits = c(0,4000)) + # Keep all y-axis scales equal
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))
    ) +
    guides(color = 'none') +
    # scale_y_log10() +
    facet_wrap(.~Continent_Region) +
    theme(
      strip.text = element_text(size = 13) 
    ) +
    labs(
      x = 'Year',
      y = 'Mt/yr',
      fill = 'Continent'
    )
  
  continent_qss_annual_ts_plotlist[[i]] <- continent_qss_annual_ts_plot
}

## Fig. 2A
# Northern landmass bar plots of annual continental sediment flux
northern_hemisphere_qss_panels <- ggpubr::ggarrange(plotlist = list(continent_qss_annual_ts_plotlist$`North America`,
                                                            continent_qss_annual_ts_plotlist$`Europe/Eurasia`,
                                                            continent_qss_annual_ts_plotlist$Asia),
                                            nrow = 1, align = 'hv')

## Fig. 2C
# Southern landmass bar plots of annual continental sediment flux
southern_hemisphere_qss_panels <- ggpubr::ggarrange(plotlist = list(continent_qss_annual_ts_plotlist$`South America`,
                                                            continent_qss_annual_ts_plotlist$Africa,
                                                            continent_qss_annual_ts_plotlist$Oceania),
                                            nrow = 1, align = 'hv')

## Fig. 2
# Combine bar plots with SSC change map
continental_change <- ggpubr::ggarrange(
  northern_hemisphere_qss_panels,
  ssc_change_map,
  southern_hemisphere_qss_panels,
  ncol = 1,
  heights = c(0.7, 1.3, 0.7),
  labels = c('A','B','C')
)
ggsave(continental_change, filename = paste0(wd_figures, 'Fig2_continental_change_panel_plot_2022.pdf'),
       width = 12, height = 13, useDingbats = F)
ggsave(continental_change, filename = paste0(wd_figures, 'Fig2_continental_change_panel_plot_2022.png'),
       width = 12, height = 13)

#### 10. COLOR AND CLUSTER METADATA VISUALIZATION ####
# Need discharge data at this point
# Use color for each month at each river
# Add discharge to this table, to scale bars
site_monthly_band_median_color <- site_monthly_band_median_color[
  avg_monthly_Q_Qss, on = c('site_no','month')
][,':='(Continent_Region = ifelse(Continent_Region %in% c('Europe','Eurasia'),
                                  'Europe/Eurasia', Continent_Region))]

# Labels for month plots
month_labels <- c('J','F','M','A','M','J','J','A','S','O','N','D')

#### 10A. WHEEL PLOT EXAMPLES FOR EACH CONTINENT ####
# Set up plots
# Hemispheres (north and south)
hemispheres <- c('Northern','Southern')
# Continent names, in plot order
continent_order <- c('North America','Europe/Eurasia', 'Asia',
                     'South America', 'Africa','Oceania')
# Drainage area and discharge for each continent
# (Adding a hemisphere column to table as well)
continent_area_Q <- continent_area_Q[,':='(hemisphere = ifelse(
  Continent_Region %in% c('South America', 'Africa','Oceania'), 'Southern','Northern'))]


# Select rivers for example wheel plots
# Three from each landmass
color_wheel_plot_rivers <- c(
  'Amazon','Magdalena','Paraná', # S. America
  'Mississippi','MacKenzie','Potomac', # N. America
  'Amur','Huanghe','Mekong', # Asia
  'Ikopa', 'Orange','Senegal', # Africa
  'Danube','Po','Rhone', #Europe/Eurasia
  'Murray-Darling','Burdekin','Purari') # Oceania

# Set up blank wheel plots and name them by landmass
wheel_plots_continent <- vector('list', nrow(continent_area_Q))
names(wheel_plots_continent) <- continent_order

# Loop through each landmass: make a wheel plot triple for each
# triples is best
for(i in 1:6){
  # Select a landmass from list of 6
  continent_sel <- continent_area_Q[i]
  continent_sel_region <- continent_sel$Continent_Region
  # Select three rivers from that landmass
  site_monthly_band_median_color_sel <- site_monthly_band_median_color[
    RiverName %in% color_wheel_plot_rivers & 
      Continent_Region %in% continent_sel_region]
  
  # Make a wheel plot for three selected rivers from that landmass
  wheel_plots_indiv_continental <- ggplot(site_monthly_band_median_color_sel) +
    geom_bar(width = 1, color = 'grey50', stat = 'identity', lwd = 0.25,
             aes(x = month, y = Q_cms_norm, fill = rgb(B3/2100, B2/2100, B1/2100))) +
    season_facet +
    facet_wrap(dir = 'v', 
               .~reorder(RiverName, Longitude_mil),
               # Continent_Region~RiverName, 
               nrow = 1, drop = T
    ) +
    scale_fill_identity() +
    scale_y_sqrt(limits = c(0,1)) +
    scale_x_continuous(breaks = c(1:12), labels = month_labels, expand = expansion(mult = 0)) +
    # theme(
    #   panel.grid.major = element_line(color = 'black')
    # ) +
    coord_polar() +
    theme_void() +
    theme(strip.text.x=element_text(margin=margin(b=3), size = 10),
          title = element_text(size = 8, face = 'bold'),
    ) +
    labs(title = continent_sel_region)
  
  which_wheel_plot <- which(names(wheel_plots_continent) == continent_sel_region)
  wheel_plots_continent[[which_wheel_plot]] <- wheel_plots_indiv_continental
}

## Fig. 1
# Combine example color wheels (Fig. 1A, Fig. 1C)
# with Avg. SSC Map (Fig. 1B)
ssc_map_with_wheels <- ggpubr::ggarrange(plotlist = list(
  ggpubr::ggarrange(plotlist = list(
    wheel_plots_continent$`North America`,
    wheel_plots_continent$`Europe/Eurasia`,
    wheel_plots_continent$Asia),
    nrow = 1),
  avg_ssc_map,
  NULL,
  ggpubr::ggarrange(plotlist = list(
    wheel_plots_continent$`South America`,
    wheel_plots_continent$Africa,
    wheel_plots_continent$Oceania),
    nrow = 1)
),
ncol = 1, 
heights = c(1,4,-0.02,1),
labels = c('A','B','','C')
)

ggsave(ssc_map_with_wheels, filename = paste0(wd_figures, 'Fig1_ssc_map_with_wheels.pdf'),
       width = 11, height = 8, useDingbats = F)
ggsave(ssc_map_with_wheels, filename = paste0(wd_figures, 'Fig1_ssc_map_with_wheels.png'),
       width = 11, height = 8)

#### 10B. BY-RIVER CLUSTER PARALLEL PLOT VISUALIZATION ####
# Fig. S2A
# Plot cluster for each outlet river station
outlet_river_cluster_plot <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = river_cluster_assignment_single,
             aes(x = Longitude, y = Latitude, 
                 fill = as.factor(cluster)), size = 3, pch = 21, stroke = 0.2) + 
  scale_fill_manual(values = c('#1f78b4', '#33a02c','#e31a1c','#ff7f00','#6a3d9a','#a6cee3')) +
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
    fill = '**River**<br>**Cluster**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        # legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

# Convert by-site band data (each band average value is a column)
# from wide to long
# Now each river has a site for each band value
# Used to visualize cluster data in parallel plots
cluster_visualization_example <- melt(
  merge(
    site_monthly_band_median_color[,.(
      B1 = mean(B1, na.rm = T),
      B2 = mean(B2, na.rm = T),
      B3 = mean(B3, na.rm = T),
      B4 = mean(B4, na.rm = T),
      B5 = mean(B5, na.rm = T),
      B7 = mean(B7, na.rm = T)),
      by = .(site_no)],
    river_cluster_assignment_single[, .(site_no, cluster)],
    by = 'site_no'),
  id.vars = c('site_no', 'cluster'))

## Fig S2B
# Plot cluster parallel plots
# Each river is a single line, each cluster is a panel
outlet_river_cluster_parallel_plot <- ggplot(
  cluster_visualization_example,
  aes(x = variable, y = value/10000, 
      # alpha = ifelse(RiverName == 'Congo',1,0.2),
      # color = ifelse(RiverName %in% c('Congo','Nottaway'),'Of interest','Other'))) +
      color = factor(cluster))) +
  geom_line(aes(group = paste0(site_no,cluster)),
            alpha = 0.4) +
  # scale_alpha_identity() +
  scale_color_manual(values = c('#1f78b4', '#33a02c','#e31a1c','#ff7f00','#6a3d9a','#a6cee3')) +
  season_facet +
  facet_wrap(.~paste0('Cluster ', cluster)) +
  labs(
    x = 'Landsat band',
    y = 'Surface reflectance'
  )

## Fig. S2C
# Plot of river color for each cluster
# Broken out by SSC range in color palette
max_reflectance <- max(cluster_ssc_cat_colors[,.(blue, red, green)])*1.5
# Set up plot for either true color or false color
raster_color_types <- c(geom_raster(aes(fill = rgb(red/max_reflectance,green/max_reflectance,blue/max_reflectance))), # true color
                        geom_raster(aes(fill = rgb(B4/4000,B3/4000,B2/4000))) # false color)
)

cluster_ssc_category_color_plot <- 
  ggplot(cluster_ssc_cat_colors, aes(x = cluster, y = SSC_cat)) +
  # ggplot(ssc_category_color, aes(x = reorder(paste0(cluster_sel, ' ', site_no), cluster_sel), y = ssc_category)) + # for by site
  raster_color_types[1] +
  scale_fill_identity() +
  season_facet + 
  # scale_x_continuous(expand_scale(add = c(0,0))) + 
  # scale_y_discrete(expand_scale(mult = c(0,0))) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(
    y = 'SSC range (mg/L)',
    x = 'River cluster'
  )




## Fig. S2
# Combine plots to make Fig. S2
cluster_metaanalysis_plots <- outlet_river_cluster_plot / 
  (outlet_river_cluster_parallel_plot + cluster_ssc_category_color_plot +
     plot_layout(widths = c(1.5, 1))) +
  plot_annotation(tag_levels = 'A')

ggsave(cluster_metaanalysis_plots, filename = paste0(wd_figures, 'FigS2_cluster_metaanalysis_plots.pdf'),
       width = 10, height = 11, useDingbats = F)
ggsave(cluster_metaanalysis_plots, filename = paste0(wd_figures, 'FigS2_cluster_metaanalysis_plots.png'),
       width = 10, height = 11)

ggsave(outlet_river_cluster_plot, filename = paste0(wd_figures, 'outlet_river_cluster_plot.pdf'),
       width = 11, height = 6)
ggsave(outlet_river_cluster_plot, filename = paste0(wd_figures, 'outlet_river_cluster_plot.png'),
       width = 11, height = 6)

############## III. CONTEXTUALIZE CHANGES ############## 
#### 11. INDIVIDUAL RIVER DAM IMPACTS ####
#### 11A. RESERVOIR CAPACITY VS. SEDIMENT FLUX, BY RIVER ####
# Table of each river: 
# Columns: one column for each decade's sediment flux
# One column for each decade's reservoir capacity
qss_2020_1984_w_dams <- dcast.data.table(
  na.omit(dam_ssc_comb[,.(tons_month_decadal = mean(tons_month, na.rm = T),
                          tons_month_decadal_sd = sd(tons_month, na.rm = T),
                          cap_mcm_cml = max(cap_mcm_cml, na.rm = T)),
                       by = .(RiverName, site_no, decade, Q_cms_avg, drainage_area_km2, Latitude, Longitude)],
          cols = c('Q_cms_avg','tons_month_decadal','decade', 'cap_mcm_cml')),
  RiverName + site_no + Q_cms_avg + drainage_area_km2 + Latitude + Longitude ~ decade,
  value.var = c('tons_month_decadal','tons_month_decadal_sd','cap_mcm_cml'))[
    ,':='(dams = ifelse(cap_mcm_cml_2015 == 0, 'No dams','Dams'))
  ]

# Calculate sediment flux relative change
# Also calculate dam capacity change relative to annual discharge
# Use 1984-1990 vs 2015-2020 if exists
# Otherwise, use 1990-1995 or 1995-2000
qss_2020_1984_w_dams <- qss_2020_1984_w_dams[,':='(
  tons_month_change = ifelse(!is.na(tons_month_decadal_1990),
                             (tons_month_decadal_2015 - tons_month_decadal_1990)/tons_month_decadal_1990,
                             ifelse(!is.na(tons_month_decadal_1995),
                                    (tons_month_decadal_2015 - tons_month_decadal_1995)/tons_month_decadal_1995,
                                    (tons_month_decadal_2015 - tons_month_decadal_2000)/tons_month_decadal_2000)),
  dam_cap_change = ifelse(!is.na(tons_month_decadal_1990),
                          (cap_mcm_cml_2015 - cap_mcm_cml_1990)/(Q_cms_avg * 3600 * 24 * 365.25/1e6) * 100,
                          ifelse(!is.na(tons_month_decadal_1995),
                                 (cap_mcm_cml_2015 - cap_mcm_cml_1995)/(Q_cms_avg * 3600 * 24 * 365.25/1e6) * 100,
                                 (cap_mcm_cml_2015 - cap_mcm_cml_2000)/(Q_cms_avg * 3600 * 24 * 365.25/1e6) * 100)))]


# Add pre-dam sediment flux estimate
qss_2020_1984_w_dams <- qss_2020_1984_w_dams[
  Qss_pre_byRiver[,.(RiverName, Continent_Region, tons_yr_preDam, tons_yr, Q_cms, Q_cms_preDam)],
  on = 'RiverName'][
    ,':='(tons_month_change_preDam = (tons_month_decadal_2015*12 - tons_yr_preDam)/tons_yr_preDam)
  ]

# Calculate Reservoir Capacity Index (RCI) for 2015
qss_2020_1984_w_dams <- qss_2020_1984_w_dams[,':='(RCI = cap_mcm_cml_2015/(Q_cms_avg * 3600 * 24 * 365.25/1e6) * 100)]

# Plot change in sediment flux relative to new reservoir capacity
# Since 1984
qss_2020_1984_river_change_plot <- ggpubr::ggarrange(
  ggplot(qss_2020_1984_w_dams[dam_cap_change == 0], aes(x = 'No dams', y = tons_month_change)) +
    geom_boxplot() +
    scale_y_continuous(limits = c(-1, 5)) +
    season_facet +
    labs(x = '',
         y = 'Sediment flux change\n(relative to 1990 baseline)'),
  ggplot(qss_2020_1984_w_dams[dam_cap_change > 0], aes(x = dam_cap_change, y = tons_month_change)) +
    geom_point(size = 2) +
    geom_smooth(method = 'lm', lty = 'dashed', color = 'black', formula = y~x, se = F) +
    scale_x_log10(breaks = c(0.01, 1, 100), labels = c(0.01, 1, 100)) +
    scale_y_continuous(limits = c(-1, 5)) +
    season_facet +
    theme(axis.title.y = element_blank()) +
    labs(x = 'Reservoir capacity as % of annual discharge'),
  nrow = 1, align = 'h', widths = c(0.25, 1))

ggsave(qss_2020_1984_river_change_plot, filename = paste0(wd_figures, 'qss_2020_1984_river_change_plot.pdf'),
       width = 7, height = 4, useDingbats = F)
ggsave(qss_2020_1984_river_change_plot, filename = paste0(wd_figures, 'qss_2020_1984_river_change_plot.png'),
       width = 7, height = 4)

## Fig. S10B
# Plot change in sediment flux relative to new reservoir capacity
# Relative to pre-dam estimate
# Left side: boxplot of rivers without dams
change_since_predam_no_dams <- ggplot(qss_2020_1984_w_dams[cap_mcm_cml_2015 == 0], aes(x = 'No dams', y = tons_month_change_preDam)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-1, 5)) +
  season_facet +
  labs(x = '',
       y = 'Sediment flux change\n(relative to pre-dam baseline)')
# Right side: scatterplot of rivers with major dams
change_since_predam_vs_RCI_with_dams <- ggplot(qss_2020_1984_w_dams[cap_mcm_cml_2015 > 0], 
       aes(x = cap_mcm_cml_2015/(Q_cms_avg * 3600 * 24 * 365.25/1e6) * 100, 
           y = tons_month_change_preDam)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm', lty = 'dashed', color = 'black', formula = y~x, se = F) +
  scale_x_log10(breaks = c(0.01, 1, 100), labels = c(0.01, 1, 100)) +
  scale_y_continuous(limits = c(-1, 5)) +
  season_facet +
  labs(x = 'Reservoir capacity as % of annual discharge',
       y = 'Sediment flux change\n(relative to pre-dam baseline)')


qss_2020_preDam_river_change_plot <- ggpubr::ggarrange(
  change_since_predam_no_dams,
  change_since_predam_vs_RCI_with_dams +
    theme(axis.title.y = element_blank()),
  nrow = 1, align = 'h', widths = c(0.25, 1))

ggsave(qss_2020_preDam_river_change_plot, filename = paste0(wd_figures, 'qss_2020_preDam_river_change_plot.pdf'),
       width = 7, height = 4, useDingbats = F)
ggsave(qss_2020_preDam_river_change_plot, filename = paste0(wd_figures, 'qss_2020_preDam_river_change_plot.png'),
       width = 7, height = 4)

#### 11B. CALCULATE SEDIMENT FLUX AND RCI ACROSS EPOCHS ####

# Merge sedimetnt flux by decade table 
# with trend test table for Kendall tests
qss_epochs_table <- qss_2020_1984_w_dams[
    ,':='(
      Mt_yr_decadal_1990 = tons_month_decadal_1990 * 12/1e6,
      Mt_yr_decadal_2000 = tons_month_decadal_2000 * 12/1e6,
      Mt_yr_decadal_2015 = tons_month_decadal_2015 * 12/1e6,
      Mt_yr_preDam = tons_yr_preDam/1e6,
      Mt_yr_reference = tons_yr/1e6,
      Q_cms_reference = Q_cms)
  ][,.(
    site_no, RiverName, Continent_Region, Latitude, Longitude, Q_cms_preDam, Q_cms_avg, dams, RCI,
    Mt_yr_decadal_1990, Mt_yr_decadal_2000, Mt_yr_decadal_2015, Mt_yr_preDam, Mt_yr_reference, Q_cms_reference)][
      trend_test_final_comb[,.(RiverName, site_no, trend_sign)],
      on = c('site_no','RiverName')
      ][!is.na(dams)]

# Calculate change across epochs. 
# Add columns for difference between sediment fluxes, dams
# Add a latitude group column for later partitioning
qss_epochs_table <- qss_epochs_table[
        ,':='(log10_Qss_difference = log10(Mt_yr_reference) - log10(Mt_yr_decadal_1990),
          log10_Q_cms_difference = log10(Q_cms_reference) - log10(Q_cms_avg),
          dams = ifelse(Mt_yr_reference != Mt_yr_preDam, 'Dams', dams),
          Mt_yr_decadal_1990 = ifelse(is.na(Mt_yr_decadal_1990), Mt_yr_decadal_2000, Mt_yr_decadal_1990),
          latitude_group = ifelse(abs(Latitude > 20), 'Mid- and High-Latitudes','Equatorial'))
  ]

median(sqrt(qss_epochs_table[dams == 'No dams']$log10_Qss_difference^2), na.rm = T)
median(sqrt(qss_epochs_table[log10_Q_cms_difference < 0.1]$log10_Qss_difference^2), na.rm = T)
median(sqrt(qss_epochs_table$log10_Qss_difference^2), na.rm = T)

# Are trends significantly related to amount of deforestation in watershed?
# The answer is no: minimal relationship without more sophisticated breakdown.
# What about in cases of low reservoir storage
ssc_trend_vs_deforestation_low_dam <- lm(trend_tau_ssc~cumulative_deforestation_rel,
                                         trend_test_final_comb[cumulative_deforestation_rel < 100 & dam_storage_percent < 1])
# What about in cases of high reservoir storage
ssc_trend_vs_deforestation_med_dam <- lm(trend_tau_ssc~cumulative_deforestation_rel,
                                         trend_test_final_comb[cumulative_deforestation_rel < 100 & dam_storage_percent >= 1 & dam_storage_percent < 20])

summary(ssc_trend_vs_deforestation_low_dam)
summary(ssc_trend_vs_deforestation_med_dam)

## Fig. S10A
# Plot distribution of SSC trend vs. dam + deforestation combinations
deforestation_dams_ssc_effect <- ggplot(trend_test_final_comb[cumulative_deforestation_rel < 100], 
                                        aes(x = cut(cumulative_deforestation_rel, c(-1, 20, Inf), 
                                                    labels = c('<20%', '>20%')), 
                                            y = trend_tau_ssc,
                                            fill = cut(cumulative_deforestation_rel, c(-1, 20, Inf), 
                                                       labels = c('<20% Deforest', '>20% Deforest')))) +
  geom_boxplot(notch = T) + 
  scale_fill_manual(values = c(decrease_color, increase_color)) +
  facet_wrap(.~cut(dam_storage_percent, c(-1, 0, 10, Inf), labels = c('No major dams', '< 10% RCI', '> 10% RCI'))) +
  season_facet +
  labs(
    fill = 'Deforestation',
    x = 'Watershed deforestation',
    y = 'SSC trend (1984-2021)\n(Kendall Tau statistic)'
  )

ggsave(deforestation_dams_ssc_effect, filename = paste0(wd_figures, 'deforestation_dams_ssc_effect.png'),
       width = 4.5, height = 5)
ggsave(deforestation_dams_ssc_effect, filename = paste0(wd_figures, 'deforestation_dams_ssc_effect.pdf'),
       width = 4.5, height = 5, useDingbats = F)

## Fig. S10
# A: boxplots of deforestation and dams vs. change in SSC
# B: Scatter plot of change in sediment flux 
# (relative to pre-dam)
deforest_dam_effect_plots_combined <- 
  deforestation_dams_ssc_effect /
  change_since_predam_vs_RCI_with_dams + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(1, 0.6))

ggsave(deforest_dam_effect_plots_combined, filename = paste0(wd_figures, 'FigS10_deforest_dam_effect_plots_combined.png'),
       width = 5, height = 9)
ggsave(deforest_dam_effect_plots_combined, filename = paste0(wd_figures, 'FigS10_deforest_dam_effect_plots_combined.pdf'),
       width = 5, height = 9, useDingbats = F)

#### 11C. SEDIMENT FLUX VS. DISCHARGE PLOT ####
## Fig. S13
# Sediment flux vs. discharge for each river
# Colored by latitude zone
Qss_Qcms_rating_curve <- ggplot(qss_epochs_table[
  !(RiverName %chin% c('Nueces'))
], aes(x = Q_cms_avg, y = Mt_yr_decadal_2015)) +
  geom_point(pch = 21, color = 'black', size = 2, stroke = 0.2, 
             # aes(fill = abs(log10_Q_cms_difference))) +
             aes(fill = Latitude)) +
  geom_smooth(aes(group = latitude_group,
                  color = latitude_group),
              lty = 'dashed', se = F, method = 'lm') +
  season_facet + 
  # scale_x_log10(labels = fancy_scientific_modified, limits = c(0.001, 2000)) +
  # scale_y_log10(labels = fancy_scientific_modified, limits = c(0.001, 2000)) +
  scale_color_manual(values = c('Equatorial' = '#D6AE42', 'Mid- and High-Latitudes' = 'grey40')) +
  scale_fill_gradientn(colors = c('#368096', '#4EB6D7','#D6AE42','grey70', 'white'), 
                       limits = c(-60,60), oob = squish) +
  scale_x_log10(labels = fancy_scientific_modified, 
                breaks = breaks, minor_breaks = NULL) +
  scale_y_log10(labels = fancy_scientific_modified,
                breaks = breaks, minor_breaks = NULL) + 
  annotation_logticks(sides = 'trbl') +
  coord_cartesian(clip = 'off') +
  theme(
    axis.ticks = element_line(size = 0.25),
    legend.position = 'right'
  ) +
  labs(
    x = 'Avg. Discharge (m<sup>3</sup>/s)',
    y = 'Sediment Flux (Mt/yr) from satellite estimate',
    color = 'Latitude zone'
  ) +
  theme(legend.title = element_markdown(),
        axis.title.x = element_markdown())

ggsave(Qss_Qcms_rating_curve, filename = paste0(wd_figures, 'FigS13_Qss_Qcms_rating_curve.png'),
       width = 6, height = 6)
ggsave(Qss_Qcms_rating_curve, filename = paste0(wd_figures, 'FigS13_Qss_Qcms_rating_curve.pdf'),
       width = 6, height = 6, useDingbats = F)

# Same rating curve plot, but:
# --broken out by latitude group
# --colored by reservoir capacity (RCI)
# (not used in paper)
Qss_Qcms_rating_curve_dammed <- ggplot(qss_epochs_table, 
                                       aes(x = Q_cms_avg, y = Mt_yr_decadal_2015)) +
  geom_point(pch = 21, color = 'black', size = 2, stroke = 0.2, 
             # aes(fill = abs(log10_Q_cms_difference))) +
             aes(fill = RCI + 1)) +
  geom_smooth(aes(group = latitude_group,
                  color = latitude_group),
              lty = 'dashed', se = F, method = 'lm') +
  season_facet + 
  facet_wrap(.~latitude_group) +
  scale_color_manual(values = c('Equatorial' = '#D6AE42', 'Mid- and High-Latitudes' = 'grey40')) +
  scale_fill_gradientn(colors = c('#368096', 'white'), 
                       trans = 'log10', limits = c(1, 100), oob = squish) +
  scale_x_log10(labels = fancy_scientific_modified, 
                breaks = breaks, minor_breaks = NULL) +
  scale_y_log10(labels = fancy_scientific_modified,
                breaks = breaks, minor_breaks = NULL) +
  annotation_logticks(sides = 'trbl') +
  coord_cartesian(clip = 'off') +
  theme(
    axis.ticks = element_line(size = 0.25),
    legend.position = 'right'
  ) +
  labs(
    x = 'Discharge (m3/s)',
    y = 'Sediment Flux (Mt/yr) from satellite estimate',
    color = 'Latitude zone'
  )

ggsave(Qss_Qcms_rating_curve_dammed, filename = paste0(wd_figures, 'Qss_Qcms_rating_curve_dammed.png'),
       width = 6, height = 6)
ggsave(Qss_Qcms_rating_curve_dammed, filename = paste0(wd_figures, 'Qss_Qcms_rating_curve_dammed.pdf'),
       width = 6, height = 6, useDingbats = F)

#### 12. COMPARE SEDIMENT FLUX TO REFERENCE, PRE-DAM DATA ####

# Legend data.table
legend_dt <- data.table(Mt_yr_preDam = 2, Mt_yr_decadal_1990 = 300, 
                        Mt_yr_reference = 0.1)
legend_dt_text <- data.table(text = c('Post-dam or\nno dam', 'Pre-dam'),
                             Mt_yr_decadal_1990 = c(400, 400), 
                             Mt_yr_reference = c(0.1,2))
## Fig S5A
# Compare satellite-estimates of sediment flux to reference dataset
# Show pre-dam flux as grey circles on same plot
# Arrow points to post-dam flux
qss_reference_compare <- ggplot(qss_epochs_table[
  !(RiverName %chin% c('Nueces')) &
    Mt_yr_decadal_1990 > 0.05
], aes(x = Mt_yr_reference, y = Mt_yr_decadal_1990)) +
  geom_point(pch = 21, stroke = 0.2, size = 2.5, color = 'black', aes(x = Mt_yr_preDam), fill = 'grey80') +
  geom_point(pch = 21, stroke = 0.2, size = 2.5, color = 'black', fill = 'grey20') +
  geom_point(data = legend_dt, pch = 21, stroke = 0.2, size = 2.5, color = 'black', aes(x = Mt_yr_preDam), fill = 'grey80') +
  geom_point(data = legend_dt, pch = 21, stroke = 0.2, size = 2.5, color = 'black', fill = 'grey20') +
  geom_segment(data = qss_epochs_table[Mt_yr_preDam > 0 & Mt_yr_preDam != Mt_yr_reference & Mt_yr_decadal_1990 > 0.05],
               aes(x = Mt_yr_preDam, xend = 10^(log10(Mt_yr_reference) + 0.2),
                   y = Mt_yr_decadal_1990, yend = Mt_yr_decadal_1990),
               arrow = arrow(type = 'closed', length = unit(4, 'points')),
               lwd = 0.2, linejoin = 'mitre',position = position_nudge(x = -0.02)) +
  geom_segment(data = legend_dt,
               aes(x = Mt_yr_preDam, xend = 10^(log10(Mt_yr_reference) + 0.2),
                   y = Mt_yr_decadal_1990, yend = Mt_yr_decadal_1990),
               arrow = arrow(type = 'closed', length = unit(4, 'points')),
               lwd = 0.2, linejoin = 'mitre',position = position_nudge(x = -0.02)) +
  geom_text(data = legend_dt_text, aes(label = text), vjust = 0, lineheight = 0.8) +
  geom_abline(slope = 1, intercept = 0) +
  season_facet + 
  scale_color_manual(values = c('Increase' = increase_color, 'Decrease' = decrease_color,
                                'No change' = 'grey60')) +
  scale_x_log10(labels = fancy_scientific_modified, limits = c(0.02, 2000), 
                breaks = breaks, minor_breaks = NULL) +
  scale_y_log10(labels = fancy_scientific_modified, limits = c(0.02, 2000),
                breaks = breaks, minor_breaks = NULL) +
  annotation_logticks(sides = 'trbl') +
  coord_cartesian(clip = 'off') +
  theme(
    axis.ticks = element_line(size = 0.25)
  ) +
  labs(
    x = 'Sediment Flux (Mt/yr) from Milliman and Farnsworth, 2012',
    y = 'Sediment Flux (Mt/yr) from satellite estimate'
  )

# Add a column for increasing and decreasing rivers
qss_epochs_table <- qss_epochs_table[
  ,':='(trend_sign_label = ifelse(trend_sign == 'Increase', 'Incr.', 
                                  ifelse(trend_sign == 'Decrease','Decr.','Not signif.')))
]

# Table of legend positions for arrow annotation
legend_dt_text_dammed <- data.table(Mt_yr_reference = c(0.05, 0.05, 0.05), 
                                    Mt_yr_decadal_1990 = 10^c(2.1, 2.2, 2.9),
                                    Mt_yr_decadal_2015 = 10^c(1.8, 2.5, 2.6),
                                    text_position_y = 10^c(1.95, 2.35, 2.75),
                                    dams = c('Dams', 'Dams','Dams'),
                                    trend_sign = c('Decrease','Increase','No change'))

## Fig. S5B 
# Panel plots of sediment flux compared to reference
# Upper panel: rivers with dams
# Lower panel: rivers without major dams
# Colored by whether flux has increased or decreased
qss_reference_compare_dammed <- ggplot(qss_epochs_table[
  !(RiverName %chin% c('Nueces')) &
    Mt_yr_decadal_1990 > 0.05 & 
    !is.na(dams)], 
  aes(x = Mt_yr_reference, y = Mt_yr_decadal_1990)) +
  geom_segment(aes(y = Mt_yr_decadal_1990, yend = Mt_yr_decadal_2015, xend = Mt_yr_reference,
                   color = trend_sign),
               arrow = arrow(type = 'closed', length = unit(4, 'points'))) +
  geom_segment(data = legend_dt_text_dammed, 
               aes(y = Mt_yr_decadal_1990, yend = Mt_yr_decadal_2015, xend = Mt_yr_reference,
                   color = trend_sign),
               arrow = arrow(type = 'closed', length = unit(4, 'points'))) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(data = legend_dt_text_dammed, 
            aes(x = 10^(log10(Mt_yr_reference)+0.15), y = text_position_y, label = trend_sign), 
            hjust = 0, vjust = 0.5) +
  season_facet + 
  facet_wrap(.~paste0('Post-1984 change, rivers with ', casefold(dams)), nrow = 2) +
  scale_color_manual(values = c('Increase' = increase_color, 'Decrease' = decrease_color,
                                'No change' = 'grey60')) +
  scale_x_log10(labels = fancy_scientific_modified, limits = c(0.02, 2000), 
                breaks = breaks, minor_breaks = NULL) +
  scale_y_log10(labels = fancy_scientific_modified, limits = c(0.02, 2000),
                breaks = breaks, minor_breaks = NULL) +
  annotation_logticks(sides = 'trbl') +
  coord_cartesian(clip = 'off') +
  theme(
    axis.ticks = element_line(size = 0.25)
  ) +
  labs(
    x = 'Sediment Flux (Mt/yr) from Milliman and Farnsworth, 2012',
    y = 'Sediment Flux (Mt/yr) from satellite estimate',
    color = ''
  )

# Add labels to reference plot
Qss_label_plot <- ggplot() + 
  geom_blank() +
  season_facet +
  theme_transparent() +
  theme(axis.title.y = element_text()) +
  labs(y = 'Sediment Flux (Mt/yr) from satellite estimate')

Qss_reference_label_plot <- ggplot() + 
  geom_blank() +
  season_facet +
  theme_transparent() +
  theme(axis.title.x = element_text()) +
  labs(x = 'Sediment Flux (Mt/yr) from Milliman and Farnsworth, 2012')

qss_reference_compare_panel_plot <- ggpubr::ggarrange(
  qss_reference_compare + labs(x = '',y = ''),
  qss_reference_compare_dammed + labs(x = '', y = ''),
  align = 'hv', labels = c('A','B'), widths = c(1, 0.8))

## Fig. S5
# Combine S5 plots
# A: Sediment flux vs. reference, with pre-dam annotated
# B: Sediment flux vs. reference, paneled by pre-dam, 
# (colored by direction of change)
qss_reference_compare_panel_plot_wlabels <- ggpubr::ggarrange(ggpubr::ggarrange(
  Qss_label_plot, NULL, qss_reference_compare_panel_plot, 
  widths = c(0.06, -0.03, 1), nrow = 1),
  NULL,
  Qss_reference_label_plot, nrow = 3, heights = c(1,-0.04, 0.06))

ggsave(qss_reference_compare_panel_plot_wlabels, 
       filename = paste0(wd_figures, 'FigS5_qss_reference_compare_panel_plot.png'),
       width = 8, height = 6) 
ggsave(qss_reference_compare_panel_plot_wlabels, 
       filename = paste0(wd_figures, 'FigS5_qss_reference_compare_panel_plot.pdf'),
       width = 8, height = 6, useDingbats = F) 

## CALCULATION FOR PAPER ##
# Reduction in sediment flux relative to pre-dam is significantly related to reservoir capacity
summary(qss_2020_preDam_river_change_lm) 

qss_2020_preDam_river_change_lm <- lm(
  tons_month_change_preDam ~ log10(RCI),
  qss_2020_1984_w_dams[RCI > 0])

# Reduction in sediment flux relative to 1984 is marginally related to new reservoir capacity
post_1984_rci_sed_reduction_model <- lm(tons_month_change~log10(dam_cap_change), 
                                        data = qss_2020_1984_w_dams[dam_cap_change > 0])

summary(post_1984_rci_sed_reduction_model) 


#### vvv WORKING HERE vvv ####
#### MOVE THIS BACK TO 8_FIGURES_AND_FINAL_CALCULATIONS ####
#### NEED TO CREATE precip_annual_avg ####
#### 13. THREE-WAY DOUBLE MASS ANALYSIS ####
#### 13A. COMBINE PRECIP, DISCHARGE, SEDIMENT, DAMS, DEFOREST ####
# Combine sediment and climate
river_Q_Qs_precip <- na.omit(Qss_river_annual[watershed_climate[
  ,.(year, RiverName, precip_m_yr, drainage_area_km2)], on = c('RiverName', 'year')][ 
    ,':='(Continent_Region = ifelse(Continent_Region %in% c('Europe','Eurasia'),
                                    'Europe/Eurasia', Continent_Region))
  ], cols = 'RiverName')

# Add sediment trend
river_Q_Qs_precip_deforest <- river_Q_Qs_precip[na.omit(trend_test_final_comb[
  ,.(RiverName, trend_sign, trend_tau, dam_storage_percent, cumulative_deforestation_rel)
], cols = 'RiverName'), on = 'RiverName']


# Plot increase at increasing rivers
SSC_anomaly_for_increasing_rivers <- ggplot(river_Q_Qs_precip_deforest[trend_sign == 'Increase'][
  ,':='(SSC_anomaly = SSC_mgL/mean(SSC_mgL, na.rm = T),
        Precip_anomaly = precip_m_yr/mean(precip_m_yr, na.rm = T)),
  by = .(RiverName, Continent_Region)
][,.(Precip_anomaly_mean = mean(Precip_anomaly, na.rm = T),
     SSC_anomaly_mean = mean(SSC_anomaly, na.rm = T)),
  by = .(year - year%%3, Continent_Region)],
aes(x = year, y = SSC_anomaly_mean)) + 
  geom_step() +
  season_facet +
  facet_wrap(.~Continent_Region)

# Calculate double-mass on discharge, sediment, and precip
Q_Qs_precip_triple_double_mass <- river_Q_Qs_precip[
  ,.(SSC_mgL = cumsum(SSC_mgL/mean(SSC_mgL, na.rm = T)),
     Q_cms = cumsum(Q_cms/mean(Q_cms, na.rm = T)),
     precip = cumsum(precip_m_yr/mean(precip_m_yr, na.rm = T)),
     year = year),
  by = .(RiverName, Continent_Region)
]

#### 13B. PLOT THREE-WAY DOUBLE-MASS PLOTS ####
# Set up double mass plot
# Most rapid increasing and decreasing rivers
double_mass_rivers <- c('Zhujiang','Mississippi', 'MacKenzie', 'Pindaré', 'Changjiang','Chorokhi', 'Purari', 'Rokan',
                        'Sittang', 'Song Hong', 'Suriname', 'Essequibo', 'Tarsus', 'Tugela', 'Turiaço', 'Yesil Irmak',
                        'Coco', 'Corubal', 'Curua Una', 'Guadiana', 'Huaihe', 'Huanghe',
                        'Inderagiri', 'Indus', 'Itata', 'Itajai-Acu', 'Itapecuru',
                        'Kajan', 'Kelantan', 'Kizil Irmak', 'Licungo', 'Mahakam', 'Mazaroni',
                        'Mendawai', 'Murray-Darling')
double_mass_labels <- na.omit(trend_test_final_comb[
  ,':='(cap_mcm_cml_max = cap_mcm_cml,
        dam_storage_percent_max = dam_storage_percent,
        cumulative_deforestation_rel_max = cumulative_deforestation_rel)][
          ,.(RiverName, trend_sign, trend_tau, cap_mcm_cml_max, dam_storage_percent_max, cumulative_deforestation_rel_max)
        ], cols = 'RiverName')

# Add labels to table
Q_Qs_precip_triple_double_mass <- Q_Qs_precip_triple_double_mass[double_mass_labels,
                                                                 on = c('RiverName')]

# Add dam capacity
Q_Qs_precip_triple_double_mass <- Q_Qs_precip_triple_double_mass[annual_cap_mcm_comb,
                                                                 on = c('RiverName', 'year')][
                                                                   , ':='(rci_normalized = cap_mcm_cml/cap_mcm_cml_max),
                                                                   by = .(RiverName)]
# Add deforestation
Q_Qs_precip_triple_double_mass <- merge(Q_Qs_precip_triple_double_mass,
                                        watershed_deforestation[,.(RiverName, year, cumulative_deforestation_rel)],
                                        by = c('RiverName', 'year'), all.x = T)

## Fig. S9
# Double-mass plots (normalized)
# Plot Q anomaly against SSC anomaly & precip anomaly
# A: Rivers with SSC Increase
# B: Rivers with SSC Decrease
for(i in 1:2){
  increase_decrease <- c('Increase', 'Decrease')[i]
  double_mass_labels_sel <- double_mass_labels[
    # RiverName %chin% double_mass_rivers & 
    trend_sign == increase_decrease][
      order(-abs(trend_tau))
    ][RiverName != 'Limpopo'][1:9]
  
  dt_sel <- Q_Qs_precip_triple_double_mass[
    RiverName %chin% double_mass_labels_sel$RiverName & 
      trend_sign == increase_decrease][
        ,':='(RiverName = ifelse(RiverName == 'Changjiang', 'Yangtze', RiverName))
      ]
  continent_labels <- dt_sel[,.(N = .N,
                                RCI = max(dam_storage_percent_max, na.rm = T),
                                Deforest_percent = max(cumulative_deforestation_rel, na.rm = T)), 
                             by = .(RiverName, Continent_Region)][
                               ,':='(Continent_label = factor(Continent_Region, 
                                                              levels = c('Asia','Europe/Eurasia','North America',
                                                                         'South America','Africa', 'Oceania'),
                                                              labels = c('AS','EU','NA','SA','AF','OC')))
                             ]
  double_mass_plots <- ggplot(
    dt_sel, 
    aes(x = Q_cms/37, y = SSC_mgL/37)) + 
    # geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    geom_line(color = increase_color) + 
    geom_line(aes(y = precip/37), color = decrease_color) + 
    geom_line(aes(y = rci_normalized)) +
    geom_line(aes(y = cumulative_deforestation_rel/cumulative_deforestation_rel_max), color = 'forest green') + 
    geom_text(data = continent_labels, aes(x = 1/37, y = Inf, label = Continent_label),
              hjust = 0, vjust = 1.5) +
    geom_text(data = continent_labels, aes(x = 1, y = -Inf, label = paste0('RCI: ', round(RCI, 0))),
              hjust = 1, vjust = -0.5, color = 'blue') +
    geom_text(data = continent_labels, aes(x = 1, y = 0.2, label = paste0('Defor: ', round(Deforest_percent, 0))),
              hjust = 1, vjust = -0.5, color = 'forest green') +
    # geom_vline(xintercept = yangtze_annual_double_mass[year == 2003]$Q_cms_yr) +
    season_facet +
    facet_wrap(.~RiverName, scales = 'free') + 
    scale_x_continuous(breaks = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 1)) +
    labs(
      x = 'Discharge',
      y = paste0("<span style = 'color:#FA9E2F'>**SSC**</span>",
                 ", <span style = 'color:#0076AD'>**Precipitation**</span>",
                 ", <span style = 'color:black'>**RCI**</span>",
                 ", and <span style = 'color:#228b22'>**Deforestation**</span>"),
      title = paste0('**Rivers with rapid SSC ', casefold(increase_decrease),'**')
    ) + 
    theme(axis.title.y = element_markdown(),
          plot.title = element_markdown(),
          panel.background = element_rect(fill = 'grey95'))
  
  if(i == 1){
    double_mass_plots_all <- double_mass_plots
  }else{
    double_mass_plots_all <- double_mass_plots_all / double_mass_plots + 
      plot_annotation(tag_levels = 'A')
  }
  
}

ggsave(double_mass_plots_all, filename = paste0(wd_figures, 'FigS9_double_mass_plots_all.png'),
       width = 5, height = 8.5)
ggsave(double_mass_plots_all, filename = paste0(wd_figures, 'FigS9_double_mass_plots_all.pdf'),
       width = 5, height = 8.5, useDingbats = F)

print(double_mass_plots_all)

yangtze_double_mass_plot <- ggplot(
  yangtze_annual_double_mass, aes(x = Q_cms_yr, y = tons_yr/1e6)) + 
  geom_point() +
  geom_line() + 
  geom_vline(xintercept = yangtze_annual_double_mass[year == 2003]$Q_cms_yr) +
  season_facet


#### 14.EXTRA WHEEL PLOTS ####
#### 14A. EXAMPLE WHEEL PLOTS FOR EACH CONTINENT ####
# One example river per continent
# (Not used in paper)
wheel_plots_example_continental_6 <- ggplot(site_monthly_band_median_color[RiverName %in% c('Amazon',
                                                                                            'Mississippi',
                                                                                            'Mekong',
                                                                                            'Ikopa', 
                                                                                            'Danube',
                                                                                            'Burdekin')]) +
  geom_bar(width = 1, color = 'black', stat = 'identity', lwd = 0.25,
           aes(x = month, y = Q_cms_norm, fill = rgb(B3/2300, B2/2300, B1/2300))) +
  season_facet +
  facet_wrap(dir = 'v', factor(Continent_Region, 
                               levels = c('North America','South America',
                                          'Europe/Eurasia','Africa','Asia','Oceania'))~
               RiverName, nrow = 2) +
  scale_fill_identity() +
  scale_y_sqrt(limits = c(0,1.1)) +
  scale_x_continuous(breaks = c(1:12), labels = month_labels) +
  coord_polar() +
  theme_void()

ggsave(wheel_plots_example_continental_6, filename = paste0(wd_figures, 'wheel_plots_example_continental_6.pdf'),
       width = 6, height = 4, useDingbats = F)
ggsave(wheel_plots_example_continental_6, filename = paste0(wd_figures, 'wheel_plots_example_continental_6.png'),
       width = 6, height = 4)

# Three examples per continent
# (Not used in paper)
wheel_plots_example_continental <- ggplot(site_monthly_band_median_color[RiverName %in% c('Amazon','Santa Cruz','Paraná',
                                                                                   'Mississippi','MacKenzie','Yukon',
                                                                                   'Changjiang','Huanghe','Mekong',
                                                                                   'Ikopa', 'Orange','Senegal',
                                                                                   'Danube','Maas','Rhone',
                                                                                   'Murray-Darling','Burdekin','Fly')]) +
  geom_bar(width = 1, color = 'black', stat = 'identity', lwd = 0.25,
           aes(x = month, y = Q_cms_norm, fill = rgb(B3/2500, B2/2500, B1/2500))) +
  season_facet +
  facet_wrap(dir = 'v', Continent_Region~RiverName, nrow = 3) +
  scale_fill_identity() +
  scale_y_sqrt(limits = c(0,1.1)) +
  scale_x_continuous(breaks = c(1:12), labels = month_labels) +
  # theme(
  #   panel.grid.major = element_line(color = 'black')
  # ) +
  coord_polar() +
  theme_void()


#### 14B. ALL RIVER WHEELS, ORGANIZED BY LANDMASS, CLUSTER ####
# Grid plot of wheels for each continent, cluster
## Aesthetically better to left justify? a challenge...
for(i in 1:6){
  # Select a landmass from list of 6
  continent_sel <- continent_area_Q[i]
  continent_sel_region <- continent_sel$Continent_Region
  dt_color <- merge(
    na.omit(
      plot_data[Continent_Region == continent_sel_region], 
      cols = c('B1','B2','B3')),
    river_cluster_assignment_single[,.(site_no,cluster)],
    by = 'site_no')
  
  # Get unique clusters for landmass
  unique_clusters <- sort(unique(dt_color$cluster))
  continent_cluster_plots <- vector('list', length(unique_clusters))
  names(continent_cluster_plots) <- unique_clusters
  n_unique_clusters <- length(unique_clusters)
  # For each cluster, make a facet plot with a wheel for each river in that cluster
  for(j in 1:n_unique_clusters){
    cluster_sel <- unique_clusters[j]
    dt_color_sel <- dt_color[cluster == cluster_sel]
    n_sites_sel <- length(unique(dt_color_sel$site_no))
    # Set number of plot rows based on the number of sites in each category
    nrow_sel <- ceiling(n_sites_sel/6)
    
    # Make faceted plot -- each river gets a panel
    continent_river_color_cluster_grid_plot <- ggplot(dt_color_sel) +
      geom_bar(width = 1, color = 'black', stat = 'identity', lwd = 0.2,
               aes(x = month, y = Q_cms_norm, fill = rgb(B3/2100, B2/2100, B1/2100))) +
      season_facet +
      facet_wrap(.~RiverName, nrow = nrow_sel) +
      scale_fill_identity() +
      scale_y_sqrt(limits = c(0,1.1)) +
      scale_x_continuous(breaks = c(1:12), labels = month_labels) +
      coord_polar() +
      theme_void() 
    # Combine plots for each cluster into all-landmass plot
    if(j == 1){
      continent_cluster_plots <- wrap_elements(
        set_panel_size(continent_river_color_cluster_grid_plot, 
                       width  = unit(1.15, "in"),
                       height = unit(1.15, "in")))
      nrow_total <- nrow_sel
    }else{
      continent_cluster_plots <- continent_cluster_plots / 
        wrap_elements(
          set_panel_size(continent_river_color_cluster_grid_plot, 
                         width  = unit(1.15, "in"),
                         height = unit(1.15, "in")))
      nrow_total <- c(nrow_total, nrow_sel)
    
    }
    # continent_cluster_plots[[j]] <- set_panel_size(continent_river_color_cluster_grid_plot, width  = unit(1.15, "in"),
    #                                                height = unit(1.15, "in"))
    
  }
  
  continent_cluster_plots <- continent_cluster_plots + plot_layout(heights = nrow_total/sum(nrow_total))
  # continent_river_color_grid_plot <- ggpubr::ggarrange(plotlist = continent_cluster_plots,
  #                                              title = continent_sel_region,
  #                                              ncol = 1)
  
  continent_cluster_plots <- continent_cluster_plots + 
    plot_annotation(tag_levels = 1, tag_prefix = 'Cluster ', title = continent_sel_region) & 
    theme(plot.tag = element_text(size = 11, face = 'bold'),
          plot.title = element_text(size = 14, face = 'bold'))
  # continent_river_color_grid_plot <- grid.arrange(grobs = continent_cluster_plots, ncol = 1, nrow = length(unique_clusters))
  ggsave(continent_cluster_plots, 
         filename = paste0(wd_figures, 'continent_river_color_grid_plot_', gsub('/','_',continent_sel_region), '.png'),
         width = 8, height = sum(nrow_total) * 1.6)
  ggsave(continent_cluster_plots, 
         filename = paste0(wd_figures, 'continent_river_color_grid_plot_', gsub('/','_',continent_sel_region), '.pdf'),
         width = 8, height = sum(nrow_total) * 1.75)
  
}



#### PLOTS AND ANALYSIS NOT USED IN PAPER ####
#### 1. KENDALL TREND COUNT OF RIVERS BY LATITUDE ####
# Count of trends in each direction, by latitude bin
trend_count_by_latitude <- na.omit(trend_test_final_comb, cols = 'trend_tau')[
  ,':='(Latitude_binned = Latitude - Latitude%%3,
        sign = ifelse(trend_tau > 0, 'Increasing','Decreasing'))][
          ,.(N_signed = .N),
          by = .(Latitude_binned, sign)
        ]

ggplot(trend_count_by_latitude, aes(x = Latitude_binned, y = N_signed, 
                                    group = sign,
                                    color = sign)) +
  geom_line() +
  season_facet +
  rotate()
#### 2. NEW RESERVOIR CAPACITY SINCE 1984, MAP ####

# Map of NEW reservoir capacity (since 1984)
# (Not used in paper)
watershed_dam_change_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = pre_post_dam, 
             stroke = 0.2, pch = 21,
             aes(x = Longitude, y = Latitude, 
                 fill = (area_km2_cml - area_km2_cml_1984)/drainage_area_km2 * 100,
                 size = (area_km2_cml - area_km2_cml_1984)/drainage_area_km2 * 100)) +
  geom_point(data = pre_post_dam[(area_km2_cml - area_km2_cml_1984) == 0], 
             stroke = 0.2, pch = 21,
             aes(x = Longitude, y = Latitude, 
                 size = (area_km2_cml - area_km2_cml_1984)/drainage_area_km2 * 100),
             fill = 'grey40') +
  scale_fill_gradientn(limits = c(0, 0.25), oob = squish,
                       colors = c('#2D2E7D','#7475B0','#ACCAE3','#B08874','#E6AC76')) +
  scale_size_continuous(range = c(2,10), breaks = c(100, 500, 1000)) +
  season_facet +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(-61, 86), expand = expansion(mult = c(0, 0))) +
  labs(
    fill = 'Change in\n% Watershed Reservoir Area ,\n2000-2020',
    size = '**Avg. Discharge**<br>**(100s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())
#### 3. YANGTZE RIVER DAMS AND SEDIMENT FLUX ####
#### 4. DEFOREST., DAM AREA CORRELATED WITH FLUX REDUCTION? ####
# Fraction of watershed area covered by reservoirs, deforestation
watershed_impacts <- pre_post_ssc[,.(RiverName, Latitude, Longitude, tons_month_2015_1990,Q_cms_avg)][
  pre_post_deforestation[,.(RiverName, cumulative_deforestation_rel)], on = 'RiverName'][
    pre_post_dam[,.(RiverName, drainage_area_km2, area_km2_cml)], on = 'RiverName'
  ]

# Does it correlate with reduction in sediment flux?
# Answer is weakly
# Not a robust evaluation
ggplot(watershed_impacts[area_km2_cml > 0], aes(x = log10(area_km2_cml/drainage_area_km2), y = tons_month_2015_1990)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  season_facet

# Does deforestation correlate?
# No.
ggplot(watershed_impacts[area_km2_cml > 0], aes(x = log10(cumulative_deforestation_rel), y = tons_month_2015_1990)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  season_facet
#### 5. WHEEL PLOTS FOR EVERY RIVER, MAX SSC COLOR ARROW MAP ####
# Replace high values in color data to allow for plotting of clearer colors
plot_data <- site_monthly_band_median_color[,':='(
  B1 = ifelse(B1 >= 2100, 2099, B1),
  B2 = ifelse(B2 >= 2100, 2099, B2),
  B3 = ifelse(B3 >= 2100, 2099, B3),
  B4 = ifelse(B4 >= 2100, 2099, B4)
)]

# Function to make a wheel plot given input data
# Input data needs to have columns of B1, B2, B3 (for blue, green, red)
# In reflectances scaled by 10,000
# e.g., B1 = fractional reflectance in blue band * 10000
wheel_plot_fun <- function(data){
  # NOTE: drop = F keeps every month even if missing
  ggplot(data) +
    geom_bar(width = 1, color = 'black', stat = 'identity', lwd = 0.2,
             aes(x = factor(month, levels = c(1:12)), 
                 y = Q_cms_norm, 
                 fill = rgb(B3/2100, B2/2100, B1/2100))) +
    # season_facet +
    # facet_wrap(dir = 'v', Continent_Region~RiverName, nrow = 3) +
    # scale_x_continuous(breaks = c(1:12), labels = month_labels, expand = expansion(mult = 0)) +
    scale_fill_identity(drop = F) +
    scale_y_sqrt(limits = c(0,1.1)) +
    scale_x_discrete(breaks = c(1:12), labels = month_labels, drop = F) + 
    # theme(
    #   panel.grid.major = element_line(color = 'black')
    # ) +
    coord_polar() +
    theme_void()
}

# TEST
# wheel_plot_fun(plot_data[RiverName == 'Alazeya'])

# Move plot to a particular location given by its latitude and longitude
annotation_fun <- function(data, Latitude_mil,Longitude_mil, plot_fun) {
  subplot = plot_fun(data)
  sub_grob <- annotation_custom(ggplotGrob(subplot), 
                                x = Longitude_mil-5, y = Latitude_mil-5, 
                                xmax = Longitude_mil+5, ymax = Latitude_mil+5)
}

# Remove NA values from by-site reflectance data
# Arrange them into subgrobs using longitude and latitude
subgrobs <- na.omit(plot_data, cols = c('B1','B2','B3', 'Latitude_mil','Longitude_mil')) %>% 
  nest(-Longitude_mil,-Latitude_mil)  %>%
  pmap(annotation_fun, plot_fun = wheel_plot_fun)

# Make a plot of the world map for wheel plots to get added to
p <- ggplot(data = Qss_peak_month) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'black', lwd = 0.2) +
  season_facet +
  # theme(legend.position = c(0.62, 0.08)) +
  # theme(panel.background = element_rect(fill = 'grey30')) +
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
    fill = '',
    size = ''
  )

monthly_color_map_plot <- p + subgrobs

ggsave(monthly_color_map_plot, filename = paste0(wd_figures, 'test_month_radial_plot.pdf'),
       width = 18, height = 10)

ggsave(monthly_color_map_plot, filename = paste0(wd_figures, 'monthly_color_map_plot.png'),
       width= 18, height = 10)

# Get color of highest QSS month
# For each river
peak_ssc_month <- na.omit(site_monthly_band_median_color[
  Qss_peak_month[,.(RiverName, month)], on = c('RiverName','month')],
  cols = c('Longitude_mil','Latitude_mil','B3','B2','B1'))[order(B3, decreasing = T)]

peak_ssc_month_color_map <- p +
  geom_point(data = peak_ssc_month, aes(x = Longitude_mil, y = Latitude_mil),
             size = 1, color = 'black') +
  geom_spoke(data = peak_ssc_month, 
             arrow = arrow(length = unit(0.12, 'cm'), type = 'closed'),
             lwd = 0.75, linejoin = 'mitre',
             aes(x = Longitude_mil, y = Latitude_mil, angle = -(month/12 * 3.14 * 2) + 3.14/2,
                 radius = Q_cms_norm*5, color = rgb(B3/2100, B2/2100, B1/2100))) +
  scale_color_identity()

ggsave(peak_ssc_month_color_map, 
       filename = paste0(wd_figures, 'peak_ssc_month_color_map.pdf'),
       width = 11, height = 6)

#### 6. EVERY RIVER PANEL PLOT SSC, DAM, DEFORESTATION TIMESERIES ####
# Make panel plot with every 4 rivers (alphabetical order)
# Plots show daily SSC variation, dam building, deforestation
for(i in 1:ceiling(nrow(ssc_rivers_stats)/4)){
  # Select four rivers
  # Calculate discharge and sediment flux monthly
  example_rivers <- ssc_rivers_stats[(i*4 - 3):(i*4)]$RiverName
  ssc_subset <- sed_transport_month_yr[RiverName %in% example_rivers & 
                                         SSC_mgL < 8000]
  
  ssc_monthly_subset <- na.omit(ssc_subset[
    ,.(tons_month = mean(tons_month, na.rm = T),
       Q_cms = mean(Q_cms, na.rm = T),
       SSC_mgL = mean(SSC_mgL, na.rm = T)), by = .(RiverName, decade, month)
  ], cols = 'decade')
  
  # Make a monthly plot of discharge and sediment flux
  monthly_q_qss_plot <- ggplot(ssc_monthly_subset,
                               aes(x = month, y = tons_month/1e6)) +
    geom_line(color = 'brown') +
    geom_line(aes(y = Q_cms/1000), color = 'blue') +
    facet_grid(RiverName~decade, scales = 'free_y') +
    season_facet +
    labs(x = 'Month',
         y = 'Mt/month\nAvg. Discharge (cms)')
  
  ggsave(monthly_q_qss_plot, width = 6, height = 10,
         filename = paste0(wd_exports,'monthly_Q_Qss_panel_', i,'.pdf'),
  ) 
  
  river_landsat_pred_subset <- 
    
    SSC_timeseries <- ggplot(ssc_subset,
                             aes(x = month/12 + year, y = SSC_mgL)) +
    # geom_line(aes(alpha = SSC_source, group = RiverName)) +
    geom_line(aes(color = SSC_source, group = RiverName)) +
    # geom_step(data = na.omit(dam_timeseries[
    #   RiverName %chin% example_rivers & year > 1983][
    #     ssc_rivers_stats[,.(site_no, RiverName, drainage_area_km2)], 
    #     on = c('RiverName','site_no')
    #   ]),
    #   aes(x = ymd(paste0(year,'-01-01')), 
    #       y = area_km2_cml/drainage_area_km2*10000*100),
    #       alpha = 0.5) +
    geom_smooth(color = 'orange', method = 'loess') +
    facet_wrap(.~RiverName, scales = 'free_y', nrow = 1) +
    season_facet +
    scale_x_continuous(limits = c(1984, 2021)) +
    # scale_alpha_manual(values = c('Rating curve estimate' = 0.2, 'Landsat observation' = 0.4)) +
    scale_color_manual(values = c('Rating curve estimate' = 'grey80', 'Landsat observation' = 'grey60')) +
    # scale_x_date(limits = c(ymd('1984-01-01', '2021-01-01'))) +
    # scale_y_continuous(
    #     
    #     # Features of the first axis
    #     name = "SSC (mg/L)",
    #     # Add a second axis and specify its features
    #     sec.axis = sec_axis(trans = ~./10000, 
    #                         name = 'Dam Storage (Percent watershed area)')
    #   ) +
    # theme(
    #   axis.title.y = element_text(color = "orange"),
  #   axis.title.y.right = element_text(color = "black"))
  labs(
    x = '',
    y = 'SSC (mg/L)'
  )
  
  
  dams_timeseries <- ggplot(
    na.omit(dam_timeseries[
      RiverName %chin% example_rivers & year > 1983][
        ssc_rivers_stats[,.(site_no, RiverName, drainage_area_km2)],
        on = c('RiverName','site_no')
      ]),
    aes(x = year,
        y = area_km2_cml/drainage_area_km2 * 100)) +
    geom_step() +
    geom_step(data = watershed_deforestation[
      RiverName %chin% example_rivers],
      aes(y = cumulative_deforestation_rel),
      color = 'forest green') +
    facet_wrap(.~RiverName, scales = 'free_y', nrow = 1) +
    season_facet +
    scale_x_continuous(limits = c(1984, 2021)) +
    scale_y_continuous(limits = c(0, NA),
                       sec.axis = sec_axis(~.*1, name='Deforested area\n(% watershed area)')) +
    labs(
      x = '',
      y = 'Reservoir area\n(% watershed area)'
    ) + 
    theme(axis.title.y.right = element_text(color = 'forest green'))
  
  SSC_dams_timeseries <- ggpubr::ggarrange(dams_timeseries, SSC_timeseries,
                                           nrow = 2, align = 'hv')
  # print(dam_cumulative_mcm_example_plot)
  ggsave(SSC_dams_timeseries, width = 9, height = 6,
         filename = paste0(wd_exports,'SSC_dams_ts_panel_', i,'.pdf'),
  ) 
}




#### 7. COMPARISON OF DISCHARGE FROM THIS STUDY TO REFERENCE ####
Q_cms_reference_compare <- ggplot(qss_epochs_table[
  !(RiverName %chin% c('Nueces'))
], aes(x = Q_cms_reference, y = Q_cms_avg)) +
  geom_point(pch = 21, color = 'black', size = 2, stroke = 0.2, 
             # aes(fill = abs(log10_Q_cms_difference))) +
             aes(fill = abs(log10_Qss_difference))) +
  geom_abline(slope = 1, intercept = 0) +
  season_facet + 
  # scale_x_log10(labels = fancy_scientific_modified, limits = c(0.001, 2000)) +
  # scale_y_log10(labels = fancy_scientific_modified, limits = c(0.001, 2000)) +
  scale_fill_gradientn(colors = c('royal blue', 'white'), limits = c(0, 0.75), oob = squish) +
  scale_x_log10(labels = fancy_scientific_modified, limits = c(10, 1e5), 
                breaks = breaks, minor_breaks = NULL) +
  scale_y_log10(labels = fancy_scientific_modified, limits = c(10, 1e5),
                breaks = breaks, minor_breaks = NULL) +
  annotation_logticks(sides = 'trbl') +
  coord_cartesian(clip = 'off') +
  theme(
    axis.ticks = element_line(size = 0.25)
  ) +
  labs(
    x = 'Avg. discharge (m<sup>3</sup>/s) from Milliman and Farnsworth, 2012',
    y = 'Avg. discharge (m<sup>3</sup>/s) estimated in this study'
  ) +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown()
  )
#### 8. NOT DONE! SSC TIMESERIES VS DAM BUILDING, EACH RIVER ####
river_landsat_daily <- river_landsat_pred[ssc_ann_monthly[, .(Q_cms,year,RiverName,site_no, month)], 
                                          on = c('year','RiverName','site_no', 'month')][
                                            ,':='(year = year(sample_dt))
                                          ][
                                            dam_timeseries, on = c('site_no','RiverName', 'year')
                                          ][
                                            ,':='(tons_month = SSC_mgL * Q_cms)
                                          ]


# Identify years with big dam area increases
for(i in 1:length(ssc_rivers$RiverName)){
  
  river_name_sel <- ssc_rivers$RiverName[i]
  river_landsat_daily_sel <- river_landsat_daily[
    RiverName ==  river_name_sel & SSC_mgL > 0 & SSC_mgL < 20000]
  first_yr <- min(river_landsat_daily_sel[SSC_mgL > 0]$year,na.rm = T)
  final_dam_area <- max(dam_timeseries[RiverName == river_name_sel]$area_km2_cml, na.rm = T)
  ssc_dam_diff <- dam_timeseries[RiverName == river_name_sel][
    ,':='(dam_area_increase_percent = area_km2/area_km2_cml*100)][
      (year == first_yr) | (dam_area_increase_percent > 5 & year > first_yr + 1),
      .(year, area_km2, area_km2_cml)][
        ,':='(dam_period = ifelse(year == first_yr & area_km2_cml < 0.05 * final_dam_area, 'Pre-dam',
                                  paste0('Dam period ', year)))
      ][,.(year, dam_period)]
  
  setkey(ssc_dam_diff, 'year')
  setkey(river_landsat_daily_sel, 'year')
  
  ssc_dam_periods <- ssc_dam_diff[river_landsat_daily_sel, roll = Inf
  ]
  
  ssc_dam_periods_summary <- ssc_dam_periods[
    ,.(SSC_mgL = mean(SSC_mgL,na.rm = T),
       period_start_dt = min(sample_dt, na.rm = T),
       period_end_dt = max(sample_dt, na.rm = T)),
    by = dam_period
  ]
  dam_timeseries_plot <- ggplot() +
    geom_line(data = ssc_dam_periods, 
              aes(x = sample_dt, y = SSC_mgL, color = dam_period,
                  group = dam_period)) +
    geom_linerange(data = ssc_dam_periods_summary,
                   aes(xmin = period_start_dt,
                       xmax = period_end_dt,
                       y = SSC_mgL,
                       color = dam_period)) + 
    # geom_smooth() +
    season_facet +
    labs(
      x = '',
      y = 'SSC (mg/L)',
      title = river_name_sel
    )
  
  ggsave(dam_timeseries_plot, 
         filename = paste0(wd_exports, 'Dam_SSC_timeseries_', 
                           gsub('/', '', river_name_sel), '.pdf'),
         width = 7, height = 5, useDingbats = F)
}




