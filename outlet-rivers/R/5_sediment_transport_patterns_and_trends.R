# Calculate patterns and trends in suspended sediment transport.
# For individual rivers over time and in various aggregations/subsets
# (landmass, pre-/post-dam, etc.)
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
#### 1A. IMPORT SUSPENDED SEDIMENT TRANSPORT DATA AND METADATA ####
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

#### 2. STATISTICS AND TIMESERIES OF RIVER SSC, QSS, AND Q ####
#### 2A. BY-RIVER AVERAGES ####
# Calculate averages for each river
months_with_observations_meta <- sed_transport_month_yr[SSC_source == 'Landsat observation',.(
  N = .N
), by = .(RiverName, site_no)]

# Calculate all-time averages by river
ssc_rivers_stats_full <- sed_transport_month_yr[
  ,.(tons_month_avg = mean(tons_month, na.rm = T),
     sd_comb = sd(SSC_mgL, na.rm = T)/mean(SSC_mgL, na.rm = T)/sqrt(.N),
     sd_indiv = sum(SSC_mgL_se_rel^2, na.rm = T)/.N,
     tons_month_se_rel = mean(SSC_mgL_se_rel, na.rm = T),
     Q_cms_avg = mean(Q_cms, na.rm = T),
     N_filled = .N),
  by = .(RiverName, ID, site_no, Continent_Region, Latitude, Longitude, drainage_area_km2, Climate_T, PrimGeo)][
    ,':='(SSC_mgL = tons_month_avg/Q_cms_avg * 0.391993,
          sd_total = sqrt(sd_comb^2 + sd_indiv^2)/sqrt(N_filled))
  ]

# Merge with metadata
ssc_rivers_stats_full <- merge(ssc_rivers_stats_full, 
                       months_with_observations_meta, 
                       by = c('RiverName','site_no'),
                       all.x = T)

# Make a simple table of average SSC, Q by site, all time
ssc_rivers_stats <- ssc_rivers_stats_full[,.(site_no, RiverName, ID, drainage_area_km2, Latitude, Longitude,
                                     N, Q_cms_avg, tons_month_avg
                                     )]

# Select just the River Name and site number for each river
ssc_rivers <- ssc_rivers_stats[,.(site_no, RiverName, ID)]

#### 2B. BY-RIVER ANNUAL, DECADAL, MONTHLY ####
# Avg. Qss for each river, by year
Qss_river_annual <- sed_transport_month_yr[
  ,.(tons_month_avg = mean(tons_month, na.rm = T),
     sd_comb = sd(SSC_mgL, na.rm = T)/mean(SSC_mgL, na.rm = T)/sqrt(.N),
     sd_indiv = sum(SSC_mgL_se_rel^2, na.rm = T)/.N,
     tons_month_se_rel = mean(SSC_mgL_se_rel, na.rm = T),
     Q_cms = mean(Q_cms, na.rm = T),
     N = .N),
  by = .(RiverName, ID, site_no, year, decade, Continent_Region, Latitude, Longitude, Climate_T, PrimGeo)][
    ,':='(SSC_mgL = tons_month_avg/Q_cms * 0.391993,
          sd_total = sqrt(sd_comb^2 + sd_indiv^2)/sqrt(N))
  ]

# Avg. SSC for each river, by year (direct observations only)
SSC_annual_avg <- sed_transport_month_yr[SSC_source == 'Landsat observation'][
  ,.(SSC_mgL = mean(SSC_mgL, na.rm = T),
     N = .N),
  by = .(RiverName, ID, site_no, year, decade, Continent_Region, Latitude, Longitude, Climate_T, PrimGeo)]


# Avg. Qss for each river, by decade
Qss_river_decadal <- sed_transport_month_yr[
  ,.(tons_month_avg = mean(tons_month, na.rm = T),
     sd_comb = sd(SSC_mgL, na.rm = T)/mean(SSC_mgL, na.rm = T)/sqrt(.N),
     sd_indiv = sum(SSC_mgL_se_rel^2, na.rm = T)/.N,
     tons_month_se_rel = mean(SSC_mgL_se_rel, na.rm = T),
     Q_cms = mean(Q_cms, na.rm = T),
     N = .N),
  by = .(RiverName, ID, site_no, decade, Continent_Region, Latitude, Longitude, Climate_T, PrimGeo)][
    ,':='(SSC_mgL = tons_month_avg/Q_cms * 0.391993,
          sd_total = sqrt(sd_comb^2 + sd_indiv^2)/sqrt(N))
  ]


# Avg. monthly Q, all time
avg_monthly_Q_Qss <- sed_transport_month_yr[,.(Q_cms = mean(Q_cms, na.rm = T),
                                               Qss_tons_month = mean(tons_month, na.rm = T),
                                               Qss_tons_month_sd = sd(tons_month, na.rm = T)),
                                        by = .(RiverName, ID, site_no, Continent_Region, month, Latitude_mil, Longitude_mil)][
                                          ,':='(Q_cms_norm = Q_cms/max(Q_cms, na.rm = T),
                                                Qss_norm = Qss_tons_month/max(Qss_tons_month, na.rm = T),
                                                Qss_rel_sd = sd(Qss_tons_month/max(Qss_tons_month, na.rm = T))),
                                          by = .(RiverName, ID, site_no, Continent_Region, Latitude_mil, Longitude_mil)
                                        ]


# Monthly peak discharge and which month
Qss_peak_month <- sed_transport_month_yr[
  ,.(tons_month_avg = mean(tons_month, na.rm = T),
     tons_month_se_rel_avg = mean(SSC_mgL_se_rel, na.rm = T),
     Q_cms = mean(Q_cms, na.rm = T),
     N = .N),
  by = .(RiverName, Continent_Region, Latitude, Longitude, month, Climate_T, PrimGeo)][
    ,':='(SSC_mgL = tons_month_avg)
  ]

# Which month has peak discharge
Qss_peak_month <- Qss_peak_month[Qss_peak_month[,.I[which.max(SSC_mgL)],
                                                by = RiverName]$V1]

#### 2C. BY-CONTINENT PRE-DAM, POST-DAM, BY-DECADE ####
# Calculate continent area and annual discharge
# (TO DO: calculate annual precipitation here too)
continent_area_Q <- na.omit(ssc_rivers_stats[,.(RiverName, ID, site_no, drainage_area_km2)][
  ssc_rivers_stats_full, on = c('RiverName', 'site_no', 'ID')
][
  ,':='(Continent_Region = ifelse(Continent_Region %in% c('Europe','Eurasia'),'Europe/Eurasia', Continent_Region))][
    ,.(drainage_area_km2 = sum(drainage_area_km2, na.rm = T),
       Q_km3_yr = sum(Q_cms_avg, na.rm = T) * 0.0315569),
    by = Continent_Region])[
      Continent_Region != ''
    ]

# Final river names (with Landsat and Milliman data)
final_river_IDs <- unique(sed_transport_month_yr$ID)
# For pre data, select Landsat rivers from Milliman database
# Modify columns to match Landsat analysis output
Qss_pre_byRiver <- Q_stn_metadata[
  ID %chin% final_river_IDs][
    ,':='(
      decade = 'Baseline',
      year = 'Baseline',
      Q_cms = Q * 31.6888,
      Q_cms_preDam = PreDam_Q * 31.6888,
      tons_yr = TSS * 1e6,
      tons_yr_preDam = ifelse(!is.na(PD_TSS), PD_TSS * 1e6, TSS*1e6))][
        ,.(RiverName, Continent_Region, decade, Q_cms, Q_cms_preDam, tons_yr, tons_yr_preDam)]

# Number of rivers with pre-dam estimate
n_rivers_with_pre_dam_estimate <- nrow(Qss_pre_byRiver[!is.na(tons_yr_preDam)])

# Sum continent by decade
# Propagate uncertainty
Qss_global_annual <- na.omit(sed_transport_month_yr[
  ,':='(Continent_Region = ifelse(Continent_Region %in% c('Europe','Eurasia'),'Europe/Eurasia', Continent_Region))][
    ,.(tons_month_avg = mean(tons_month, na.rm = T),
       tons_month_se_rel_avg = mean(SSC_mgL_se_rel, na.rm = T),
       N = .N),
    by = .(RiverName, Continent_Region, year, decade, month)][
      ,.(tons_yr = sum(tons_month_avg, na.rm = T),
         sd_indiv = sum(tons_month_se_rel_avg^2, na.rm = T)/.N,
         tons_month_se_rel = mean(tons_month_se_rel_avg, na.rm = T),
         N = .N),
      by = .(Continent_Region, year, decade)
    ], cols = 'tons_yr')[
      ,':='(tons_yr_se_rel = sd_indiv)
    ]


# Sum pre-dam data by continent, compute SE
Qss_pre_dam_continental <- Qss_pre_byRiver[
  ,':='(Continent_Region = ifelse(Continent_Region %in% c('Europe','Eurasia'),
                                  'Europe/Eurasia', Continent_Region))][
                                    ,.(tons_yr = sum(tons_yr_preDam, na.rm = T)),
                                    by = .(Continent_Region, decade)][
                                      ,':='(year = 1970,
                                            decade = 1970)]

# Make combined table with 
# a) a pre-dam estimate (either pre-dam or Milliman estimate if no pre-dam estimate)
# b) every year in landsat record
Qss_global_annual_ref <- na.omit(rbind(Qss_global_annual, 
                                       Qss_pre_dam_continental, 
                                       use.names = T, fill = T)[
                                         # ,':='(decade = factor(decade, levels = c('Baseline','1990','2000','2010')))
                                         # ,':='(year = factor(year, levels = c('Pre-dam','Baseline',paste0(c(1984:2020)),
                                         #                                        ordered = T)))
                                       ], cols = c('Continent_Region', 'year'))
Qss_global_annual_ref <- Qss_global_annual_ref[Continent_Region != '']
#### 3. EXPORT EXPANDED METADATA TO DRIVE ####
# Write updated tables of rivers to drive

#### 3A. River summaries ####
# Just name, site_no, and id
fwrite(ssc_rivers, file = paste0(wd_imports, 'ssc_outlet_rivers_names.csv'))

# Add discharge and area metadata
fwrite(ssc_rivers_stats, file = paste0(wd_imports, 'ssc_rivers_stats.csv'))

# Full by-river statistics and metadata
fwrite(ssc_rivers_stats_full, file = paste0(wd_imports, 'ssc_rivers_stats_full.csv'))

# Pre-dam estimates from Milliman and Farnsworth, 2012
fwrite(Qss_pre_byRiver, file = paste0(wd_imports, 'Qss_pre_byRiver.csv'))

#### 3B. River Fluxes ####
# Annual by-river flux with metadata
fwrite(Qss_river_annual, file = paste0(wd_imports, 'Qss_river_annual.csv'))

# Annual by-river flux with metadata
fwrite(Qss_river_decadal, file = paste0(wd_imports, 'Qss_river_decadal.csv'))

# Annual by-river SSC with metadata. 
# Only direct observations and no discharge, so weighted towards low SSC
fwrite(SSC_annual_avg, file = paste0(wd_imports, 'SSC_annual_avg_direct_observations_no_discharge.csv'))

#### 3C. Monthly summaries ####
# Monthly averages of discharge, sediment flux
fwrite(avg_monthly_Q_Qss, file = paste0(wd_imports, 'avg_monthly_Q_Qss.csv'))

# Timing of peak suspended sediment concentration
fwrite(Qss_peak_month, file = paste0(wd_imports, 'Qss_peak_month.csv'))

#### 3D. Continent summaries ####
# Drainage area and total discharge, by continent
fwrite(continent_area_Q, file = paste0(wd_imports, 'continent_area_Q.csv'))

# Annual sediment flux, by continent
fwrite(Qss_global_annual_ref, file = paste0(wd_imports, 'Qss_global_annual_ref.csv'))




#### TESTING ####
# Plot Avg. annual SSC, weighted by avg. monthly discharge
ggplot(sed_transport_month_yr[RiverName == 'Mississippi'], aes(x = year + month/12, y = SSC_mgL)) +
  geom_line(alpha = 0.35) +
  season_facet +
  geom_point(data = Qss_river_decadal[RiverName == 'Mississippi'], 
             aes(x = decade + 2.5, y = SSC_mgL), size = 4, color = 'orange') + 
  geom_errorbar(data = Qss_river_decadal[RiverName == 'Mississippi'], 
                aes(x = decade + 2.5, ymin = SSC_mgL - SSC_mgL * sd_total,
                    ymax = SSC_mgL + SSC_mgL * sd_total), color = 'black')

# Monthly by-river plots
ggplot(avg_monthly_Q_Qss[RiverName %chin% c('Changjiang','Huanghe','Mississippi','MacKenzie', 'Amazon', 'Danube','Rhone')], 
       aes(x = factor(month), y = Qss_norm, color = RiverName)) +
  geom_line(aes(group = RiverName)) +
  geom_errorbar(aes(group = RiverName, ymin = Qss_norm - Qss_rel_sd, ymax = Qss_norm + Qss_rel_sd)) +
  facet_wrap(.~Continent_Region)

