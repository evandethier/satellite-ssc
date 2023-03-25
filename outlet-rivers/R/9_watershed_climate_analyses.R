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
# Add discharge and area metadata
ssc_rivers_stats <- fread(file = paste0(wd_imports, 'ssc_rivers_stats.csv'))

# Full by-river statistics and metadata
ssc_rivers_stats_full <- fread(file = paste0(wd_imports, 'ssc_rivers_stats_full.csv'))

# Climate data for each watershed
watershed_climate <- fread(paste0(wd_imports, 'outlet-rivers-watershed-climate-stats.csv'))[
  ,':='(
    precip_m_yr = total_precipitation*12,
    site_no = paste0('st_',site_no))
][
  # Join with general watershed statistics
  # Could do annual statistics for more detail
  ssc_rivers_stats[,.(site_no, RiverName, Latitude, Longitude, Q_cms_avg)], on = 'site_no']


#### 2. WATERSHED CLIMATE ANALYSES AND MAPS ####
#### 2A. AVERAGE WATERSHED CLIMATE, RUNOFF RATIO ####
# Map precipitation for 2020
# (not used in paper)
avg_precip_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = watershed_climate[year == 2020], pch = 21, stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = precip_m_yr,
                 size = Q_cms_avg/100)) +
  scale_fill_gradientn(guide = guide_colorbar(order = 1), trans = 'log10', 
                       colors = rev(c('#3671BF','#044D8C','#F2BE5C','#F27405'))) +
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
    fill = '**Precipitation**<br>**(m/yr)**',
    size = '**Avg. Discharge**<br>**(100s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(avg_precip_map, filename = paste0(wd_figures, 'avg_precip_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(avg_precip_map, filename = paste0(wd_figures, 'avg_precip_map.png'),
       width = 11, height = 6)


# Map average watershed runoff ratio
runoff_ratio_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = watershed_climate[year == 2020], pch = 21, stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = (Q_cms_avg * 31.5569)/(precip_m_yr*drainage_area_km2),
                 size = Q_cms_avg/100)) +
  scale_fill_gradientn(guide = guide_colorbar(order = 1), 
                       trans = 'log10', limits = c(0.1,1), oob = squish,
                       colors = rev(c('#183152', '#ABC8E2','#E1E6FA','#E6AC76'))) +
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
    fill = '**Runoff Ratio**',
    size = '**Avg. Discharge**<br>**(100s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(runoff_ratio_map, filename = paste0(wd_figures, 'runoff_ratio_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(runoff_ratio_map, filename = paste0(wd_figures, 'runoff_ratio_map.png'),
       width = 11, height = 6)

# Runoff ratio vs. temperature
runoff_ratio_vs_temperature_plot <- ggplot(data = watershed_climate[year == 2020],
                                           aes(x = mean_2m_air_temperature, y = (Q_cms_avg * 31.5569)/(precip_m_yr*drainage_area_km2))) +
  geom_point() + 
  geom_smooth(method = 'lm', formula = y~x, se = F, lty = 'dashed') +
  season_facet +
  # scale_x_log10() +
  scale_y_log10()

#### 2B. CHANGES IN PRECIPITATION, TEMPERATURE ####
# Calculate change in temperature and precipitation since the 1980s
# Temperature change measured in degrees
# Precipitation change measured in percent change 
climate_river_decadal_change <- dcast.data.table(
  watershed_climate[,':='(decade =  ifelse(year < 1990, 1990,
                                           ifelse(year > 2019, 2015, 
                                                  year - year%%5)))],
  RiverName + Latitude + Longitude ~ decade,
  value.var = c('precip_m_yr', 'mean_2m_air_temperature'), fun.aggregate = mean)[
    ,':='(precip_change_2015_1990 = ifelse(!is.na(precip_m_yr_1990), (precip_m_yr_2015 - precip_m_yr_1990)/precip_m_yr_1990,
                                           (precip_m_yr_2015 - precip_m_yr_2000)/precip_m_yr_2000)*100,
          temperature_change_2015_1990 = ifelse(!is.na(mean_2m_air_temperature_1990), (mean_2m_air_temperature_2015 - mean_2m_air_temperature_1990),
                                                (mean_2m_air_temperature_2015 - mean_2m_air_temperature_2000)))
  ][
    # If the difference between epochs is significant, get the sign of the change
    ,':='(precip_sign = ifelse(abs(precip_change_2015_1990) < 5, 'No change', 
                               ifelse(precip_change_2015_1990 > 0, 'Increase','Decrease')),
          temperature_sign = ifelse(abs(temperature_change_2015_1990) < 0.5, 'No change', 
                                    ifelse(temperature_change_2015_1990 > 0, 'Increase','Decrease')))]

## Fig. S7
# Change in precipitation since the 1980s
precip_change_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = climate_river_decadal_change[ssc_rivers_stats_full[,.(RiverName, Q_cms_avg)], on = 'RiverName'], stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = precip_change_2015_1990,
                 size = Q_cms_avg/1000,
                 shape = precip_sign)) +
  scale_fill_gradientn(limits = c(-30, 30), colors = c('red','grey80','blue'), oob = squish,
                       guide = guide_colorbar(order = 1)) +
  scale_size_continuous(range = c(2,10), breaks = c(1, 10, 100),
                        guide = guide_legend(order = 2)) +
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
    fill = '**% Change**<br>**Precipitation**<br>(2020 - 1990)',
    size = '**Avg. Discharge**<br>**(1,000s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())

ggsave(precip_change_map, filename = paste0(wd_figures, 'FigS7_precip_change_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(precip_change_map, filename = paste0(wd_figures, 'FigS7_precip_change_map.png'),
       width = 11, height = 6)

# Change in temperature since the 1980s
# (not used in paper)
temp_change_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = 'grey95', color = 'grey10', lwd = 0.2) +
  geom_point(data = na.omit(climate_river_decadal_change[ssc_rivers_stats_full[,.(RiverName, Q_cms_avg)], on = 'RiverName'],
                            cols = 'temperature_sign'), stroke = 0.2,
             aes(x = Longitude, y = Latitude, 
                 fill = temperature_change_2015_1990,
                 size = Q_cms_avg/100,
                 shape = temperature_sign)) +
  scale_fill_gradientn(limits = c(-3, 3), colors = rev(c('red','grey80','blue')), oob = squish,
                       guide = guide_colorbar(order = 1)) +
  scale_size_continuous(range = c(2,10), breaks = c(100, 500, 1000),
                        guide = guide_legend(order = 2)) +
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
    fill = '**Temperature**<br>**Change (°C)**<br>(2020 - 1990s)',
    size = '**Avg. Discharge**<br>**(100s m<sup>3</sup>/s)**'
  ) +
  theme(legend.position = c(0.02, 0.4),
        legend.justification = 'left',
        legend.key.width = unit(0.4,"cm"),
        legend.title = element_markdown())
ggsave(temp_change_map, filename = paste0(wd_figures, 'temp_change_map.pdf'),
       width = 11, height = 6, useDingbats = F)
ggsave(temp_change_map, filename = paste0(wd_figures, 'temp_change_map.png'),
       width = 11, height = 6)
