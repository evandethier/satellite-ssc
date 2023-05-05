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
wd_imports <- paste0(wd_root,"/outlet_river_controls/imports/")
# Exports folder (save tables and reports here. May remove if replaced by wd_figures)
wd_exports <- paste0(wd_root,"/outlet_river_controls/exports/")

wd_figures <- paste0(wd_root, "/outlet_river_controls/figures/")
wd_visuals <- paste0(wd_root, "/outlet_river_controls/visuals/")

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_imports, wd_exports, wd_figures, wd_visuals)
for(i in 1:length(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}

#### IMPORT ANNUAL AND MONTHLY DATA ####
qss_river_annual <- fread(paste0(wd_imports, 'Qss_river_annual.csv'))

qss_river_monthly <- fread(paste0(wd_imports, 'sediment_transport_month_yr_data.csv'))



large_river_names <- c('Amazon', 'Chao Phraya', 'Colorado', 'Danube', 'Ebro', 'Fraser', 'Ganges',
                  'Irrawaddy', 'Indus', 'Kolyma', 'Lena', 'MacKenzie', 'Magdalena', 'Mekong', 
                  'Mississippi','Niger', 'Nile', 'Orinoco', 'Paraná', 'Po', 'Rhone', 'Song Hong', 
                  'Yukon', 'Zhujiang', 'Salween', 'Congo', 'Rhine')

# 'Red': 'Song Hong'
# 'Ganges-Brahmaputra': 'Ganges'
# 'Parana': 'Paraná'

# Filter to just large rivers
qss_river_annual_large <- qss_river_annual[RiverName %chin% large_river_names]
qss_river_monthly_large <- qss_river_monthly[RiverName %chin% large_river_names]

# Remove NA values
qss_river_annual_large <- qss_river_annual_large[!is.na(tons_month_avg)]
qss_river_monthly_large <- qss_river_monthly_large[!is.na(tons_month_avg)]

# Evaluate missing names from large river list
# All rivers
unique(qss_river_annual$RiverName)
# Rivers included from list
large_river_names_included <- unique(qss_river_annual_large$RiverName)
# Rivers missing from list
large_river_names_missing <- large_river_names[!(large_river_names %in% large_river_names_included)]
# Which rivers are missing
print(large_river_names_missing)
length(large_river_names_missing)


# Write large river timeseries to file
fwrite(qss_river_annual_large, file = paste0(wd_imports, 'Qss_large_river_annual.csv'))
fwrite(qss_river_monthly_large, file = paste0(wd_imports, 'Qss_large_river_monthly.csv'))

#### ANNUAL FLUX ANALYSIS ####

annual_timeseries_plots <- ggplot(qss_river_annual_large, aes(x = year, y = tons_month_avg/1e6)) +
  geom_errorbar(aes(ymin = tons_month_avg/1e6 - tons_month_avg/1e6*tons_month_se_rel,
                    ymax = tons_month_avg/1e6 + tons_month_avg/1e6*tons_month_se_rel)) +
  geom_line() +
  geom_point() + 
  facet_wrap(.~RiverName, scales = 'free_y') +
  season_facet +
  labs(
    x = 'Year',
    y = 'Qss (tons/month)'
  )

ggsave(annual_timeseries_plots, filename = paste0(wd_figures, 'qss_large_river_annual_timeseries.png'),
       width = 10, height = 9)

#### MONTHLY FLUX ANALYSIS ####
monthly_timeseries_plots <- ggplot(qss_river_monthly_large, 
                                  aes(x = year + month/12, y = tons_month/1e6)) +
  geom_errorbar(aes(ymin = tons_month/1e6 - tons_month_se/1e6,
                    ymax = tons_month/1e6 + tons_month_se/1e6)) +
  geom_line() +
  geom_point() + 
  facet_wrap(.~RiverName, scales = 'free_y', ncol = 1) +
  season_facet +
  labs(
    x = 'Year',
    y = 'Qss (tons/month)'
  )


ggsave(monthly_timeseries_plots, filename = paste0(wd_figures, 'qss_large_river_monthly_timeseries.png'),
       width = 10, height = 14)
