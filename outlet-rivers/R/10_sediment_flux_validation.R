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
#### 1A. MISSISSIPPI RIVER ####
# Mississippi River from Tarbert, HISTORICAL
# (1950-2008)
mississippi_historic_insitu <- fread(paste0(wd_imports,'mississippi_tarbert_historical_flux.csv'))

# Import Mississippi data from USGS, MORE RECENT
# (1985-present)
tarbert <- data.table(readNWISqw(siteNumbers = '07295100', parameterCd = '80154'))
tarbert_Q_SSC <- data.table(
  readNWISqw(
    startDate = '1985-01-01', 
    siteNumbers = '07295100', 
    parameterCd = c('80154', '00061')
  ))

#### 1B. YANGTZE RIVER ####
# From Yang et. al, 2018 Annual Yangtze River SSC and Qs
# (1956-2015)
yangtze_insitu <- fread(paste0(wd_imports,'yangzte-1956-2015-sediment.csv'))

#### 1C. ZHUJIANG RIVER ####
# From Yang et. al, 2018 Annual Zhujiang River SSC and Qs
# zhujiang_insitu <- fread(paste0(wd_imports,'zhujiang-1957-2007-sediment.csv'))
# (1954-2016)
zhujiang_insitu <- fread(paste0(wd_imports,'zhujiang-1954-2016-sediment.csv'))

#### 1D. HUAIHE RIVER ####
# From Yang et. al, 2018 Annual Huaihe River SSC and Qs
# (1954-2016)
huaihe_insitu <- fread(paste0(wd_imports,'huaihe-1954-2016-sediment-Li-2018.csv'))

#### 1E. VALIDATION DATA ####
# From satellite estimates
# Monthly
validation_data_monthly <- fread(file = paste0(wd_imports, 'satellite_validation_data_monthly.csv'))
# Annual
# Annual by-river flux with metadata
Qss_river_annual <- fread(file = paste0(wd_imports, 'Qss_river_annual.csv'))

#### 2. INDIVIDUAL RIVER VALIDATION ####
#### 2A. MISSISSIPPI VALIDATION ####
## Historic
# Calculate monthly SSC averages for Mississippi data
tarbert_monthly_table <- tarbert[
  ,':='(
    month = month(sample_dt),
    year = year(sample_dt)
  )][
    ,.(SSC_mgL = mean(result_va, na.rm = T)),
    by = .(month, year)
  ]
# Calculate monthly Q averages for Mississippi data
tarbert_monthly_Q <- tarbert_Q_SSC[parm_cd == '00061'][
  ,':='(
    month = month(sample_dt),
    year = year(sample_dt)
  )][
    ,.(Q_cms = mean(result_va, na.rm = T)*0.02831),
    by = .(month, year)
  ]

# Combine Q and SSC data to compute monthly sediment flux for Mississippi
tarbert_monthly_with_Q <- tarbert_monthly_table[tarbert_monthly_Q, on = c('month', 'year')][
  ,':='(tons_month = 0.0850354 * Q_cms * SSC_mgL * 30.5)
]

# Average annual sediment flux for Mississippi
# (not used)
# tarbert_annual_with_Q <- tarbert_monthly_with_Q[
#   ,.(tons_month = mean(tons_month, na.rm = T)),
#   by = year
# ][,.(tons_yr = mean(tons_month)*12),
#   by = year]

# Average sediment flux by 2-yr periods for Mississippi
tarbert_annual_with_Q_2yr <- tarbert_monthly_with_Q[
  ,.(tons_month = mean(tons_month, na.rm = T),
     tons_month_se = sd(tons_month, na.rm = T)/sqrt(.N),
     N = .N),
  by = year - year%%2
][,.(tons_yr = mean(tons_month, na.rm = T)*12,
     tons_yr_se = mean(tons_month_se)*12,
     N = sum(N, na.rm = T)),
  by = year]

# Combine with more historic data
tarbert_annual_with_Q_2yr <- rbind(
  tarbert_annual_with_Q_2yr[year > 1989 & year < 2022], 
  mississippi_historic_insitu[year < 1990][,.(tons_yr = mean(Mt_yr, na.rm = T)*1e6,
                             tons_yr_se = sd(Mt_yr, na.rm = T)*1e6/sqrt(.N)),
                          by = .(year - year%%2)],
  use.names = T, fill = T)

# From Landsat: Annual Mississippi River SSC and Q
mississippi_annual <- Qss_river_annual[RiverName == 'Mississippi']

# From Landsat: Average Mississippi sediment flux for 2-yr periods
mississippi_annual_2yr <- mississippi_annual[
  ,.(tons_month_avg = mean(tons_month_avg, na.rm = T),
     tons_month_se = sqrt(sum(sd_comb^2 + sd_indiv^2))/sqrt(.N)),
  by = .((year - year%%2))
]

mississippi_error_calc <- data.table(
  satellite_Mt_yr = mississippi_annual_2yr$tons_month_avg*12/1e6,
  insitu_Mt_yr = tarbert_annual_with_Q_2yr[year > 1982][order(year)]$tons_yr/1e6)[
    ,':='(dev = satellite_Mt_yr - insitu_Mt_yr,
          avg = (insitu_Mt_yr + satellite_Mt_yr)/2)
  ]
mississippi_percent_rmse <- round(
  sqrt(sum(mississippi_error_calc$dev^2)/nrow(mississippi_error_calc))/
  mean(mississippi_error_calc$avg)*100,1)

mississippi_percent_cumulative_error <- round(
  (sum(mississippi_error_calc$satellite_Mt_yr) - sum(mississippi_error_calc$insitu_Mt_yr))/
  ((sum(mississippi_error_calc$satellite_Mt_yr) + sum(mississippi_error_calc$insitu_Mt_yr))/2)*100,1)

## Fig. S8A
# Plot sediment flux comparison for satellite and in situ data sources
mississippi_flux_satellite_insitu_comparison_plot <- ggplot() +
  geom_ribbon(data = tarbert_annual_with_Q_2yr,
              aes(x = year + 0.5,
                  ymin = tons_yr/1e6 - tons_yr_se/1e6,
                  ymax = tons_yr/1e6 + tons_yr_se/1e6),
              fill = 'red', alpha = 0.3) +
  geom_ribbon(data = mississippi_annual_2yr,
              aes(x = year + 0.5,
                  ymin = tons_month_avg/1e6*12 - tons_month_avg/1e6*12*tons_month_se,
                  ymax = tons_month_avg/1e6*12 + tons_month_avg/1e6*12*tons_month_se),
              fill = 'navy', alpha = 0.3) +
  geom_line(data = tarbert_annual_with_Q_2yr, 
                            aes(x = year + 0.5, y = tons_yr/1e6), color = 'red') + 
  geom_point(data = tarbert_annual_with_Q_2yr, 
                            aes(x = year + 0.5, y = tons_yr/1e6), color = 'black') + 
  geom_line(data = mississippi_annual_2yr, 
             aes(x = year + 0.5, y = tons_month_avg/1e6*12),
             color = 'navy') +
  geom_point(data = mississippi_annual_2yr, 
             aes(x = year + 0.5, y = tons_month_avg/1e6*12),
             color = 'navy') +
  geom_text_repel(data = tarbert_annual_with_Q_2yr[order(year)][4], aes(
    x = year + 0.5, y = tons_yr/1e6, label = 'Mississippi R. in situ data\nfrom US Geological Survey'),
    nudge_x = 5, nudge_y = 50, hjust = 0) +
  geom_text_repel(data = mississippi_annual_2yr[11], aes(
    x = year + 0.5, y = tons_month_avg/1e6*12, label = 'Satellite data\nfrom this study'),
    nudge_x = -5, nudge_y = 100, hjust = 1) +
  season_facet +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = '',
    y = 'Mt/yr'
  )

tarbert_summary <- tarbert[year(sample_dt) > 1985][
  ,.(SSC_mgL = mean(result_va, na.rm = T),
     N = .N,
     SSC_mgL_se = sd(result_va, na.rm = T)/sqrt(.N)),
  by = year(sample_dt) - year(sample_dt)%%2]

mississippi_outlet_summary <- validation_data_monthly[RiverName == 'Mississippi'][
  ,.(SSC_mgL = mean(SSC_mgL, na.rm = T),
     N = .N,
     SSC_mgL_se = sqrt(sd(SSC_mgL, na.rm = T)^2 + (0.73*mean(SSC_mgL, na.rm = T))^2)/sqrt(.N)),
  by = year - year%%2]


# Avg. annual SSC, weighted by avg. monthly discharge for Mississippi from satellite and in situ
mississippi_ssc_satellite_insitu_comparison_plot <- 
  ggplot(validation_data_monthly[RiverName == 'Mississippi'], aes(x = year + month/12, y = SSC_mgL)) +
  geom_line(alpha = 0.25, color = 'navy') +
  geom_line(data = tarbert_monthly_table, alpha = 0.25, color = 'red') +
  # geom_smooth() +
  season_facet +
  geom_point(data = mississippi_outlet_summary, 
             aes(x = year + 1, y = SSC_mgL), size = 2, color = 'navy') + 
  geom_line(data = mississippi_outlet_summary, 
            aes(x = year + 1, y = SSC_mgL), size = 0.5, color = 'navy') + 
  geom_linerange(data = mississippi_outlet_summary,
                 aes(x = year + 1, ymin = SSC_mgL - SSC_mgL_se,
                     ymax = SSC_mgL + SSC_mgL_se), color = 'navy') +
  geom_point(data = tarbert_summary, aes(x = year + 1, y = SSC_mgL)) + 
  geom_line(data = tarbert_summary, 
            aes(x = year + 1, y = SSC_mgL), size = 0.5, color = 'black') + 
  geom_linerange(data = tarbert_summary, 
                 aes(x = year + 1, ymin = SSC_mgL - SSC_mgL_se,
                     ymax = SSC_mgL + SSC_mgL_se), color = 'red') +
  scale_y_continuous(limits = c(0, 500)) +
  theme(
    legend.position = c(0.8, 0.9)
  ) +
  labs(
    x = '',
    y = 'SSC (mg/L)'
  )


#### 2B. YANGTZE VALIDATION ####
# From Landsat: Annual Mississippi River SSC and Q
yangtze_annual <- Qss_river_annual[RiverName == 'Changjiang']

yangtze_insitu_2yr <- yangtze_insitu[
  ,.(SSC_mgL = mean(SSC_mgL, na.rm = T),
     SSC_mgL_se = sd(SSC_mgL, na.rm = T)/sqrt(.N),
     Mt_yr = mean(Mt_yr, na.rm = T),
     Mt_yr_se = sd(Mt_yr, na.rm = T)/sqrt(.N)),
  by = .(year - year%%2)
]

# From Landsat: Average Yangtze sediment flux for 2-yr periods
yangtze_annual_2yr <- yangtze_annual[
  ,.(tons_month_avg = mean(tons_month_avg, na.rm = T),
     tons_month_se = sqrt(sum(sd_comb^2 + sd_indiv^2))/.N/sqrt(.N)),
  by = .(year - year%%2)
]


yangtze_error_calc <- data.table(
  satellite_Mt_yr = yangtze_annual_2yr[year < 2016]$tons_month_avg*12/1e6,
  insitu_Mt_yr = yangtze_insitu_2yr[year > 1982][order(year)]$Mt_yr)[
    ,':='(dev = satellite_Mt_yr - insitu_Mt_yr,
          avg = (insitu_Mt_yr + satellite_Mt_yr)/2)
  ]
yangtze_percent_rmse <- round(
  sqrt(sum(yangtze_error_calc$dev^2)/nrow(yangtze_error_calc))/
    mean(yangtze_error_calc$avg)*100,1)

yangtze_percent_cumulative_error <- round(
  (abs(sum(yangtze_error_calc$satellite_Mt_yr) - sum(yangtze_error_calc$insitu_Mt_yr)))/
    ((sum(yangtze_error_calc$satellite_Mt_yr) + sum(yangtze_error_calc$insitu_Mt_yr))/2)*100,1)



yangtze_annual_double_mass <- yangtze_annual[
  ,.(tons_yr = cumsum(tons_month_avg * 12),
     Q_cms_yr = cumsum(Q_cms * 12),
     year = year)
]

## Fig. S8B
# Plot sediment flux comparison for satellite and in situ data sources
yangtze_flux_satellite_insitu_comparison_plot <- ggplot() +
  geom_ribbon(data = yangtze_insitu_2yr,
              aes(x = year + 0.5,
                  ymin = Mt_yr - Mt_yr_se,
                  ymax = Mt_yr + Mt_yr_se),
              fill = 'red', alpha = 0.3) +
  geom_ribbon(data = yangtze_annual_2yr,
              aes(x = year + 0.5,
                  ymin = tons_month_avg/1e6*12 - tons_month_avg/1e6*12*tons_month_se,
                  ymax = tons_month_avg/1e6*12 + tons_month_avg/1e6*12*tons_month_se),
              fill = 'navy', alpha = 0.3) +
  geom_line(data = yangtze_insitu_2yr, 
            aes(x = year + 0.5, y = Mt_yr), color = 'red') + 
  geom_point(data = yangtze_insitu_2yr, 
             aes(x = year + 0.5, y = Mt_yr), color = 'black') + 
  geom_line(data = yangtze_annual_2yr, 
            aes(x = year + 0.5, y = tons_month_avg/1e6*12),
            color = 'navy') +
  geom_point(data = yangtze_annual_2yr, 
             aes(x = year + 0.5, y = tons_month_avg/1e6*12),
             color = 'navy') +
  geom_text_repel(data = yangtze_insitu_2yr[5], aes(
    x = year + 0.5, y = Mt_yr, label = 'Yangtze R. in situ data\nfrom Yangtze Water Resources Committee'),
    nudge_x = 5, nudge_y = 50, hjust = 0) +
  geom_text_repel(data = yangtze_annual_2yr[5], aes(
    x = year + 0.5, y = tons_month_avg/1e6*12, label = 'Satellite data\nfrom this study'),
    nudge_x = -5, nudge_y = -50, hjust = 1) +
  season_facet +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = '',
    y = 'Mt/yr'
  )

yangtze_outlet_summary <- validation_data_monthly[RiverName == 'Changjiang'][
  ,.(SSC_mgL = mean(SSC_mgL, na.rm = T),
     N = .N,
     SSC_mgL_se = sqrt(sd(SSC_mgL, na.rm = T)^2 + (0.73*mean(SSC_mgL, na.rm = T))^2)/sqrt(.N)),
  by = year - year%%2]

# Avg. annual SSC, weighted by avg. monthly discharge
yangtze_ssc_satellite_insitu_comparison_plot <- 
  ggplot(validation_data_monthly[RiverName == 'Changjiang'], aes(x = year + month/12, y = SSC_mgL)) +
  # geom_line(alpha = 0.25, color = 'navy') +
  geom_line(data = validation_data_monthly[RiverName == 'Changjiang' & month %in% c(1,2)][
    ,.(SSC_mgL = mean(SSC_mgL, na.rm = T)), by = .(year - year%%2)
  ], aes(x = year), alpha = 0.25, color = 'navy') +
  geom_line(data = validation_data_monthly[RiverName == 'Changjiang' & month %in% c(7,8,9,10)][
    ,.(SSC_mgL = mean(SSC_mgL, na.rm = T)), by = .(year - year%%2)
  ], aes(x = year), alpha = 0.25, color = 'navy') +
  geom_line(data = yangtze_insitu, aes(x = year, y = SSC_mgL_dry_season), alpha = 0.25, color = 'red') +
  geom_line(data = yangtze_insitu, aes(x = year, y = SSC_mgL_wet_season), alpha = 0.25, color = 'red') +
  # geom_smooth() +
  season_facet +
  geom_point(data = yangtze_outlet_summary, 
             aes(x = year + 1, y = SSC_mgL), size = 2, color = 'navy') + 
  geom_line(data = yangtze_outlet_summary, 
            aes(x = year + 1, y = SSC_mgL), size = 0.5, color = 'navy') + 
  geom_linerange(data = yangtze_outlet_summary,
                 aes(x = year + 1, ymin = SSC_mgL - SSC_mgL_se,
                     ymax = SSC_mgL + SSC_mgL_se), color = 'navy') +
  geom_point(data = yangtze_insitu_2yr, aes(x = year + 1, y = SSC_mgL)) + 
  geom_line(data = yangtze_insitu_2yr, 
            aes(x = year + 1, y = SSC_mgL), size = 0.5, color = 'black') + 
  geom_linerange(data = yangtze_insitu_2yr, 
                 aes(x = year + 1, ymin = SSC_mgL - SSC_mgL_se,
                     ymax = SSC_mgL + SSC_mgL_se), color = 'red') +
  scale_y_continuous(limits = c(0, NA)) +
  theme(
    legend.position = c(0.8, 0.9)
  ) +
  labs(
    x = '',
    y = 'SSC (mg/L)'
  )





#### 2C. HUAIHE VALIDATION ####
# From Landsat: Annual Huaihe River SSC and Q
huaihe_annual <- Qss_river_annual[RiverName == 'Huaihe']


huaihe_insitu_2yr <- huaihe_insitu[
  ,.(Mt_yr = mean(Mt_yr, na.rm = T),
     Mt_yr_se = sd(Mt_yr, na.rm = T)/sqrt(.N)),
  by = .(year - year%%2)
]

# From Landsat: Average Yangtze sediment flux for 2-yr periods
huaihe_annual_2yr <- huaihe_annual[
  ,.(tons_month_avg = mean(tons_month_avg, na.rm = T),
     tons_month_se = sqrt(sum(sd_comb^2 + sd_indiv^2))/.N/sqrt(.N)),
  by = .(year - year%%2)
]


huaihe_error_calc <- data.table(
  satellite_Mt_yr = huaihe_annual_2yr[year > 1990 & year < 2018]$tons_month_avg*12/1e6,
  insitu_Mt_yr = huaihe_insitu_2yr[year > 1990 & year < 2018][order(year)]$Mt_yr)[
    ,':='(dev = satellite_Mt_yr - insitu_Mt_yr,
          avg = (insitu_Mt_yr + satellite_Mt_yr)/2)
  ]
huaihe_percent_rmse <- round(
  sqrt(sum(huaihe_error_calc$dev^2)/nrow(huaihe_error_calc))/
    mean(huaihe_error_calc$avg)*100,1)

huaihe_percent_cumulative_error <- round(
  (abs(sum(huaihe_error_calc$satellite_Mt_yr) - sum(huaihe_error_calc$insitu_Mt_yr)))/
    ((sum(huaihe_error_calc$satellite_Mt_yr) + sum(huaihe_error_calc$insitu_Mt_yr))/2)*100,1)



huaihe_annual_double_mass <- huaihe_annual[
  ,.(tons_yr = cumsum(tons_month_avg * 12),
     Q_cms_yr = cumsum(Q_cms * 12),
     year = year)
]

## Fig. S8C
# Plot sediment flux comparison for satellite and in situ data sources
huaihe_flux_satellite_insitu_comparison_plot <- ggplot() +
  geom_ribbon(data = huaihe_insitu_2yr,
              aes(x = year + 0.5,
                  ymin = Mt_yr - Mt_yr_se,
                  ymax = Mt_yr + Mt_yr_se),
              fill = 'red', alpha = 0.3) +
  geom_ribbon(data = huaihe_annual_2yr[year > 1986],
              aes(x = year + 0.5,
                  ymin = tons_month_avg/1e6*12 - tons_month_avg/1e6*12*tons_month_se,
                  ymax = tons_month_avg/1e6*12 + tons_month_avg/1e6*12*tons_month_se),
              fill = 'navy', alpha = 0.3) +
  geom_line(data = huaihe_insitu_2yr,
            aes(x = year + 0.5, y = Mt_yr), color = 'red') +
  geom_point(data = huaihe_insitu_2yr,
             aes(x = year + 0.5, y = Mt_yr), color = 'black') +
  geom_line(data = huaihe_annual_2yr[year > 1986], 
            aes(x = year + 0.5, y = tons_month_avg/1e6*12),
            color = 'navy') +
  geom_point(data = huaihe_annual_2yr[year > 1986], 
             aes(x = year + 0.5, y = tons_month_avg/1e6*12),
             color = 'navy') +
  geom_text_repel(data = huaihe_insitu_2yr[5], aes(
    x = year + 0.5, y = Mt_yr, label = 'Huai R. in situ data\nfrom Li et al, 2018'),
    nudge_x = 7, nudge_y = 8, hjust = 0) +
  geom_text_repel(data = huaihe_annual_2yr[year > 1986][4], aes(
    x = year + 0.5, y = tons_month_avg/1e6*12, label = 'Satellite data\nfrom this study'),
    nudge_x = 5, nudge_y = 6, hjust = 0) +
  season_facet +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = '',
    y = 'Mt/yr'
  ) 

#### 2D. ZHUJIANG VALIDATION ####
# From Landsat: Annual Zhugiang River SSC and Q
zhujiang_annual <- Qss_river_annual[RiverName == 'Zhujiang']

zhujiang_insitu_2yr <- zhujiang_insitu[
  ,.(Mt_yr = mean(Mt_yr, na.rm = T),
     Mt_yr_se = sd(Mt_yr, na.rm = T)/sqrt(.N)),
  by = .(year - year%%2)
]

# From Landsat: Average Yangtze sediment flux for 2-yr periods
zhujiang_annual_2yr <- zhujiang_annual[
  ,.(tons_month_avg = mean(tons_month_avg, na.rm = T),
     tons_month_se = sqrt(sum(sd_comb^2 + sd_indiv^2))/.N/sqrt(.N)),
  by = .(year - year%%2)
]


zhujiang_error_calc <- data.table(
  satellite_Mt_yr = zhujiang_annual_2yr[year > 1990 & year < 2018]$tons_month_avg*12/1e6,
  insitu_Mt_yr = zhujiang_insitu_2yr[year > 1990][order(year)]$Mt_yr)[
    ,':='(dev = satellite_Mt_yr - insitu_Mt_yr,
          avg = (insitu_Mt_yr + satellite_Mt_yr)/2)
  ]
zhujiang_percent_rmse <- round(
  sqrt(sum(zhujiang_error_calc$dev^2)/nrow(zhujiang_error_calc))/
    mean(zhujiang_error_calc$avg)*100,1)

zhujiang_percent_cumulative_error <- round(
  (abs(sum(zhujiang_error_calc$satellite_Mt_yr) - sum(zhujiang_error_calc$insitu_Mt_yr)))/
    ((sum(zhujiang_error_calc$satellite_Mt_yr) + sum(zhujiang_error_calc$insitu_Mt_yr))/2)*100,1)



zhujiang_annual_double_mass <- zhujiang_annual[
  ,.(tons_yr = cumsum(tons_month_avg * 12),
     Q_cms_yr = cumsum(Q_cms * 12),
     year = year)
]

## Fig. S8D
# Plot sediment flux comparison for satellite and in situ data sources
zhujiang_flux_satellite_insitu_comparison_plot <- ggplot() +
  geom_ribbon(data = zhujiang_insitu_2yr,
              aes(x = year + 0.5,
                  ymin = Mt_yr - Mt_yr_se,
                  ymax = Mt_yr + Mt_yr_se),
              fill = 'red', alpha = 0.3) +
  geom_ribbon(data = zhujiang_annual_2yr[year > 1990],
              aes(x = year + 0.5,
                  ymin = tons_month_avg/1e6*12 - tons_month_avg/1e6*12*tons_month_se,
                  ymax = tons_month_avg/1e6*12 + tons_month_avg/1e6*12*tons_month_se),
              fill = 'navy', alpha = 0.3) +
  geom_line(data = zhujiang_insitu_2yr,
            aes(x = year + 0.5, y = Mt_yr), color = 'red') +
  geom_point(data = zhujiang_insitu_2yr,
             aes(x = year + 0.5, y = Mt_yr), color = 'black') +
  geom_line(data = zhujiang_annual_2yr[year > 1990], 
            aes(x = year + 0.5, y = tons_month_avg/1e6*12),
            color = 'navy') +
  geom_point(data = zhujiang_annual_2yr[year > 1990], 
             aes(x = year + 0.5, y = tons_month_avg/1e6*12),
             color = 'navy') +
  geom_text_repel(data = zhujiang_insitu_2yr[5], aes(
    x = year + 0.5, y = Mt_yr, label = 'Zhujiang R. in situ data\nfrom Li et al, 2018'),
    nudge_x = 5, nudge_y = -20, hjust = 0) +
  geom_text_repel(data = zhujiang_annual_2yr[year > 1990][5], aes(
    x = year + 0.5, y = tons_month_avg/1e6*12, label = 'Satellite data\nfrom this study'),
    nudge_x = 5, nudge_y = 20, hjust = 0) +
  season_facet +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = '',
    y = 'Mt/yr'
  )


#### 3. COMBINED VALIDATION PLOTS AND STATISTICS ####
#### 3A. COMBINED VALIDATION PLOTS ####
## Fig. S8
# Combine plots:
# A:Mississippi, B:Yangtze, C:Zhujiang, D:Huaihe
combined_plot_four_rivers <- ((
  (mississippi_flux_satellite_insitu_comparison_plot + 
     theme(plot.title = element_markdown(size = 11),
           axis.title.x = element_blank()) + 
     labs(title = '**Mississippi R.**')) +
    yangtze_flux_satellite_insitu_comparison_plot +
    theme(plot.title = element_markdown(size = 11),
          axis.title.x = element_blank()) + 
    labs(title = '**Yangtze R.**'))) / 
  ((
    zhujiang_flux_satellite_insitu_comparison_plot +
      theme(plot.title = element_markdown(size = 11),
            axis.title.x = element_blank()) + 
      labs(title = '**Zhujiang R.**')) +
     (huaihe_flux_satellite_insitu_comparison_plot + 
        theme(plot.title = element_markdown(size = 11),
              axis.title.x = element_blank()) + 
        labs(title = '**Huai R.**'))) +
  plot_annotation(tag_levels = 'A')

ggsave(combined_plot_four_rivers, filename = paste0(wd_figures, 'FigS8_combined_four_river_validation_plots.png'),
       width = 9.5, height = 9.5)
ggsave(combined_plot_four_rivers, filename = paste0(wd_figures, 'FigS8_combined_four_river_validation_plots.pdf'),
       width = 9.5, height = 9.5, useDingbats = F)

# SSC and SSL combined plot for just Yangtze and Mississippi
# (not used in paper)
combined_plot_ssc_ssl <- (yangtze_ssc_satellite_insitu_comparison_plot + 
                    labs(title = '**Yangtze R**') + 
                    theme(plot.title = element_markdown(size = 9),
                          axis.title.x = element_blank()) +
                    mississippi_ssc_satellite_insitu_comparison_plot + 
                    labs(title = '**Mississippi R.**') + 
                    theme(plot.title = element_markdown(size = 9),
                          axis.title.x = element_blank()))/
  (yangtze_flux_satellite_insitu_comparison_plot +
     theme(axis.title.x = element_blank()) +
     mississippi_flux_satellite_insitu_comparison_plot + 
     theme(axis.title.x = element_blank())) +
  plot_annotation(tag_levels = 'A')

combined_plot_yangtze_mississippi <- (yangtze_flux_satellite_insitu_comparison_plot +
                    theme(plot.title = element_markdown(size = 9),
                          axis.title.x = element_blank()) + 
                    labs(title = '**Yangtze R.**')) /
  (mississippi_flux_satellite_insitu_comparison_plot + 
     theme(plot.title = element_markdown(size = 9),
           axis.title.x = element_blank()) + 
     labs(title = '**Mississippi R.**')) +
  plot_annotation(tag_levels = 'A')

ggsave(combined_plot_yangtze_mississippi, filename = paste0(wd_figures, 'mississippi_yangtze_combined_validation_plots.png'),
       width = 5, height = 9.5)
ggsave(combined_plot_yangtze_mississippi, filename = paste0(wd_figures, 'mississippi_yangtze_combined_validation_plots.pdf'),
       width = 5, height = 9.5, useDingbats = F)








   
#### 3B. COMBINED STATISTICS ####
combined_validation_statistics <- 
  data.table(RiverName = c('Mississippi','Yangtze','Zhujiang','Huaihe'),
             rmse = c(mississippi_percent_rmse, yangtze_percent_rmse,
                                     zhujiang_percent_rmse, huaihe_percent_rmse),
             percent_cumul_error = c(mississippi_percent_cumulative_error, yangtze_percent_cumulative_error,
                                     zhujiang_percent_cumulative_error, huaihe_percent_cumulative_error)
             
  )

# Write to folder
fwrite(combined_validation_statistics, paste0(wd_exports, 'combined_validation_statistics.csv'))


