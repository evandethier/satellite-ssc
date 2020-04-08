#### LIBRARY IMPORTS ####
library(dataRetrieval)
library(tidyhydat)

library(readr)
library(readxl)

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)
library(gstat)
library(ggspatial)
library(svglite)
library(plotly)

library(data.table)

library(dplyr) # hoping to move away from dplyr
library(tidyverse)
library(tidyquant)
library(tidyr)
library(broom)
library(modelr)

library(scales)
library(kdensity)
library(NbClust)
library(zoo)
library(segmented)
library(lubridate)
library(reshape2)
library(matrixStats)
library(smoother)


library(glmnet)
library(boot)
library(kernelboot)
library(np)

library(automap)
library(sp)
library(USAboundaries)
library(sf)
library(rgeos)
library(raster)
library(rgdal)
library(maptools)
library(PBSmapping)

#### SET DIRECTORIES ####
  # Set root directory
  wd_root <- "~/satellite-ssc"
  
  # Imports folder (store all import files here)
  wd_imports <- paste0(wd_root,"/ssc-imports/")
  # Exports folder (save all figures, tables here)
  wd_exports <- paste0(wd_root,"/ssc-exports/")
  
  wd_figures <- paste0(wd_exports, "ssc-figures/")
  wd_exports_gc <- paste0(wd_exports,"ssc-gc-plots/")
  wd_station_standalone <- paste0(wd_exports, "ssc-station-vs-standalone-models/")
  wd_standalone_models <- paste0(wd_exports, "ssc-standalone-models/")
  wd_standalone_figures <- paste0(wd_standalone_models, "ssc-standalone-figures/")
  wd_autocorrelation <- paste0(wd_exports, "ssc-autocorrelation/")

  # Create folders within root directory to organize outputs if those folders do not exist
  export_folder_paths <- c(wd_exports, wd_figures, wd_exports_gc,wd_station_standalone, 
                           wd_standalone_models, wd_standalone_figures, wd_autocorrelation)
  for(i in 1:length(export_folder_paths)){
    path_sel <- export_folder_paths[i]
    if(!dir.exists(path_sel)){
      dir.create(path_sel)}
  }
  
#### INITIALIZE MAP DATA FOR N.AMERICA ####
setwd(wd_imports)
# get states shapefile for clipping/display
us_states <-  us_boundaries(map_date = NULL, type = c("state"), resolution = c("low"), states = NULL)
# us_states <- us_states[us_states$state_abbr != "AK" & us_states$state_abbr != "HI" & us_states$state_abbr != "PR",]

# convert shapefile to Spatial class
us_states <- as(us_states, 'Spatial')

# plot(us_states)
# plot(us_states_merge)

# set projection for states shapefile
projection <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# import canadian provinces
canada_prov <- read_sf(dsn = "Canada", layer = "Canada")
# plot(canada_prov)
canada_geom <- st_geometry(canada_prov)
# attributes(canada_geom)

# do conversions and projections for canadian provinces to match us states
canada_prov <- as(canada_prov, 'Spatial')

canada_prov <- spTransform(canada_prov, projection)
# coordinates(canada_prov) <- ~long+lat
proj4string(canada_prov) <- proj4string(us_states)

names(canada_prov) <- c('name','frenchName') # rename columns labeling canadian provinces
canada_prov <- canada_prov[,c('name')] # only select column with province name
us_states <- us_states[,c('name')] # only select column with province name
us_states <- rbind(us_states,canada_prov) # combine canadian and us shapefiles

# fortify state shapefile
us_ca <- fortify(us_states)

#### --- ####
#### IMPORT AND CLEAN -- LANDSAT DATA ####
set.seed(1)
    # raw continuous data for regression
    
    # Import landsat spectral data from each site of interest
    ls_raw <- fread('landsat57_v1.dat')
    
    ls_raw_1 <- na.omit(ls_raw[,':='(site_no = ifelse(name == "",as.character(site_no),gsub('usgs|qw|dv',"",name)),
                                     # Rename columns for simplicity
                                     B1 = B1_median,
                                     B2 = B2_median,
                                     B3 = B3_median,
                                     B4 = B4_median,
                                     B5 = B5_median,
                                     B6 = B6_median,
                                     B7 = B7_median,
                                     nd52 = nd_median,
                                     num_pix = B2_count,
                                     sample_dt = ymd(date),
                                     landsat_dt = ymd(date)
    )], cols = c('B1','B2','B3','B4','B5','B7'))[
      B1 > 0 & B2 > 0 & B3 > 0 & B4 > 0 & B5 > 0 & B7 > 0][
        ,':='( 
          # add new columns with band ratios
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
          Latitude = lat,
          Longitude = lon,
          station_nm = site_no,
          sensor = ifelse(grepl('LT',`system:index`),'Landsat 5','Landsat 7'))][ 
            # select only columns of interest
            ,.(station_nm, sensor, site_no, Latitude,Longitude,sample_dt, num_pix, landsat_dt,
               B1,B2,B3,B4,B5,B6,B7,B2.B1,B2.B1,B3.B1,B4.B1,B5.B1,B7.B1,B3.B2,B4.B2,B5.B2,
               B7.B2,B4.B3,B5.B3,B7.B3,B5.B4,B7.B4,B7.B5,nd52,cloud_cover,cloud_qa_count,cloud_qa_3km,snow_ice_qa_count, snow_ice_qa_3km,
               solar_az, solar_zen,sr_atmos_opacity_median,sr_cloud_qa_median
            )]
    
    ## Identify sites with too few Landsat measurements to be reliable
    
    # Calculate number of satellite samples, per site
    n_sat_samples <- ls_raw_1[,.(N_samples = .N), 
                              by = .(site_no, Latitude, Longitude)]
    
    
    stns_too_narrow <- fread('ssc_stns_too_narrow.dat')
    # Select sites with > 100 satellite measurements, remove sites that are too narrow
    site_no_n100 <- n_sat_samples[N_samples >= 100 &
                                  !(site_no %chin% stns_too_narrow$station_nm), 
                              site_no]
    
    # Filter landsat data by sites with > 100 satellite measurements & wide enough river
    ls_raw_1 <- ls_raw_1[site_no %chin% site_no_n100]
    
    n_sat_samples_n100 <- n_sat_samples[site_no %chin% site_no_n100]
    # Plot number of satellite samples per site as a histogram
    n_sat_samples_histogram <- 
      ggplot(n_sat_samples_n100, aes(x = N_samples)) + 
      geom_histogram(color = 'black', lwd = 0.25, binwidth = 100) +
      geom_vline(data = n_sat_samples_n100[,.(N_samples = mean(N_samples))], aes(xintercept = N_samples),
                 color = 'orange', lty = 'dashed') + 
      geom_text(data = n_sat_samples_n100[,.(N_samples = mean(N_samples))], 
                aes(label = paste0('mean = ',round(N_samples), ' samples'), x = N_samples + 40, y = 120),
                hjust = 0, size = 3) +
      season_facet + 
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
      labs(
        x = 'Number of cloud-free satellite samples/site',
        y = 'Number of sites'
      )
    # Save satellite images/site histogram
    ggsave(n_sat_samples_histogram, filename = paste0(wd_figures,'n_sat_samples_histogram.pdf'), width = 4, height = 4, useDingbats = F)
    

#### IMPORT AND CLEAN -- IN SITU DATA ####
    
    # From USGS
    # Sampling type codes
    usgs_method_codes <- fread('usgs_sampling_type_codes_82398.csv')[,':='(Value = as.character(Value), Code = NULL)]
    # Import raw in situ data from USGS, explicitly make site_no character
    usgs_insitu_raw <- fread(colClasses =list(character = 'site_no'), 'qwdv_data_3km.csv')[,':='(
      # Rename and calculate columns
      agency_cd = agency_cd.x,
      sample_dt = date(sample_dt),
      site_no = as.character(site_no),
      # sample_tm = as.character(sample_tm),
      SSC_mgL = ifelse(!is.na(p80154),as.numeric(p80154),as.numeric(p00530)),
      POC_mgL = as.numeric(p00689),
      p63 = ifelse(as.numeric(p63)>100,NA,as.numeric(p63)), # numeric sand fraction
      sample_depth_m = ifelse(!is.na(p00003), as.numeric(p00003) * 0.3048, as.numeric(p00098)),
      width_m = as.numeric(p00004)*0.3048,
      sample_method = as.character(p82398),
      sampler = NA,
      Latitude = dec_lat_va,
      Longitude = dec_long_va,
      begin_date = date(begin_date),
      end_date = date(end_date),
      Slope = NA,
      Intercept = NA,
      Reference = NA
    )][!is.na(SSC_mgL),.(agency_cd,station_nm, site_no, sample_dt, 
                         SSC_mgL, POC_mgL, p63, sample_method, sampler, sample_depth_m, width_m,
                         Latitude, Longitude, drainage_area_km2, alt_m,
                         data_type,begin_date,end_date, Slope, Intercept, Reference)]
    # Match sample method code to description of method
    setDT(usgs_insitu_raw)[usgs_method_codes, c('Method') := .(i.Method), on = .(sample_method == Value)][
      ,':='(sample_method = Method, Method = NULL)]
    
    # From HYDAT
    # Import raw in situ data from HYDAT (Canadian data)
    ca_insitu_raw <- fread('hydat_sed_data_3km.csv', colClasses =list(character = 'site_no'))[,':='(
      site_no = as.character(site_no),
      sample_dt = date(sample_dt),
      agency_cd = 'HYDAT',
      # sample_tm = as.character(sample_tm),
      sample_depth_m = NA, 
      POC_mgL = NA,
      width_m = NA,
      alt_m = NA,
      Slope = NA,
      Intercept = NA,
      Reference = NA
    )][,':='(
      begin_date = min(sample_dt),
      end_date = max(sample_dt)), by = site_no][
        !is.na(SSC_mgL),.(agency_cd,station_nm, site_no, sample_dt, 
                          SSC_mgL, POC_mgL, p63, sample_method, sampler, sample_depth_m, width_m,
                          Latitude, Longitude, drainage_area_km2, alt_m,
                          data_type,begin_date,end_date, Slope, Intercept, Reference)]
    # From HYBAM
    # Import raw in situ data from South America (HYBAM)
    sa_insitu_raw <- fread('hybam_sed_data_3km.csv')[,'site_no' := as.character(site_no)]
    sa_insitu_site_nos <- unique(sa_insitu_raw[!is.na(site_no),site_no])
    # Import site info, add and rename columns
    hybam_sts <- fread('hybam_sts.csv',
                       colClasses =list(character = 'site_no'))[
                         , .(site_no = as.character(site_no), 
                             river_nm = paste0(`River Name`, ' at ', `Fluviometric Station`), 
                             drainage_area_km2 = as.numeric(`Drainage Area (KmÂ² )`)*1000, 
                             Latitude, Longitude)][site_no %chin% sa_insitu_site_nos]
    
    setkey(sa_insitu_raw,site_no)
    setkey(hybam_sts,site_no)
    # Join raw in situ data with site info, add and rename columns
    sa_insitu_raw_1 <- sa_insitu_raw[hybam_sts][,':='(agency_cd = 'HYBAM',
                                                      site_no = as.character(site_no),                                                  
                                                      sample_depth_m = 'Surface', 
                                                      sample_datetime = ymd_hms(Date),
                                                      sample_dt = date(ymd_hms(Date)),
                                                      # sample_tm = as.character(hms(ymd_hms(Date))),
                                                      Q_cms = NA,
                                                      data_type = 'qw',
                                                      POC_mgL = NA,
                                                      p63 = NA,
                                                      width_m = NA,
                                                      alt_m = NA,
                                                      begin_date = min(date(ymd_hms(Date))),
                                                      end_date = max(date(ymd_hms(Date)))
    )][!is.na(SSC_mgL),.(agency_cd,station_nm, site_no, sample_dt, 
                         SSC_mgL, POC_mgL, p63, sample_method, sampler, sample_depth_m, width_m,
                         Latitude, Longitude, drainage_area_km2, alt_m,
                         data_type,begin_date,end_date)]
    # Import depth integration calculations
    hybam_depth_integration <- fread('amazon-depth-integration.csv')
    setkey(hybam_depth_integration,station_nm)
    setkey(sa_insitu_raw_1,station_nm)
    # Adjust surface measurements to depth integrated for deep SA rivers
    # Add Depth integrated, calculated as the sample method and sample depth
    sa_insitu_raw_2 <- merge(sa_insitu_raw_1,hybam_depth_integration, 
                             all = T, by = 'station_nm')[!is.na(Slope),
                                                         ':='(SSC_mgL = SSC_mgL * Slope + Intercept,
                                                              sample_method = 'Depth integrated, calculated',
                                                              sample_depth_m = 'Depth integrated, calculated')]
    
    # Import raw in situ data from WRA (Taiwan)
    # Import topographic data from SRTM for approximating Taiwan site drainage area
    all_sts_topo <- fread('ls_training_station_topo.dat')[
      , .(agency_cd = agency_cd, 
          site_no = gsub('usgs|qw|dv',"",name),
          name = name, Latitude = lat, Longitude = lon, drainage_area_km2 = drainage_area_km2, alt_m = elev_min)
      ]
    # Only use stations with drainage area > 800 km2 (for Taiwan; threshold for other regions is 3000 km2)
    wra_sts <- all_sts_topo[drainage_area_km2 > 800 & agency_cd == 'WRA']
    wra_site_names <- wra_sts[,name]
    wra_insitu_raw <- fread('taiwan-1999-2014-raw-stationInfo.dat')[SiteEngName %in% wra_site_names]
    
    setkey(wra_sts,name)
    setkey(wra_insitu_raw,SiteEngName)
    
    # Join raw in situ data with site info, add and rename columns
    wra_insitu_raw_2 <- wra_insitu_raw[wra_sts][,':='(
      SSC_mgL = ssc_mgL,
      sample_depth_m = NA,
      station_nm = SiteEngName,
      sample_method = 'Depth Integrated',
      sampler = 'USGS sampler',
      sample_datetime = mdy(datetime_final), # column 'datetime' has 1999 samples with date 2099
      sample_dt = date(mdy(datetime_final)),
      # sample_tm = as.character(hms(ymd_hms(Date))),
      Q_cms = NA,
      data_type = 'qw',
      POC_mgL = NA,
      p63 = NA,
      width_m = NA,
      Slope = NA,
      Intercept = NA,
      Reference = NA
    )][,':='(
      begin_date = min(sample_dt),
      end_date = max(sample_dt)), by = site_no][
        !is.na(SSC_mgL),.(agency_cd,station_nm, site_no, sample_dt, 
                          SSC_mgL, POC_mgL, p63, sample_method, sampler, sample_depth_m, width_m,
                          Latitude, Longitude, drainage_area_km2, alt_m,
                          data_type,begin_date,end_date, Slope, Intercept, Reference)]

#### --- ####
#### FINAL DATA PREPARATION: REMOVE SITES WITH INSUFFICIENT/UNSUITABLE DATA ####
## Bind together in situ data from different agencies
all_acy_insitu_raw <- rbind(usgs_insitu_raw, ca_insitu_raw, sa_insitu_raw_2, wra_insitu_raw_2)[
  is.na(sample_method), sample_method  := 'Unknown'][,
                                                     match_name := ifelse(agency_cd %chin% c('HYDAT','HYBAM'),
                                                                          station_nm,site_no)]

# Select usable (>> drainage area threshold, no canyons) sites with sufficient Landsat images
all_acy_sts <- unique(all_acy_insitu_raw[match_name %chin% site_no_n100,match_name])
site_no_ls_insitu_n100 <- 
  # Remove in situ data at stations with insufficient landsat images
  all_acy_insitu_raw[match_name %chin% all_acy_sts]

## Remove Landsat data from unsuitable sites (insufficient drainage, canyons)
ls_clean <- ls_raw_1[site_no %chin% all_acy_sts]

# Calculate daily in situ mean value of continuous variables (SSC_mgL, , POC_mgL, P63, sample_depth_m, width_m)
# Do this by setting every value of those variables to be the mean daily value at that site and date
# Then remove duplicates by site and date

all_acy_insitu_daily_mean <- all_acy_insitu_raw[, ':='(
                                            SSC_mgL = mean(SSC_mgL, na.rm = T), 
                                            POC_mgL = mean(POC_mgL, na.rm = T), 
                                            p63 = mean(p63, na.rm = T), 
                                            sample_depth_m = mean(as.numeric(sample_depth_m), na.rm = T),
                                            width_m = mean(width_m, na.rm = T)), 
                                                by = c('sample_dt', 'site_no','data_type', 'sample_method','sampler')]

duplicate_insitu_rows <- which(duplicated(all_acy_insitu_daily_mean[,.(site_no, sample_dt, data_type, sample_method, sampler)]))
all_acy_insitu_daily_mean <- all_acy_insitu_daily_mean[-duplicate_insitu_rows]

# Plot sampling method types
sample_method_bar_chart <- ggplot(all_acy_insitu_raw[,.(N = .N),by = sample_method], 
                                  aes(x = reorder(sample_method,N), y = N)) + 
  geom_bar(stat = 'identity') + 
  season_facet +
  scale_y_log10(labels = fancy_scientific) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    x = 'Sampling Method',
    y = 'N in situ samples'
  )


#### JOIN LANDSAT AND IN SITU DATA, WITH LAG OF UP TO 10 DAYS, RESTRICT TO < 3 DAYS ####

## ls_clean: data.table, all Landsat data at sites with > 100 images
## all_acy_insitu_raw: data.table, all in situ data from USGS, HYDAT, HYBAM, WRA

## Join Landsat data with in situ data, allowing for as much as a 10-day lead/lag

# Prepare join
# Set match name
ls_clean[,':='(match_name = site_no)]

# Set match keys to be the match name, sample date
setkeyv(ls_clean,c('match_name'))
setkeyv(all_acy_insitu_daily_mean,c('match_name'))


key(ls_clean)
key(all_acy_insitu_daily_mean)

# Join Landsat data with in situ data
lag_days <- 8
# ls_insitu_raw <- setDT(ls_clean)[all_acy_insitu_daily_mean, roll = lag_days][ # uni-directional join
ls_insitu_raw <- setDT(ls_clean)[, ':='(match_dt_start = sample_dt - lag_days,
                                       match_dt_end = sample_dt + lag_days,
                                       sample_dt = NULL)][
  all_acy_insitu_daily_mean[,':='(match_dt = sample_dt)], 
  on = .(match_name == match_name, match_dt_start <= match_dt, match_dt_end >= match_dt)][ # bi-directional join
  !is.na(B1)][ # removes days with in situ SSC but no Landsat image
    ,lag_days := as.numeric(difftime(sample_dt,landsat_dt),'days')][ # calculate the lead/lag between in situ and satellite sample
      # Remove columns that are duplicated 
      ,':='(site_no = i.site_no, i.site_no = NULL, i.Latitude = NULL, i.Longitude = NULL, i.station_nm = NULL)][
        # Add Log10 SSC, squared cols, and remove no values, too cold (B6 < 269) values
        SSC_mgL > 0,':='(log10_SSC_mgL = log10(SSC_mgL),
                         B1.2 = B1^2,
                         B2.2 = B1^2,
                         B3.2 = B1^2,
                         B4.2 = B1^2,
                         B5.2 = B1^2,
                         B7.2 = B1^2
        )][!is.na(log10_SSC_mgL) & B6 > 268]

# Select minimum lead/lag row
setkey(ls_insitu_raw[,abs_lag_days := abs(lag_days)], abs_lag_days)
ls_insitu_raw <- ls_insitu_raw[, .SD[1], .(site_no, sample_dt, data_type, sample_method, sampler)]

fwrite(ls_insitu_raw,paste0(wd_exports,'ls_insitu_match.csv'))


#### REGRESSION VARIABLES ####
# Select explanatory variables (regressors) for mulitple regression model
regressors_no_site <- c('B1', 'B2', 'B3', 'B4', 'B5', 'B7', # raw bands
                        'B2.B1', 'B3.B1', 'B4.B1', 'B5.B1', 'B7.B1', # band ratios
                        'B3.B2', 'B4.B2', 'B5.B2', 'B7.B2',
                        'B4.B3', 'B5.B3', 'B7.B3',
                        'B5.B4', 'B7.B4', 'B7.B5')
regressors_all <- c('B1', 'B2', 'B3', 'B4', 'B5', 'B7', # raw bands
                    'B1.2', 'B2.2', 'B3.2', 'B4.2', 'B5.2', 'B7.2', # squared bands
                    'site_no', # no clear way to add categorical variables
                    'B2.B1', 'B3.B1', 'B4.B1', 'B5.B1', 'B7.B1', # band ratios
                    'B3.B2', 'B4.B2', 'B5.B2', 'B7.B2',
                    'B4.B3', 'B5.B3', 'B7.B3',
                    'B5.B4', 'B7.B4', 'B7.B5')


# Select regressors - bands and band ratios
regressors_primary <- c(regressors_no_site, 'B4.B3.B1') # all regressors


#### --- ####

#### RUN REGRESSION FOR CLUSTERING WITH 1-7 CLUSTERS -- TAKES ~45 MINS ####
# https://en.wikipedia.org/wiki/Color_quantization something to check out
set.seed(1)
site_band_quantiles_all <- ls_clean[
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

vis_nir_bands <- c('B1','B2','B3','B4','B2.B1','B3.B1','B4.B1','B3.B2','B4.B2','B4.B3','B4.B3.B1')

site_band_scaling_all <- scale(site_band_quantiles_all[,..vis_nir_bands])
cluster_var_combinations <- Map(as.data.frame, sapply(seq_along(vis_nir_bands), function(k) t(combn(vis_nir_bands,k))))
for(i in 4:length(vis_nir_bands)){
  cluster_var_k_sel <- cluster_var_combinations[[i]]
  for(k in 1:nrow(cluster_var_k_sel)){
    print(paste0(i, " ", k))
    cluster_var_sel <- c(as.matrix(cluster_var_k_sel[k,]))
    
  cluster_var_label <- paste(cluster_var_sel, collapse = "_")
  ccc_result <- data.table(cbind(cluster_var_label, i, c(4:10), NbClust(site_band_scaling_all[,cluster_var_sel],
          min.nc=4, max.nc=10, index="ccc", method="kmeans")$All.index))
  
  colnames(ccc_result) <- c('variables','nvars','nclusters','ccc')
  if(k == 1 & i == 4){
    ccc_master <- ccc_result
  }else{
  ccc_master <- rbind(ccc_master, ccc_result)
  }
  if(k%%100 == 0){
    print(ccc_result)
  }
  }}

ccc_analysis <- ccc_master[,':='(nvars = as.numeric(nvars),
                                 nclusters = as.numeric(nclusters),
                                 ccc = as.numeric(ccc))]

ccc_best <- ccc_analysis[nclusters > 5 & nvars < 6][, .(mean_ccc = mean(ccc, na.rm = T)), by = variables][order(-mean_ccc)]

ccc_plot <- ggplot(ccc_analysis, aes(x = factor(nclusters), y = ccc, color = factor(nvars))) + 
  geom_boxplot() +
  # geom_point() + 
  # scale_color_fivethirtyeight() +
  season_facet + 
  theme(legend.position = 'right') + 
  labs(
    x = 'Number of clusters',
    y = 'Cubic clustering criterion',
    color = 'Number of variables'
  )

ggsave(ccc_plot, filename = paste0(wd_exports,'ccc_optimize_plot.pdf'), width = 7, height = 7)

# Calculate k-means cluster based on all regressors at all sites
# # Using raw band and band ratio values
# Select colors for plotting
cl_colors <- brewer.pal(name = 'Paired',n=12)

# Select variables to use for clustering
# clustering_vars <- c('B1','B4','B2.B1','B3.B1', 'B4.B3.B1')
clustering_vars <- unlist(strsplit(as.character(ccc_best[1,'variables']),'_')) # based on optimal cluster vars from ccc analysis
# clustering_vars <- c('B1','B4','B2.B1','B3.B1')

# Compute number of in situ-landsat pairs per station
n_insitu_samples_bySite <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3,.(N_insitu_samples = .N), by = .(site_no, agency_cd)]
# Compute band median at each site for clustering variables
# setkey(n_insitu_samples_bySite,site_no)
setkey(ls_clean,site_no)

write_csv(site_band_quantiles_all, paste0(wd_exports,'all_sts_band_medians.csv'))

## Prepare data for cluster analysis
clusters_calculated_list <- rep(list(NA), 10)
ssc_model_cl_list <- rep(list(NA), 10)
ssc_cluster_color_plot_list <- rep(list(NA), 10)
ssc_cluster_false_color_plot_list <- rep(list(NA), 10)

for(i in c(1:10)){ # test different cluster numbers
# for(i in 5){ # test different cluster numbers
  
  n_centers <- c(1:10)[i]
  
  cluster_col_name <- paste0('cluster_n',n_centers)
  # Calculate k-means cluster based on all regressors at all sites
  # # Using raw band and band ratio values
  site_band_scaling <- scale(site_band_quantiles_all[,..clustering_vars])
  clusters_calculated <- kmeans(site_band_scaling, centers = n_centers,
                                nstart = 20, iter.max = 50)
  
  clusters_calculated_list[[i]] <- clusters_calculated
  # , algorithm = 'MacQueen'
  
  # Compute cluster centers
  cluster_centers <- clusters_calculated$centers
  
  # Assign cluster to each site
  site_band_quantiles_all$cluster <- clusters_calculated$cluster
  
  clustered_sites <- site_band_quantiles_all[,.(site_no,cluster)]
  
  write_csv(site_band_quantiles_all,paste0(wd_exports,'site_band_quantiles_n',i,'.csv'))
  
  
  # TYPICAL RIVER SEDIMENT COLOR AT DIFFERENT CLUSTERS
  
  # Select SSC categories for plotting
  # ssc_categories <- c(0,25,50,100,250,500,750,1000,1500, 1e6)
  ssc_categories <- c(0,50,100,250,500,750,1e6)
  # ssc_categories <- c(0,50,100,200,500,1e6)
  # ssc_categories <- c(0,10,25,50,75,100,150,200,250,300,350, 400, 450, 500,600, 700, 800,900,1000,1100,1500, 1e6)
  
  # Generate SSC labels as 'low value' - 'high value'
  ssc_category_labels <- paste0(ssc_categories[-length(ssc_categories)],'-',c(ssc_categories[-1]))
  # Make highest SSC category "> highest value"
  ssc_category_labels[length(ssc_category_labels)] <- paste0('> ', ssc_categories[length(ssc_category_labels)])
  
  ## Add cluster group as column to ls-insitu matched data.table
  setkey(ls_insitu_raw,match_name)
  setkey(clustered_sites,site_no)
  ls_insitu_cl <- ls_insitu_raw[clustered_sites][
    ,':='(cluster_sel = cluster,
          # # Categorize SSC value as one of selected categories
          ssc_category = cut(10^log10_SSC_mgL, 
                             breaks = ssc_categories,
                             labels = ssc_category_labels))][]
  # Select cluster for analysis
  
  # # Generate median B,G,R, near-infrared for each SSC category and each cluster or site
  ssc_category_color <- ls_insitu_cl[abs(lag_days) < 8, 
                                     keyby = .(cluster_sel, ssc_category),
                                     lapply(.SD, median,na.rm = T),
                                     .SDcols = c('B1','B2','B3', 'B4')]
  
  # Create true-color and false-color plots of 'typical' river color for each SSC category at each cluster group
  for(j in 1:2){
    color_sel <- c('true_color','false_color')[j]
    # raster_color_types <- c(geom_raster(aes(fill = rgb(B3_median/3000,B2_median/3000,B1_median/3000))), # true color
    #                         geom_raster(aes(fill = rgb(B4_median/4000,B3_median/4000,B2_median/4000))) # false color)
    # )
    # data.table version
    raster_color_types <- c(geom_raster(aes(fill = rgb(B3/3000,B2/3000,B1/3000))), # true color
                            geom_raster(aes(fill = rgb(B4/4000,B3/4000,B2/4000))) # false color)
    )
    cluster_ssc_category_color_plot <- 
      ggplot(ssc_category_color, aes(x = cluster_sel, y = ssc_category)) +
      # ggplot(ssc_category_color, aes(x = reorder(paste0(cluster_sel, ' ', site_no), cluster_sel), y = ssc_category)) + # for by site
      raster_color_types[j] +
      scale_fill_identity() +
      season_facet + 
      # scale_x_continuous(expand_scale(add = c(0,0))) + 
      # scale_y_discrete(expand_scale(mult = c(0,0))) +
      theme(axis.text.x = element_text(angle = 90)) + 
      labs(
        y = 'SSC range (mg/L)',
        x = 'River grouping'
      )
    
    ggsave(cluster_ssc_category_color_plot, filename = paste0(wd_figures, cluster_col_name,'_ssc_category_color_plot_',color_sel,'.pdf'),
           width = 3, height = 3.4, useDingbats = F)
    ggsave(cluster_ssc_category_color_plot, filename = paste0(wd_figures, cluster_col_name,'_ssc_category_color_plot_',color_sel,'.png'),
           width = 3, height = 3.4)
    if(j == 1){
      ssc_cluster_color_plot_list[[i]] <- cluster_ssc_category_color_plot
    }else{
      ssc_cluster_false_color_plot_list[[i]] <- cluster_ssc_category_color_plot
    }
  }
  
  # Establish a holdout set for testing statistics
  
  ls_insitu_cl <- getHoldout(ls_insitu_cl)
  # Generate calibration model for each cluster
  ssc_model_cl_iterate <- getModels_lasso(ls_insitu_cl[abs_lag_days < 3 & 
                                                         site_no %chin% n_insitu_samples_bySite[
                                                           N_insitu_samples > 5]$site_no],
                                                      #  .SD[sample(x = .N, 
                                                      #               size = min(.N,min(.N, 500)))], 
                                                      #  by = site_no], 
                                          regressors_all)
  
  ssc_model_cl_list[[i]] <- ssc_model_cl_iterate
  # Predicted values from calibration model
  ssc_model_cl_iterate_pred <- ssc_model_cl_iterate[[1]]
  # Calculate relative error and relative station bias for holdout set using selected cluster model
  # Saves a histogram of all errors
  # Saves a panel boxplot of station bias at each cluster, with each the model type on a different panel
  ssc_model_cl_iterate_rerr_bias <- getErrorBias(ssc_model_cl_iterate_pred, paste0('ssc_model_rerr_', cluster_col_name))
  ssc_model_cl_iterate_rerr <- ssc_model_cl_iterate_rerr_bias[[1]]
  
  # ssc_model_cl_iterate_rmse <- getRMSE(ssc_model_cl_iterate_pred)
  
  # Prepare relative error for plotting
  rel_error_annotate <- data.frame(rel_error = 
                            paste0('Rel. err = ', 
                              round(ssc_model_cl_iterate_rerr[,.(mape_gl_ind,mape_cl_ind,mape_st_ind)],2))[2]) %>% 
    mutate(holdout25 = c('holdout'), SSC_mgL = 1, pred = 39000)
  # Plot actual vs. predicted for holdout. Annotate with RMSE.
  ssc_cluster_iterate_plot_holdout <- get_sscPlot(ssc_model_cl_iterate_pred,"byCluster",'no','no') +
    geom_text(data = rel_error_annotate, 
              aes(x = SSC_mgL, y = pred, label = rel_error), 
              hjust = 0,
              vjust = 0)
  
  # SAVE FIGURE
  # ggsave(ssc_cluster_iterate_plot_holdout, filename = paste0('ssc_', cluster_col_name, '_iterate_plot_holdout.pdf'), useDingbats = F, 
  #        width = 6, height = 7)
  ggsave(ssc_cluster_iterate_plot_holdout, filename = paste0(wd_exports, 'ssc_', cluster_col_name, '_iterate_plot_holdout.png'), 
         width = 6, height = 7)
  
  
  # Calculate model statistics
  n_clusters_df <- data.frame(n_clusters = n_centers)
  
  if(i == 1){
    
    cl_stats <- cbind(data.frame(ssc_model_cl_iterate_rerr),n_clusters_df)
  }else{
    cl_stats <- rbind(cl_stats,
                      cbind(data.frame(ssc_model_cl_iterate_rerr),n_clusters_df))
  }
}

#### CLUSTER STATS -- COMPUTE AFTER INITIAL RUN OF MULTIPLE CLUSTERS, GENERATE FUNCTION FOR ID-ING CLUSTER ####
    colnames(cl_stats) <- c('mape_gl_ind','mape_cl_ind','mape_st_ind',
                            'bias_gl','bias_cl','bias_st',
                            'n_clusters')
    # SAVE TABLE
    # If you don't want to save prior results to Environment
    # write_csv(cl_stats, path = paste0(wd_exports,'ssc_cluster_statistics_nclusters.csv'))
    # cl_stats <- read_csv('ssc_cluster_statistics_nclusters.csv')
    
    cl_stats_melt <- cl_stats %>% melt(measure.vars = c('bias_cl','mape_cl_ind')) %>% 
      mutate(error_term = factor(variable, 
                                 levels = c('bias_cl','mape_cl_ind'),                  
                                 labels = c('Median rel. station bias','Median rel. error'),
                                 ordered = T))
    cluster_optimize_plot <- ggplot(cl_stats_melt, 
                                    aes(x = n_clusters, y = value, group = error_term, color = error_term)) +
      scale_color_manual(values = c('#F26E50','#405173')) +
      geom_point() + 
      geom_line() + 
      scale_y_continuous(limits = c(0.3, 1.1)) + 
      scale_x_continuous(limits = c(1, 10), breaks = c(1:10)) +
      season_facet +
      theme(legend.position = c(0.8,0.8)) +
      labs(
        x = 'Number of cluster groups',
        y = 'Median rel. station bias, Median rel. error',
        color = ''
      )
    
    # SAVE FIGURE
    ggsave(cluster_optimize_plot, filename = paste0(wd_exports,'cluster_optimize_plot.pdf'), width = 4, height = 4, useDingbats = F)

# GENERATE CLUSTERING FUNCTION BASED ON OPTIMAL CLUSTER BREAKDOWN #
    getCluster <- function(df,clustering_vars,n_centers, kmeans_object){
      # Compute band median at each site for clustering variables
      site_band_quantiles_all <- df[
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



#### SELECT MODEL & MAKE ERROR AND UNCERTAINTY CALCULATIONS ####
# Select number of clusters based on best fit
cluster_n_best <- 6
ssc_model_cl_sel <- ssc_model_cl_list[[cluster_n_best]]
# Get SSC prediction for all in situ data
ssc_model_cl_iterate_pred <- ssc_model_cl_sel[[1]]
# ssc_pred_all <- ssc_model_cl_sel[[1]]
# Get prediction function
ssc_cluster_funs <- ssc_model_cl_sel[[2]]

# Export table of function coefficients
for(n in 1:length(ssc_cluster_funs)){
  cv.opt <- coef(ssc_cluster_funs[[n]], s = "lambda.1se")
coef_ex <- cbind(rownames(cv.opt),as.numeric(cv.opt))
colnames(coef_ex) <- c('variable', 'value')

write.table(coef_ex, sep = ",", file = paste0(wd_exports,'cluster_selected_cl', n, '_lasso_fit_coeff.csv'), row.names = F)
}

# Get table of matched in situ-landsat with cluster

# Make cluster parallel plots
# Get median band values per cluster
cluster_band_medians <- ssc_model_cl_iterate_pred[, lapply(.SD, median, na.rm=TRUE), 
                                                  by=.(cluster, ssc_category), .SDcols=paste0('B',c(1:5,7))][
                                                    ssc_category %in% c('0-50','100-250','250-500')
                                                    ]

# cluster_parallel_plot_avg <- ggplot(reshape2::melt(clusters_calculated_list[[cluster_n_best]]$centers, measure.vars = c('B1','B4','B2.B1','B3.B1','B4.B3.B1')), 
#                                     aes(x = Var2, y = value, color = as.factor(Var1), group = Var1)) + 
cluster_parallel_plot_avg <- ggplot(reshape2::melt(cluster_band_medians, measure.vars = paste0('B',c(1:5,7))), 
                                    aes(x = variable, y = value/10000, color = as.factor(cluster), 
                                        group = ssc_category, linetype = paste0(ssc_category, ' mg/L'))) + 
  geom_line(size = 0.25) + 
  geom_point(pch = 21, stroke = 0.25) +
  scale_color_brewer(palette = 'Paired') +
  scale_linetype_manual(values=c("dashed", 'dotdash', 'longdash')) +
  facet_wrap(.~paste0('Cluster ',cluster), ncol = 2) +
  season_facet +
  theme(legend.position = c(0.37, 0.97), legend.background = element_blank()) +
  guides(color = F, linetype = guide_legend(title = element_blank(), 
                                            nrow = 3, keyheight = 0.2, keywidth = 0.7, reverse = T)) +
  labs(
    x = '',
    y = 'Reflectance', 
    color = 'River grouping',
    linetype = 'SSC range (mg/L)'
  )
ggsave(cluster_parallel_plot_avg, filename = paste0(wd_figures,'cluster_parallel_plot_avg.pdf'), width = 4, height = 4, useDingbats = F)

 # Calculate relative error
ssc_model_all_errorbias <- getErrorBias(ssc_model_cl_iterate_pred, 'all_models')

# Save table of match in situ-landsat observations (with cluster)
write_csv(ssc_model_cl_iterate_pred, path = paste0(wd_exports,'ls_insitu_wCluster_sel.csv'))

# Plot SSC actual vs. predicted


# Calculate station bias
setkey(n_insitu_samples_bySite,site_no)
setkey(ssc_model_cl_iterate_pred,site_no)
station_bias <- ssc_model_cl_iterate_pred[
  n_insitu_samples_bySite][
    # N_samples > 12
    ] %>% 
  group_by(cluster_sel, site_no) %>% 
  summarise_at(c('log10_SSC_mgL', 'pred_gl','pred_cl','pred_st'), mean.geometric, na.rm = T) %>%
  mutate(
    bias_gl = (10^median(abs(log10(10^pred_gl/10^log10_SSC_mgL)), na.rm = T)-1),
    bias_cl = (10^median(abs(log10(10^pred_cl/10^log10_SSC_mgL)), na.rm = T)-1),
    bias_st = (10^median(abs(log10(10^pred_st/10^log10_SSC_mgL)), na.rm = T)-1))
 

# Make plots of other variables relative to residuals
# Plot sample depth vs. model error
sample_depth_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred[!is.na(as.numeric(sample_depth_m))], 
                                     aes(x = as.numeric(sample_depth_m), y = pred_cl - log10_SSC_mgL)) + 
  geom_point(alpha = 0.2) + season_facet + 
  geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(as.numeric(sample_depth_m)), .(N_samples = .N)], 
            aes(x = 1, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  scale_y_continuous(limits = c(-4,4)) +
  labs(x = 'Sample depth (m)',
       y = 'Model residual \n(Log10[predicted SSC] - Log10[in situ SSC])')

ggsave(sample_depth_vs_error_plot, filename = paste0(wd_figures,'sample_depth_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

# Sampling methods sorted into categories
all_sample_methods = unique(ssc_model_cl_iterate_pred[,.(sample_method)])
depth_integrated_methods = c('Equal width increment (ewi)','Multiple verticals', 'Equal discharge increment (edi)',
                             'Single vertical', 'SuspSed Box-single ver, depth-int, attached to structure',
                             'Depth Integration','S','Discharge integrated, equal transit rate (etr)',
                             'Point Integration', 'Depth integrated, calculated','Depth Integrated')
surface_methods = c('Grab sample  (dip)', 'Weighted bottle','Submersible pump','L')

single_depth_methods = c('Peristaltic pump','Point sample','SuspSed Pumping - stream sample using a pumping mechanism',
                         'Suction lift peristaltic pump','Thief sample','Composite - Multiple point samples',
                         'SuspSed Partial Depth,depth integrated,part of single vert.','Van Dorn sampler')
                         
miscellaneous_methods = c('Other','Timed sampling interval','E','Bedload, single equal width increment (SEWI)')                         
unknown_methods = c('Unknown')                         
ssc_model_cl_iterate_pred[,':='(depth_category = 
            ifelse(sample_method %in% depth_integrated_methods, 'Depth integrated', 
            ifelse(sample_method %in% surface_methods, 'Surface',
            ifelse(sample_method %in% single_depth_methods, 'Single depth',
            'Unknown or\n Miscellaneous'))))]

sample_method_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred[
              !is.na(depth_category)], 
                                      aes(x = depth_category, y = pred_cl - log10_SSC_mgL)) + 
  geom_boxplot() + 
  season_facet + 
  geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(depth_category), .(N_samples = .N), by = depth_category], 
            aes(x = depth_category, y = 4.1,label = paste0('N = ',N_samples)), hjust = 0.5) +
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(depth_category), .(N_samples = .N)], 
            aes(x = 0.5, y = -4.1,label = paste0('N samples, total = ',N_samples)), hjust = 0) +
  scale_y_continuous(limits = c(-4.2,4.2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = 'Sample method',
       y = 'Model residual \n(Log10[predicted SSC] - Log10[in situ SSC])')
ggsave(sample_method_vs_error_plot, filename = paste0(wd_figures,'sample_method_vs_error_plot.pdf'), 
       width = 5, height = 5, useDingbats = F)

station_bias_density_plot <- ggplot(setDT(station_bias)[!is.na(cluster_sel)][,
                                        cluster_label:=paste0('Cluster ', cluster_sel)], 
                                    aes(x = pred_gl - log10_SSC_mgL)) + 
  geom_vline(aes(xintercept = pred_gl - log10_SSC_mgL), size = 0.1) +
  geom_vline(xintercept = 0, lty = 'dashed',size = 0.75, color = 'orange') +
  geom_density(color = 'red') +
  facet_wrap(.~cluster_label, nrow = 1) + 
  season_facet + 
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(0, 4), expand = expand_scale(add = c(0,0))) + 
  labs(x = 'Avg. model residual, by station')

station_bias_byCluster_density_plot <- ggplot(setDT(station_bias)[!is.na(cluster_sel)][,
                                        cluster_label:=paste0('Cluster ', cluster_sel)], 
                                    aes(x = pred_cl - log10_SSC_mgL)) + 
  geom_vline(aes(xintercept = pred_cl - log10_SSC_mgL), size = 0.1) +
  geom_vline(xintercept = 0, lty = 'dashed',size = 0.75, color = 'orange') +
  geom_density(color = 'red') +
  facet_wrap(.~cluster_label, nrow = 1) + 
  season_facet + 
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(0, 4), expand = expand_scale(add = c(0,0))) + 
  labs(x = 'Avg. model residual, by station')

station_bias_byStation_density_plot <- ggplot(setDT(station_bias)[!is.na(cluster_sel)][,
                                    cluster_label:=paste0('Cluster ', cluster_sel)], 
                                    aes(x = pred_st - log10_SSC_mgL)) + 
  geom_vline(aes(xintercept = pred_st - log10_SSC_mgL), size = 0.1) +
  geom_vline(xintercept = 0, lty = 'dashed',size = 0.75, color = 'orange') +
  geom_density(color = 'red') +
  facet_wrap(.~cluster_label, nrow = 1) + 
  season_facet + 
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(0, 4), expand = expand_scale(add = c(0,0))) + 
  labs(x = 'Avg. model residual, by station')


ggsave(ggarrange(station_bias_density_plot,station_bias_byCluster_density_plot,station_bias_byStation_density_plot,
          nrow = 3, align = 'hv'), filename = paste0(wd_exports, 'station_bias_density_all_plots.pdf'), width = 5, height = 6, useDingbats = F)

p63_vs_error_plot_byCluster <- ggplot(ssc_model_cl_iterate_pred[!is.na(p63)][,
                            cluster_label:=paste0('Cluster ', cluster_sel)], 
                            aes(x = as.numeric(p63)/100, y = pred_cl - log10_SSC_mgL)) + 
  geom_point(alpha = 0.2) + 
  facet_wrap(.~cluster_label) +
  season_facet + 
  geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  geom_text(data = ssc_model_cl_iterate_pred[,
                                             cluster_label:=paste0('Cluster ', cluster_sel)][
                                               !is.na(p63), .(N_samples = .N), by = cluster_label], 
            aes(x = 0, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  scale_y_continuous(limits = c(-4,4)) +
  labs(x = 'Fraction < sand-sized',
       y = 'Model residual \n(Log10[predicted SSC] - Log10[in situ SSC])')

ggsave(p63_vs_error_plot_byCluster, filename = paste0(wd_figures, 'p63_vs_error_byCluster_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

p63_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred[!is.na(p63)], 
                            aes(x = as.numeric(p63)/100, y = pred_cl - log10_SSC_mgL)) + 
  geom_point(alpha = 0.2) + 
  season_facet + 
  geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(p63), .(N_samples = .N)], 
            aes(x = 0, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  scale_y_continuous(limits = c(-4,4)) +
  labs(x = 'Fraction < sand-sized',
       y = 'Model residual \n(Log10[predicted SSC] - Log10[in situ SSC])')

ggsave(p63_vs_error_plot, filename = paste0(wd_figures,'p63_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

POC_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred[!is.na(POC_mgL)], 
                            aes(x = as.numeric(POC_mgL + 1), y = pred_cl - log10_SSC_mgL)) + 
  geom_point(alpha = 0.2) + season_facet + 
  geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(POC_mgL), .(N_samples = .N)], 
            aes(x = 1, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  scale_y_continuous(limits = c(-4,4)) +
  scale_x_log10(labels = fancy_scientific) + 
  labs(x = 'POC (mg/L)',
       y = 'Model residual \n(Log10[predicted SSC] - Log10[in situ SSC])')

ggsave(POC_vs_error_plot, filename = paste0(wd_exports,'POC_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

POC_fraction_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred[!is.na(POC_mgL)], 
                            aes(x = as.numeric(POC_mgL + 1)/SSC_mgL, y = pred_cl - log10_SSC_mgL)) + 
  geom_point(alpha = 0.2) + season_facet + 
  geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(POC_mgL), .(N_samples = .N)], 
            aes(x = 0.001, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  scale_y_continuous(limits = c(-4,4)) +
  scale_x_log10(labels = fancy_scientific, limits = c(0.001, 1)) + 
  labs(x = 'Fraction POC',
       y = 'Model residual \n(Log10[predicted SSC] - Log10[in situ SSC])')

ggsave(POC_fraction_vs_error_plot, filename = paste0(wd_figures,'POC_fraction_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

numpix_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred[!is.na(num_pix)], 
                            aes(x = cut(as.numeric(num_pix), 
                                        breaks = c(0,1,5,10,20,50,100,1e5),
                                        labels = c('1','2-5','6-10','11-20','21-50','50-100','>100')), 
                                        y = abs(pred_cl - log10_SSC_mgL))) + 
  # geom_point(alpha = 0.2) + 
  geom_boxplot(outlier.shape = NA) +
  # geom_smooth(lty = 'dashed', size = 0.75, color = 'orange') + 
  season_facet + 
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(num_pix), .(N_samples = .N), 
                                             by = cut(as.numeric(num_pix), 
                                                      breaks = c(0,1,5,10,20,50,100,1e5),
                                                      labels = c('1','2-5','6-10','11-20','21-50','50-100','>100'))],
            aes(x = cut, y = 1.5,label = paste0('N =\n',N_samples)), size = 3, hjust = 0.5) +
  scale_y_continuous(limits = c(0,1.65)) +
  # scale_x_log10(labels = fancy_scientific) + 
  labs(x = 'Number of pixels sampled',
       y = 'Model error \nAbs(Log10[predicted SSC] - Log10[in situ SSC])')

ggsave(numpix_vs_error_plot, filename = paste0(wd_figures,'numpix_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

drainage_area_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred[!is.na(drainage_area_km2)], 
                                      aes(x = drainage_area_km2, y = pred_cl - log10_SSC_mgL)) + 
  geom_point(alpha = 0.1) + 
  # facet_wrap(.~cluster_sel) +
  season_facet + 
  geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  geom_text(data = ssc_model_cl_iterate_pred[!is.na(drainage_area_km2), .(N_samples = .N)], 
            aes(x = 1, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  scale_y_continuous(limits = c(-4,4)) +
  scale_x_log10(limits = c(500, 1e7), labels = fancy_scientific) + 
  labs(x = 'Drainage area (km2)',
       y = 'Model residual \n(Log10[predicted SSC] - Log10[in situ SSC])')

ggsave(drainage_area_vs_error_plot, filename = paste0(wd_exports, 'drainage_area_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

cloud_cover_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred, 
                                      aes(x = cloud_cover, y = abs(pred_cl - log10_SSC_mgL))) + 
  geom_point(alpha = 0.1) + 
  # facet_wrap(.~cluster_sel) +
  season_facet + 
  # geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  # geom_text(data = ssc_model_cl_iterate_pred[!is.na(drainage_area_km2), .(N_samples = .N)], 
  #           aes(x = 1, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  # scale_y_continuous(limits = c(-4,4)) +
  scale_x_continuous(limits = c(0, 100), expand = expand_scale(mult = c(0, 0.1))) + 
  labs(x = 'Cloud cover (%)',
       y = 'Abs(Model error) \n(|Log10[predicted SSC] - Log10[in situ SSC]|)')

ggsave(cloud_cover_vs_error_plot, filename = paste0(wd_exports,'cloud_cover_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)

cloud_pix_perc_vs_error_plot <- ggplot(ssc_model_cl_iterate_pred, 
                                      aes(x = cloud_qa_count/num_pix, 
                                          y = abs(pred_cl - log10_SSC_mgL))) + 
  geom_point(alpha = 0.1) + 
  # facet_wrap(.~cluster_sel) +
  season_facet + 
  # geom_hline(yintercept = 0, lty = 'dashed', size = 0.75, color = 'orange') +
  # geom_text(data = ssc_model_cl_iterate_pred[!is.na(drainage_area_km2), .(N_samples = .N)], 
  #           aes(x = 1, y = -3.5,label = paste0('N samples = ',N_samples)), hjust = 0) +
  # scale_y_continuous(limits = c(-4,4)) +
  scale_x_continuous(limits = c(0, 100), expand = expand_scale(mult = c(0, 0.1))) + 
  labs(x = 'Local cloud cover ratio',
       y = 'Abs(Model error) \n(|Log10[predicted SSC] - Log10[in situ SSC]|)')

ggsave(cloud_pix_perc_vs_error_plot, filename = paste0(wd_exports,'cloud_pix_perc_vs_error_ncl5.pdf'), 
       width = 5, height = 5, useDingbats = F)










#### --- ####


#### STATION DATA SUMMARY ANALYSES ####
#### STATION DATA SUMMARY -- COMPUTE MEANS, ERROR, BIAS, AT EACH STATION FOR BELOW ANALYSES ####
  # Summarize station data
  ssc_station_summary <- ssc_model_cl_iterate_pred[,.(
                                         log10_SSC_mgL = mean(log10_SSC_mgL, na.rm = T),
                                         N_samples = .N,
                                         POC_mgL = mean(POC_mgL, na.rm = T),
                                         p63 = mean(p63, na.rm = T),
                                         num_pix = mean(num_pix, na.rm = T),
                                         pred_gl = mean(pred_gl, na.rm = T),
                                         pred_cl = mean(pred_cl, na.rm = T),
                                         pred_st = mean(pred_st, na.rm = T),
                                         bias_gl = 10^(median(abs(log10(10^pred_gl/SSC_mgL)))) - 1,
                                         bias_cl = 10^(median(abs(log10(10^pred_cl/SSC_mgL)))) - 1,
                                         bias_st = 10^(median(abs(log10(10^pred_st/SSC_mgL)))) - 1),
                                      keyby = .(agency_cd, site_no, cluster_sel, Latitude, Longitude, 
                                                drainage_area_km2, begin_date, end_date)
                                      ][year(end_date) > 2030, ':='(
                                          begin_date = end_date - years(100), 
                                          end_date = begin_date + years(14))]
   

#### STATION DATA SUMMARY -- CARBON ANALYSIS ####
# Import Carbon data for each station available

# Import carbon data: generated by script: usgs-data-downloads.R
carbon_site_summary <- fread('Carbon_site_summary.csv')[,site_no:=as.character(site_no)][
              nchar(site_no) == 7,site_no := paste0(0,site_no)]

taiwan_carbon <- read_csv('taiwan_poc_hilton2010.dat') 

# # Old dplyr version
# taiwan_carbon_site_summary <- taiwan_carbon %>% mutate(POCFraction = `Organic_carbon_concentration (%)`/100) %>%
#   group_by(site_no) %>% summarise_at('POCFraction', mean, na.rm = T)

# Calculate average POC at each station
taiwan_carbon_site_summary <- setDT(taiwan_carbon)[,POCFraction := `Organic_carbon_concentration (%)`/100][
  , .(POCFraction = mean(POCFraction, na.rm = T),
      POC_N = .N), by = .(River, site_no)
]
carbon_cols <- colnames(taiwan_carbon_site_summary)[which(colnames(taiwan_carbon_site_summary) %in% 
                                                            colnames(carbon_site_summary))]

# Merge carbon site summary from USGS with taiwan data from Hilton et al., 2010
carbon_site_summary_basic <- rbind(carbon_site_summary[,..carbon_cols], taiwan_carbon_site_summary[,..carbon_cols])

# # Summarise in situ SSC data for each site, retain cluster info. (avg. SSC, modify site number to USGS format)
# ssc_station_summary <- ssc_model_cl_iterate_pred %>% group_by(site_no, cluster_sel) %>%
#   summarise_at('log10_SSC_mgL', mean, na.rm = T) %>% ungroup() %>% rowwise() %>%
#   mutate(site_no = getUSGS_site_no(site_no))
# 



# Merge carbon percentage and station cluster/avg. SSC information
# Just USGS data will include all Q data
ssc_station_summary_usgs_carbon <- merge(ssc_station_summary, carbon_site_summary, by = 'site_no')
# Including Taiwan data, which sacrifices Q data
ssc_station_summary_carbon <- merge(ssc_station_summary, carbon_site_summary_basic, by = 'site_no')

# Count number of carbon samples at each cluster group
ssc_station_summary_carbon_count <- ssc_station_summary_carbon %>% group_by(cluster_sel) %>%
  summarise_at('POCFraction',list(count_nu = length))

# Plot average carbon (POC) fraction of SSC (POC (mg/L)/SSC (mg/L)) for each cluster group (boxplot)
carbon_cluster_relation_plot <- ggplot(ssc_station_summary_carbon, 
                                       aes(x = as.factor(cluster_sel), y = POCFraction)) + 
  geom_boxplot(aes(fill = as.factor(cluster_sel))) + 
  geom_text(data = ssc_station_summary_carbon_count, aes(y = 0.2, label = paste0('n = ', count_nu)), size = 2.5) +
  scale_fill_brewer(palette = 'Paired') +
  season_facet + 
  # scale_y_log10() +
  labs(
    x = 'River grouping',
    y = 'Avg. POC Fraction'
  )

# Save carbon by cluster boxplot
ggsave(carbon_cluster_relation_plot, filename = paste0(wd_exports,'carbon_cluster_relation_plot.pdf'),
       useDingbats = F, width = 4, height = 5)

cluster_color_carbon_combined_plot <- ggarrange(ssc_cluster_color_plot_list[[cluster_n_best]],
          carbon_cluster_relation_plot, 
          nrow = 2, align = 'hv', heights = c(2,1.5),
          labels = c('c','d'))

# ggsave(cluster_color_carbon_combined_plot, filename = paste0(wd_figures,'cluster_color_carbon_combined_plot.pdf'), 
#        useDingbats = F, width = 4, height = 5)

cl_model_err_bias <- getErrorBias(ssc_model_cl_iterate_pred, paste0('ssc_model_rerr_', cluster_col_name))

cl_model_err_bias_plot <- ggarrange(cluster_optimize_plot + 
                                      theme(legend.position = c(0.7,0.9), legend.background = element_blank(),
                                            legend.text = element_text(size = 8)),
                                    cl_model_err_bias[[2]], nrow = 2, labels = c('a','b'))

cluster_combined_fig <- ggarrange(cl_model_err_bias_plot, cluster_color_carbon_combined_plot,
                                  ncol = 2)

print(ggplotly(cluster_optimize_plot))

ggsave(cluster_combined_fig, filename = paste0(wd_exports,'cluster_combined_fig.pdf'), width = 7, height = 6, useDingbats = F)
ggsave(cluster_combined_fig, filename = paste0(wd_exports,'cluster_combined_fig.png'), width = 7, height = 6)
ggsave(cluster_combined_fig, filename = paste0(wd_figures,'cluster_combined_fig.jpeg'), width = 7, height = 6)

 ssc_cluster_relation_plot <- ggplot(ssc_station_summary_carbon, 
                                       aes(x = as.factor(cluster_sel), y = 10^log10_SSC_mgL)) + 
  geom_boxplot(aes(fill = as.factor(cluster_sel))) + 
  geom_text(data = ssc_station_summary_carbon_count, aes(y = 0.2, label = paste0('n = ', count_nu))) +
  scale_fill_brewer(palette = 'Paired') +
  season_facet + 
  scale_y_log10() +
  labs(
    x = 'River grouping',
    y = 'Avg. SSC (mg/L)'
  )

# Save SSC by cluster boxplot
ggsave(ssc_cluster_relation_plot, filename = paste0(wd_figures,'ssc_cluster_relation_plot.pdf'), 
       useDingbats = F,width = 4, height = 5)

# Plot POC as a function of SSC, by cluster
ssc_vs_POC_fraction_plot <- ggplot(ssc_station_summary_carbon, 
                              aes(x = 10^log10_SSC_mgL, y = POCFraction, 
                              color = as.factor(cluster_sel))) + 
  geom_point() +
  geom_smooth(method = 'lm', formula = y~exp(-x), color = 'grey20', lty = 'dashed', size = 0.5) +
  scale_color_brewer(palette = 'Paired') +
  facet_wrap(.~paste0('Group ',as.factor(cluster_sel)), nrow = 1) +
  season_facet + 
  scale_x_log10(labels = fancy_scientific) +
  # scale_y_log10() +
  labs(
    color = 'River grouping',
    x = 'Avg. SSC (mg/L)',
    y = 'Avg. POC Fraction'
  )
# Save POC vs SSC plot
ggsave(ssc_vs_POC_fraction_plot, filename = paste0(wd_figures,'ssc_vs_POC_fraction_plot.pdf'), useDingbats = F,
       width = 4, height = 3)

#### STATION DATA SUMMARY -- GRAIN SIZE ####

p63_colors <- brewer.pal(name = 'PuOr',n=11)

grain_size_fractions <- c('0-0.25','0.25-0.5','0.5-0.75','0.75-1')

p63_plot <- ggplot(data = ssc_station_summary[p63 >= 0 & p63 <= 100], 
                   aes(x = Longitude, y = Latitude, 
                       fill = cut(p63, breaks = c(-1,25,50, 75, 101),
                                       labels = grain_size_fractions)
                       )) + 
  geom_map(data = us_ca, map = us_ca, aes(map_id = id, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) +
  geom_point(aes(size = cut(N_samples, breaks = c(0,20, 100, 500, 1e5), 
                            labels = c('1-20', '21-100','101-500', '> 500'))),
                            pch = 21, color = 'black', stroke = 0.2) +
  scale_size_discrete(range = c(1.5,4)) +
  theme_bw() + 
  scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
  # coordinates need to extend to full extent of polygon data or it gets cut off
  scale_y_continuous(limits = c(24.5, 71)) + scale_x_continuous(limits = c(-160, -65)) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.5),
    text = element_text(size=12),
    axis.text = element_text(size = 12)
  ) + 
  labs(
    x = "Longitude",
    y = "Latitude",
    size = 'N samples',
    fill = 'Fraction \n< sand-sized'
  )
ggsave(p63_plot, filename = paste0(wd_figures,'usa_grain_size_map.pdf'), width = 10, height = 5)

#### STATION DATA SUMMARY  -- SAMPLING TIMEFRAME ####
    # Plot begin date to end date of USGS stations
    # Order by begin date
    # Could subset by sampling type with some extra work in the n_sites_param
    sampling_date_range_plot <- ggplot(ssc_station_summary, aes(x = reorder(site_no, begin_date, min, order = T), color = agency_cd)) + 
      geom_linerange(aes(ymin = begin_date, ymax = end_date), lwd = 0.1) +
      scale_color_brewer(palette = 'Paired') + 
      # geom_text(data = n_sites_param, aes(y = begin_date, label = paste0('n = ',`Number of stations`)), hjust = 0) +
      season_facet + 
      # facet_wrap(.~agency_cd) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      labs(
        x = 'Stations (ordered by start date)',
        y = 'Sampling period',
        color = 'Agency'
      ) + rotate()
    
    ggsave(sampling_date_range_plot, filename = paste0(wd_figures,'usgs_sampling_date_range_plot.pdf'), width = 5, height = 5, useDingbats = F)
    
    
    # Get USGS sample dates
    usgs_all_dates_insitu <- fread('param_data_all_wSiteInfo.csv')
    usgs_all_dates_3k_insitu <- usgs_all_dates_insitu[,site_no:= substr(as.character(duplicate_site_check), 1, nchar(as.character(duplicate_site_check))-2)][
      site_no %in% unique(ls_insitu_cl$site_no)
    ]
    usgs_all_dates_3k_insitu <- usgs_all_dates_3k_insitu[,year:=year(sample_dt)]
    setkey(usgs_all_dates_3k_insitu,year)
    usgs_sample_dates <- usgs_all_dates_3k_insitu[, .(agency_cd.x, station_nm, site_no, year)][
      ,':='(agency_cd=agency_cd.x, agency_cd.x = NULL)
    ]
    n_usgs_all_dates_3k_insitu <- usgs_all_dates_3k_insitu[,.N, by = year]
    # Get HYDAT sample dates
    hydat_sample_dates <- ca_insitu_raw[,year:=year(sample_dt)][
                                        ,.(agency_cd, station_nm, site_no, year)]
    
    # Get WRA sample dates
    wra_sample_dates <- wra_insitu_raw_2[,year:=year(sample_dt)][
      ,.(agency_cd, station_nm, site_no, year)][year > 2020, year:=year-100]
    # Get HYBAM/ANA sample dates
    sa_sample_dates <- sa_insitu_raw_2[,year:=year(sample_dt)][
      ,.(agency_cd, station_nm, site_no, year)]
    # Combine all sample data information
    all_agency_sample_dts <- rbind(usgs_sample_dates, hydat_sample_dates, wra_sample_dates, sa_sample_dates)[year < 2019][
      , agency_cd:=gsub('HYDAT','WSC',agency_cd)
    ]
    
    all_agency_samples_per_station_yr <- all_agency_sample_dts[,.N, by = .(year, agency_cd, station_nm, site_no)]
    # Count the number of stations with > X measurements/yr for each agency
    all_agency_samples_per_station_yr_summary <- all_agency_samples_per_station_yr[, 
                    .("N = 1" = sum(N > 0),
                      'N = 10' = sum(N > 10),
                      'N = 50' = sum(N > 50),
                      'N = 300' = sum(N > 300)),
                        by = .(year)]
    
    insitu_sampling_stations_per_yr <- ggplot(
      reshape2::melt(all_agency_samples_per_station_yr_summary, measure.vars = c('N = 1','N = 10','N = 50','N = 300')), 
                                                aes(x = year, y = value, color = variable, group = variable)) + 
      # geom_bar(stat = 'identity', lwd = 0.1, color = 'grey20') + 
      geom_line() +
      # scale_color_manual(values = c('#245473','#F2B366')) +
      scale_color_manual(values = c('#5377A6','#024873','#025959','#027368')) +
      scale_y_continuous(expand = expand_scale(mult = c(0,0.1))) +
      scale_x_continuous(expand = expand_scale(mult = c(0,0))) +
      season_facet + 
      theme(legend.position = c(0.25,0.7)) +
      labs(
        x = '',
        y = 'N stations',
        color = expression(paste('N ', italic("in situ"), " measurements/station"))
      )
    ggsave(insitu_sampling_stations_per_yr, filename = paste0(wd_exports,'insitu_sampling_stations_per_yr.pdf'), width = 5, height = 4, useDingbats = F)
    
    
    all_agency_samples_per_agency_yr <- all_agency_sample_dts[,.N, by = .(year, agency_cd)]
    
    insitu_sampling_all_dates_barplot <- ggplot(all_agency_samples_per_agency_yr, 
                                                aes(x = year, y = N, fill = factor(agency_cd, 
                                                                                   levels = rev(c('USGS','WSC','HYBAM','WRA'))))) + 
      geom_bar(stat = 'identity', lwd = 0.1, color = 'grey20') + 
      scale_fill_manual(values = c('#D96725','#F2B366','#93BBBF','#245473')) +
      scale_y_continuous(expand = expand_scale(mult = c(0,0.1))) +
      scale_x_continuous(expand = expand_scale(mult = c(0,0))) +
      season_facet + 
      theme(legend.position = c(0.15,0.7)) +
      guides(fill = guide_legend(keyheight = 1, keywidth = 1)) +
      labs(
        x = '',
        y = expression(paste('N ', italic("in situ"), " measurements")),
        fill = 'Agency'
      )
    ggsave(insitu_sampling_all_dates_barplot, filename = paste0(wd_exports,'insitu_sampling_all_dates_barplot.pdf'), width = 4, height = 3.5, useDingbats = F)

    ggsave(ggarrange(insitu_sampling_stations_per_yr, insitu_sampling_all_dates_barplot, nrow = 2), 
           filename = paste0(wd_figures,'insitu_sampling_timeframe_combined_plots.pdf'),
           width = 5, height = 6)
#### RUN REGRESSION WITH VARIOUS SUBSETS, COMPUTE STATS ####
# For each subset: 1) take the sample
# Run calibration model on the sample both 
  # a) with the additional variable and 
  # b) with just the subset but the regular variables
# 2) Calculate and plot relative error and station bias for both models
# 3) Calculate model improvement with additional variable

# Subsets:
# The basline data.table is all the matchup days with cluster calculated at the selected N clusters
# Joined with the table summarizing SSC, p63, POC, 
setkeyv(ls_insitu_cl, c('agency_cd', 'match_name', 'Latitude', 'Longitude', 
        'drainage_area_km2', 'begin_date', 'end_date'))

setkeyv(ssc_station_summary, c('agency_cd', 'site_no', 'Latitude', 'Longitude', 
        'drainage_area_km2', 'begin_date', 'end_date'))

ls_insitu_cl_for_sampling <- na.omit(getHoldout(getCluster(ls_insitu_cl[match_name %chin% ssc_station_summary$site_no][
                                                                ssc_station_summary][
                                                                  ,cluster:=NULL
                                                                ], 
                                                   clustering_vars,cluster_n_best, 
                                                   clusters_calculated_list[[cluster_n_best]])), cols = regressors_all)
# i) Using a maximum of N samples/site to see effect on bias (overweighting with high N sites?)
ls_insitu_cl_site_maxN <- ls_insitu_cl_for_sampling[
  abs_lag_days < 1, .SD[sample(x = .N, 
                           size = min(.N,min(.N, 500)))], 
  by = site_no] # sample max 500 samples at all sites
site_maxN_model <- getModels_lasso(ls_insitu_cl_site_maxN, regressors_all)
site_maxN_model_err <- getErrorBias(site_maxN_model[[1]], 'site_maxN')[[1]][
  ,':='(var = 'site_maxN', model_type = c('sample'))
  ]
# ii) Using p63 to see effect of grain size
ls_insitu_cl_gs <- na.omit(ls_insitu_cl_for_sampling[abs_lag_days < 3], cols = 'p63')
gs_model <- getModels_lasso(ls_insitu_cl_gs, c(regressors_all, 'p63'))
gs_normal <- getModels_lasso(ls_insitu_cl_gs, regressors_all)
gs_error <- getErrorBias(gs_model[[1]], 'gs')
gs_normal_error <- getErrorBias(gs_normal[[1]], 'gs_norm')
gs_model_err <- cbind(rbind(gs_error[[1]], gs_normal_error[[1]])[
  ,':='(var = 'gs')
  ], data.table(model = c('with variable', 'without variable')))
# iii) Using site avg. p63 to see effect of grain size with higher N
ls_insitu_cl_gs_avg <- na.omit(ls_insitu_cl_for_sampling[abs_lag_days < 3], cols = 'i.p63')
gs_avg_model <- getModels_lasso(ls_insitu_cl_gs_avg, c(regressors_all, 'i.p63'))
gs_avg_normal <- getModels_lasso(ls_insitu_cl_gs_avg, regressors_all)
gs_avg_error <- getErrorBias(gs_avg_model[[1]], 'gs')
gs_avg_normal_error <- getErrorBias(gs_avg_normal[[1]], 'gs_norm')
gs_avg_model_err <- cbind(rbind(gs_avg_error[[1]], gs_avg_normal_error[[1]])[
  ,':='(var = 'gs_avg')
  ], data.table(model = c('with variable', 'without variable')))
# iv) Using fraction POC (careful, since POC is correlated with SSC)
ls_insitu_cl_POC <- na.omit(ls_insitu_cl_for_sampling[abs_lag_days < 3][,POC_fraction := POC_mgL/SSC_mgL], cols = 'POC_fraction')
POC_model <- getModels_lasso(ls_insitu_cl_POC, c(regressors_all, 'POC_fraction'))
POC_normal <- getModels_lasso(ls_insitu_cl_POC, regressors_all)
POC_error <- getErrorBias(POC_model[[1]], 'POC_fraction')
POC_normal_error <- getErrorBias(POC_normal[[1]], 'POC_norm')
POC_model_err <- cbind(rbind(POC_error[[1]], POC_normal_error[[1]])[
  ,':='(var = 'POC_fraction')
  ], data.table(model = c('with variable', 'without variable')))

# Make table of error and bias
model_error_bias_table <- cbind(data.table('Model' = c('Base','Cluster','Station', 'Grain size')),
                                data.table('Rel. err' = round(c(as.matrix(ssc_model_all_errorbias[[1]])[1:3],
                                                                as.matrix(gs_avg_error[[1]])[2]), 2),
                                           'Rel. bias' = round(c(as.matrix(ssc_model_all_errorbias[[1]])[4:6],
                                                               as.matrix(gs_avg_error[[1]])[5]), 2)))

# v) Using site avg. fraction POC (careful, since POC is correlated with SSC)
ls_insitu_cl_POC_avg <- na.omit(ls_insitu_cl_for_sampling[abs_lag_days < 3][,POC_fraction := i.POC_mgL/(10^i.log10_SSC_mgL)], cols = 'POC_fraction')
POC_avg_model <- getModels_lasso(ls_insitu_cl_POC_avg, c(regressors_all, 'POC_fraction'))
POC_avg_normal <- getModels_lasso(ls_insitu_cl_POC_avg, regressors_all)
POC_avg_error <- getErrorBias(POC_avg_model[[1]], 'POC_fraction')
POC_normal_avg_error <- getErrorBias(POC_avg_normal[[1]], 'POC_fraction_norm')
POC_avg_model_err <- cbind(rbind(POC_avg_error[[1]], POC_normal_avg_error[[1]])[
  ,':='(var = 'POC_fraction_avg')
  ], data.table(model = c('with variable', 'without variable')))
# vi) Using only obs. with lag of 0 days
ls_insitu_cl_lag0 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 1]
lag0_model <- getModels_lasso(ls_insitu_cl_lag0, regressors_all)
# vii) Using only obs. with lag of < 2 days
ls_insitu_cl_lag1 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 2]
lag1_model <- getModels_lasso(ls_insitu_cl_lag1, regressors_all)
# viii) Using obs. with lag of up to 8 days
ls_insitu_cl_lag2 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 3]
lag2_model <- getModels_lasso(ls_insitu_cl_lag2, regressors_all)

ls_insitu_cl_lag3 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 4]
lag3_model <- getModels_lasso(ls_insitu_cl_lag3, regressors_all)

ls_insitu_cl_lag4 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 5]
lag4_model <- getModels_lasso(ls_insitu_cl_lag4, regressors_all)

ls_insitu_cl_lag5 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 6]
lag5_model <- getModels_lasso(ls_insitu_cl_lag5, regressors_all)

ls_insitu_cl_lag6 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 7]
lag6_model <- getModels_lasso(ls_insitu_cl_lag6, regressors_all)

ls_insitu_cl_lag7 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 8]
lag7_model <- getModels_lasso(ls_insitu_cl_lag7, regressors_all)

ls_insitu_cl_lag8 <- ls_insitu_cl_for_sampling[
  abs_lag_days < 9]
lag8_model <- getModels_lasso(ls_insitu_cl_lag8, regressors_all)

lag0_error <- getErrorBias(lag0_model[[1]], 'lag0')
lag1_error <- getErrorBias(lag1_model[[1]], 'lag1')
lag2_error <- getErrorBias(lag2_model[[1]], 'lag2')
lag3_error <- getErrorBias(lag3_model[[1]], 'lag3')
lag4_error <- getErrorBias(lag4_model[[1]], 'lag4')
lag5_error <- getErrorBias(lag5_model[[1]], 'lag5')
lag6_error <- getErrorBias(lag6_model[[1]], 'lag6')
lag7_error <- getErrorBias(lag7_model[[1]], 'lag7')
lag8_error <- getErrorBias(lag8_model[[1]], 'lag8')
lag_model_err <- cbind(rbind(
                  lag0_error[[1]], 
                  lag1_error[[1]], 
                  lag2_error[[1]], 
                  lag3_error[[1]], 
                  lag4_error[[1]], 
                  lag5_error[[1]], 
                  lag6_error[[1]], 
                  lag7_error[[1]], 
                  lag8_error[[1]])[,
                  ,':='(var = 'Lag (d)')],
                  data.table('Lead/Lag (d)' = c(0:8)))
          
model_error_lag_scaling_plot <- ggplot(reshape2::melt(lag_model_err, measure.vars = c('mape_cl_ind','bias_cl')), 
                                                      aes(x = `Lead/Lag (d)`, y = value, color = factor(variable, labels = c('Rel. err', 'Rel. station bias')))) +
    geom_line() + 
    geom_point() + 
  season_facet + 
  theme(legend.position = 'top') +
  scale_color_manual(values = c('#F26E50','#405173')) +
  labs(
    x = 'Lead/Lag (d)',
    y = 'Relative model error/Relative station bias',
    color = ''
  )

ggsave(model_error_lag_scaling_plot, filename = paste0(wd_figures,'model_error_lag_scaling_plot.pdf'), useDingbats = F, width = 3, height = 4)

#### --- ####
#### GRAND CANYON CASE STUDY ####
#### GRAND CANYON -- IMPORT AND CLEAN DATA ####
# IMPORT DATA
    # Import GCMRC manual sampling data as data.table
    # Import GCMRC manual sampling data
    for(i in c(1:8)){
      import_file_name <- c('gcmrc.tsv',paste0('gcmrc-',c(2:8),'.tsv'))[i]
      if(i == 1){
        gcmrc_manual_raw <- fread(import_file_name, colClasses = c('USGS Station #' = 'character'))
      }else{
        gcmrc_manual_raw <- rbind(gcmrc_manual_raw,fread(import_file_name, colClasses = c('USGS Station #' = 'character')))
      }
    }

    # Subset to only include calcuations certified as 'Load' worthy by the GCMRC
    gcmrc_manual_raw <- gcmrc_manual_raw[`Use in load calculations` == 'YES'][,':='
             (SSC_mgL = `Cross-section silt&clay concentration (mg/L)` + `Cross-section sand concentration (mg/L)`,
             date = date(`start time (MST)`),
             month = as.factor(month(`start time (MST)`)),
             site_no = `USGS Station #`, 
             sampling_type = 'Manual',
             Q_cms = NA,
             time_round = round_date(ymd_hms(`start time (MST)`), '15 mins'))][
      !is.na(SSC_mgL)]
    gcmrc_manual_cols <- colnames(gcmrc_manual_raw)
    
    # Import station info
    gc_stns <- gcmrc_manual_raw[,.(avg.SSC_mgL = mean(SSC_mgL),
                                   avg.p63 = mean(`Cross-section silt&clay concentration (mg/L)`/SSC_mgL)), 
                                keyby = .(`USGS Station #`, `Station name`)][
                                  , ':='(site_no = `USGS Station #`,
                                         GCMRC_Station_nm = `Station name`,
                                         `USGS Station #` = NULL,
                                         `Station name` = NULL)
                                  ]
    
  # Import GCMRC automatic sampling data as data.table
    # Import GCMRC continuous sampling data
    gcmrc_cont_raw <- fread('glen_canyon_gcmrc_rawdata_continuous.dat',colClasses = c('station' = 'character'))[
      # Remove NA rows
      !is.na(`Suspended Sand Concentration(mg/L)`) & !is.na(`Suspended Silt-and-Clay Concentration(mg/L)`)
    ][ # Add/modify columns
      ,':='(SSC_mgL = ssc_mgL, # Change SSC column name
            site_no = paste0('0',station), # Add leading zero to station number
            date = date(mdy(date)),
            time_round = round_date(ymd_hms(`time (MST)`), '15 mins')) # Round datetime to 15 minutes
    ][ # summarize by every 15 minutes, taking the median concentration and flow
      , .(Q_cms = median(Discharge_cfs)*0.02832,
          `Cross-section silt&clay concentration (mg/L)` = median(`Suspended Silt-and-Clay Concentration(mg/L)`),
          `Cross-section sand concentration (mg/L)` = median(`Suspended Sand Concentration(mg/L)`),
          SSC_mgL = median(SSC_mgL)),
      keyby = .(site_no, time_round, date)
    ][, ':='( # Add columns for sample method
      `Sampling method` = NA,
       sampling_type = 'Automated')
    ]
    
    
# COMBINE MANUAL AND AUTOMATED SAMPLING FILES BY COMMON COLUMNS    
    # Select common column names for Manual and Automated samples
    gcmrc_colnames <- colnames(gcmrc_cont_raw)[which(colnames(gcmrc_cont_raw) %in% gcmrc_manual_cols)]
    
    # Combine Manual and Automated samples
    gcmrc_union <- rbind(gcmrc_cont_raw[,..gcmrc_colnames],
                         gcmrc_manual_raw[,..gcmrc_colnames])
    
    grain_size_fractions <- c('0-25','25-50','50-75','75-100')
    # Calculate sand percentage
    gcmrc_union <- gcmrc_union[, ':='(
                               silt_clay_percentage = cut(`Cross-section silt&clay concentration (mg/L)`/SSC_mgL, 
                                                          breaks = c(0,0.25,0.5,0.75,1), labels = grain_size_fractions),
                              sample_dt = date(gcmrc_union$time_round))]
    
# PREPARE DATA COLUMNS, MAKE LANDSAT MATCHUP, CALCULATE SSC ESTIMATE
    # Import Landsat data from GCRMC sites
    
    # Landsat data do have station information
    # They also have latitude and longitude
    gcmrc_landsat_raw <- na.omit(fread('glen-canyon-transect-b71000.dat', colClasses = c('station_no' = 'character'))[,
                                                             ':='(site_no = station_no,
                                                                  # Rename columns for simplicity
                                                                  B1 = B1_median,
                                                                  B2 = B2_median,
                                                                  B3 = B3_median,
                                                                  B4 = B4_median,
                                                                  B5 = B5_median,
                                                                  B7 = B7_median,
                                                                  num_pix = B2_count,
                                                                  sample_dt = mdy(date),
                                                                  landsat_dt = mdy(date)
                                                             )]
                                 , cols = c('B1','B2','B3','B4','B5','B7'))[
                                   B1 > 0 & B2 > 0 & B3 > 0 & B4 > 0 & B5 > 0 & B7 > 0][
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
                                       Latitude = lat,
                                       Longitude = lon,
                                       station_nm = paste0(0,station_no),
                                      site_no = paste0(0,site_no))][ 
                                         # select only columns of interest
                                         ,.(site_no, station_nm, distance, Latitude,Longitude,sample_dt, num_pix, landsat_dt,
                                            B1,B2,B3,B4,B5,B7,B2.B1,B3.B1,B4.B1,B5.B1,B7.B1,B3.B2,B4.B2,B5.B2,
                                            B7.B2,B4.B3,B5.B3,B7.B3,B5.B4,B7.B4,B7.B5, B1.2,B2.2,B3.2,B4.2,B5.2,B7.2
                                         )][site_no != "0"]
    
    # Get Landsat matchup data for same location, date as in situ data
    setkeyv(gcmrc_landsat_raw, c('site_no', 'sample_dt'))
    setkeyv(gcmrc_union, c('site_no', 'sample_dt'))
    
    ls_insitu_gc <- gcmrc_union[gcmrc_landsat_raw]
  
    
    # Add/convert columns to those used in the SSC prediction algorithm
    # Add additional column of B4/B3/B1 for cluster analysis
    ls_insitu_gc[,':='(B4.B3.B1 = B4/B3.B1,
                      log10_SSC_mgL = log10(SSC_mgL))]
    
    # Get cluster for each site based on typical spectral profile
    ls_insitu_gc <- getCluster(ls_insitu_gc, 
                                        clustering_vars,cluster_n_best, 
                                        clusters_calculated_list[[cluster_n_best]])
    # Run SSC prediction algorithm to get cluster (and global?) prediction for SSC
    # Calculate SSC prediction for all GCMRC data
    gcmrc_all_pred <- getSSC_pred(na.omit(ls_insitu_gc, cols = c(regressors_all, 'cluster_sel', 'log10_SSC_mgL')), 
                                        regressors_all, ssc_cluster_funs)[,':='(
                            month = month(time_round))]
    
    # ssc_types <- c('Cross.section.silt.clay.concentration..mg.L.', 'Cross.section.sand.concentration..mg.L.','SSC_mgL')
    ssc_types <- c('Cross-section silt&clay concentration (mg/L)', 'Cross-section sand concentration (mg/L)','SSC_mgL')
    ssc_types_name <- c('silt_clay', 'sand','combined')
    
#### GRAND CANYON -- TEST HOW GRAIN SIZE, TIME LAG AFFECTS PREDICTION PERFORMANCE ####
    # Run sensitivity analysis for time lag between in situ and Landsat sample
    
    # Loop through different time gaps between sampling and landsat image acquisition
    max_hours <- 10
    
    co_pred_rmse_info_cols <- c('Hours','Size fraction predict', 'Size fraction')
    
    co_pred_rmse_stats_cols <- c('RMSE_hold','StD_hold','RMSE_in','StD_in')
    
    co_pred_rmse_master <- data.frame(matrix(nrow = max_hours*length(ssc_types)*(length(grain_size_fractions) + 1), 
                                             ncol = (length(co_pred_rmse_info_cols) + length(co_pred_rmse_stats_cols))))
    co_pred_manual_rmse_master <- data.frame(matrix(nrow = max_hours*length(ssc_types)*(length(grain_size_fractions) + 1), 
                                                    ncol = (length(co_pred_rmse_info_cols) + length(co_pred_rmse_stats_cols))))
    
    colnames(co_pred_rmse_master) <- c(co_pred_rmse_info_cols,co_pred_rmse_stats_cols)
    colnames(co_pred_manual_rmse_master) <- c(co_pred_rmse_info_cols,co_pred_rmse_stats_cols)
    
    for(i in 1:3){
      # Select ssc grain size type (silt/clay, sand, combinded)
      ssc_type <- ssc_types[i]
      # Give name to grain size for saving files
      ssc_type_name <- ssc_types_name[i]
      # Select dataset of in situ samples and Landsat spectral information
      co_test_model <- na.omit(ls_insitu_gc, cols = c(regressors_all, 'cluster_sel'))
      
      # Set response variable (in situ SSC) to be only SSC from size-fraction class set by ssc_type
      co_test_model$log10_SSC_mgL <- log10(co_test_model[,..ssc_type])
      
      # Use prediction model to predict SSC at Colorado R. sites (manual and automated samples)
      co_test_pred <- getSSC_pred(co_test_model, regressors_all, ssc_cluster_funs)
      
      # co_test_pred$silt_clay_percentage <- (co_test_pred$`Cross-section sand concentration (mg/L)`/co_test_pred$SSC_mgL) %>%
      #   cut(breaks = c(0,0.25,0.5,0.75,1), labels = c('0-25','25-50','50-75','75-100'))
      
      
      for(k in 1:max_hours){
        # Set window of hours around 10 AM (Landsat acquisition time)
        landsat_hours_sel <- c((10-k):(10+k))
        
        # Set total hour band as a variable for naming/saving
        landsat_hours_count <- length(landsat_hours_sel)
        
        # Subset all SSC data based on on time relative to 10 AM using
        co_pred_nhours <- getHoldout(co_test_pred[
                    hour(time_round) %in% landsat_hours_sel &
                    !is.na(silt_clay_percentage) & 
                    num_pix > 10])
        
        # Subset all SSC data within time window to only include manual sample data
        co_pred_nhours_manual <- co_pred_nhours[sampling_type == 'Manual']
        
        # Calculate RMSE for all data and manual subset, separated by size fraction
        
        # for(l in 1:(length(grain_size_fractions)+1)){
        #   grain_size_sel <- c('All',grain_size_fractions)[l]
        #   if(l == 1){
        #     co_pred_nhour_rmse <- getErrorBias_simple(co_pred_nhours)
        #     co_pred_nhour_manual_rmse <- getErrorBias_simple(co_pred_nhours_manual)
        #   }else{
        #     co_pred_nhour_rmse <- getErrorBias_simple(co_pred_nhours[silt_clay_percentage == grain_size_sel])
        #     co_pred_nhour_manual_rmse <- getErrorBias_simple(co_pred_nhours_manual[silt_clay_percentage == grain_size_sel])
        #   }
        #   
        #   # Each row has rmse statistics for a given time window, grain-size fraction prediction, grain size subset of data
        #   # There are 10 time windows, 3 grain-size fraction predictions, and 5 grain-size subsets: 10*3*5 = 150 rows
        #   # Each row has seven columns.
        #   # Row index: size fraction prediction (i) * (the number of time windows (10) - 1) * time window select (k) - 1
        #   
        #   rmse_index_row <- l+(i-1)*50+5*(k-1)
        #   print(rmse_index_row)
        #   co_pred_rmse_master[rmse_index_row,co_pred_rmse_info_cols] <- c(landsat_hours_count,ssc_type_name, grain_size_sel)
        #   co_pred_rmse_master[rmse_index_row,co_pred_rmse_stats_cols] <- co_pred_nhour_rmse[,1]
        #   
        #   co_pred_manual_rmse_master[rmse_index_row,co_pred_rmse_info_cols] <- c(landsat_hours_count,ssc_type_name, grain_size_sel)
        #   co_pred_manual_rmse_master[rmse_index_row,co_pred_rmse_stats_cols] <- co_pred_nhour_manual_rmse[,1]
        #   print(co_pred_rmse_master[rmse_index_row,])
        # }
        
        
        # Plot in situ SSC (manual and automated sampling) vs predicted SSC, separated by percent sand
        co_reg_plot <- ggplot(co_test_pred[hour(time_round) %in% landsat_hours_sel &
                                         !is.na(silt_clay_percentage) & num_pix > 10,
                                         .(pred_cl=median(pred_cl, na.rm = T)),
                                            by = .(station_nm, date, silt_clay_percentage,log10_SSC_mgL, sampling_type)], 
                              aes(x = 10^log10_SSC_mgL, y = 10^pred_cl, color = silt_clay_percentage)) + 
          geom_point(pch = 16, alpha = 0.2) +
          geom_abline(slope = 1, intercept = 0) +
          facet_wrap(.~silt_clay_percentage) +
          scale_x_log10(limits = c(1,100000), labels = fancy_scientific) +
          scale_y_log10(limits = c(1,100000), labels = fancy_scientific) +
          scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
          season_facet +
          theme(legend.position = 'right') +
          labs(
            x = 'Actual SSC (mg/L)',
            y = 'Estimated SSC (mg/L)',
            color = 'Percent < sand size'
            # color = 'Sand percentage'
            # color = 'Station'
          )
        
        # Plot in situ SSC (manual sampling only) vs predicted SSC, separated by percent sand
        co_reg_plot_manual <- ggplot(co_test_pred[sampling_type == 'Manual' & 
                                                    hour(time_round) %in% landsat_hours_sel &
                                                    !is.na(silt_clay_percentage) & num_pix > 10,
                                                  .(pred_cl=median(pred_cl, na.rm = T)),
                                                  by = .(station_nm, date, silt_clay_percentage,log10_SSC_mgL, sampling_type)],
                                     aes(x = 10^log10_SSC_mgL, y = 10^pred_cl, color = silt_clay_percentage)) + 
          geom_point(pch = 16, alpha = 0.5) +
          geom_abline(slope = 1, intercept = 0) +
          facet_wrap(.~silt_clay_percentage) +
          scale_x_log10(limits = c(1,100000), labels = fancy_scientific) +
          scale_y_log10(limits = c(1,100000), labels = fancy_scientific) +
          scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
          season_facet +
          theme(legend.position = 'right') +
          labs(
            x = 'Actual SSC (mg/L)',
            y = 'Estimated SSC (mg/L)',
            color = 'Percent < sand size'
            # color = 'Sand percentage'
            # color = 'Station'
          )
        
        # Save plots of in situ SSC vs predicted SSC, and raw dataset
        ggsave(co_reg_plot, filename = paste0(wd_exports_gc,'co_holdout_ssc_pred_', ssc_type_name,'_',landsat_hours_count,'hours', '.pdf'), 
               useDingbats = F, width = 7, height = 5) 
        
        ggsave(co_reg_plot_manual, filename = paste0(wd_exports_gc,'co_holdout_manual_ssc_pred_', ssc_type_name,'_',landsat_hours_count,'hours', '.pdf'), 
               useDingbats = F, width = 7, height = 5)
      }
    }
#### GRAND CANYON -- PLOT BIAS DISTRIBUTION FOR GRAIN SIZE CATEGORIES ####
    # Plot bias with grain size fractions
    # Loop to select manual and automatic sampling types
    for(i in 1:2){
      sampling_type_sel <- c('Automated','Manual')[i]
      
      # Subset all in situ samples to only include manual or automated
      gcmrc_bias_subset <- gcmrc_all_pred %>% 
        subset(sampling_type == sampling_type_sel) %>%
        subset(!is.na(silt_clay_percentage)) %>% 
        mutate(bias_cl = ifelse(pred_cl > log10_SSC_mgL, 
                             (10^abs(log10(10^pred_cl/10^log10_SSC_mgL))-1),
                             -1*(10^abs(log10(10^pred_cl/10^log10_SSC_mgL))-1)))
      
      # To reduce number of automated samples, limit to +/- 1 hr
      if(sampling_type_sel == 'Automated'){
        gcmrc_bias_subset <- gcmrc_bias_subset %>% subset(hour(time_round) %in% c(9,10,11))
      }
      
      # Calculate number of samples in each size fraction group
      gcmrc_bias_subset_count <- gcmrc_bias_subset %>% group_by(silt_clay_percentage) %>%
        summarise_at('pred_cl', list(count_sand_fraction = length))
      
      co_grain_size_bias_plot <- ggplot(gcmrc_bias_subset, 
                                        aes(x = silt_clay_percentage, 
                                            y = bias_cl,
                                            fill = silt_clay_percentage)) + 
        geom_boxplot(outlier.shape = NA) + 
        scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
        season_facet +
        # theme(legend.position = 'right') +
        labs(
          x = 'Percent < sand size',
          y = 'Relative error'
          # fill = 'Percent < sand size'
          # color = 'Station'
        )
      if(sampling_type_sel == 'Manual'){
        co_grain_size_bias_plot <- co_grain_size_bias_plot + 
          geom_text(data = gcmrc_bias_subset_count, aes(y = 6, label = paste0('n = ', count_sand_fraction)), size = 3) + 
          scale_y_continuous(limits = c(-10,10))
      }else{
        co_grain_size_bias_plot <- co_grain_size_bias_plot + 
          geom_text(data = gcmrc_bias_subset_count, aes(y = 3.5, label = paste0('n = ', count_sand_fraction)), size = 3) + 
          scale_y_continuous(limits = c(-4,4))
        co_grain_size_bias_plot_auto <- co_grain_size_bias_plot
      }
      
      ggsave(co_grain_size_bias_plot, filename = paste0(wd_figures,'grain_size_bias_plot_',sampling_type_sel,'.pdf'), 
             useDingbats = F, width = 4, height = 5)
    }
    
    co_pred_rmse_master$sampling_type <- 'All'
    co_pred_manual_rmse_master$sampling_type <- 'Manual'
    
    co_pred_rmse_all <- rbind(co_pred_rmse_master,co_pred_manual_rmse_master)
    
    # Plot RMSE for different subcategories of grain size and hour window
    for(i in 1:2){
      sampling_type_sel <- c('All','Manual')[i]
      sampling_type_sel_2 <- c('Automated', 'Manual')[i]
      
      # Select only rows of RMSE calculation that used the sampling type of interest
      co_pred_rmse_sel <- co_pred_rmse_all %>% subset(sampling_type == sampling_type_sel)
      
      # Plot whole model RMSE as a function of sampling lag from Landsat acquisition
      co_rmse_hour_plot <- ggplot(co_pred_rmse_sel %>% subset(`Size fraction predict` == 'combined' & 
                                                                sampling_type == sampling_type_sel), 
                                  aes(x = as.numeric(Hours), y = RMSE_hold, color = `Size fraction`, group = `Size fraction`)) +
        geom_line() + 
        scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
        # facet_wrap(.~`Size fraction predict`) + 
        season_facet +
        theme(legend.position = 'right') +
        labs(
          x = 'Sample lead/lag (hrs)',
          y = 'Relative RMSE',
          color = 'Fraction < sand size'
          # color = 'Sand percentage'
          # color = 'Station'
        )
      ggsave(co_rmse_hour_plot, filename = paste0(wd_exports_gc,'co_rmse_hour_window_',sampling_type_sel,'.pdf'),
             useDingbats = F, width = 5, height = 5)
      
      # Plot residual vs lag hour
      co_residual_hour_lag_plot <- ggplot(gcmrc_all_pred %>% subset(!is.na(silt_clay_percentage) & 
                                                                      sampling_type == sampling_type_sel_2), 
                                          aes(x = as.factor(round(abs(hour(time_round) - 10),1)), 
                                              y = abs(10^pred_cl - SSC_mgL)/SSC_mgL,
                                              fill = silt_clay_percentage)) + 
        geom_boxplot(outlier.shape = NA) + 
        scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
        facet_wrap(.~silt_clay_percentage) +
        scale_y_continuous(limits = c(0,2)) +
        season_facet +
        theme(legend.position = 'right') +
        labs(
          x = 'Sample lead/lag (hrs)',
          y = 'Relative error',
          fill = 'Fraction < sand size'
          # color = 'Station'
        )
      ggsave(co_residual_hour_lag_plot, filename = paste0(wd_exports_gc,'co_residual_hour_',sampling_type_sel_2,'.pdf'),
             useDingbats = F, width = 5, height = 5)
      
      # Plot residual vs. count
      co_residual_count_plot <- ggplot(gcmrc_all_pred %>% subset(!is.na(silt_clay_percentage) & 
                                                                   sampling_type == sampling_type_sel_2), 
                                       aes(x = cut(num_pix, breaks = c(0,25,50,100,150,200, 1000),
                                                   labels = c('0-25','26-50','51-100','101-150','150-200','>200')), 
                                           y = abs(10^pred_cl - SSC_mgL)/SSC_mgL)) + 
        # fill = silt_clay_percentage
        geom_boxplot(outlier.shape = NA, fill = '#D9583B') + 
        scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
        # facet_wrap(.~silt_clay_percentage) +
        scale_y_continuous(limits = c(0,2)) +
        season_facet +
        theme(legend.position = 'right') +
        labs(
          x = 'Landsat pixel count',
          y = 'Relative error'
          # fill = 'Percent < sand size'
          # color = 'Station'
        )
      ggsave(co_residual_count_plot, filename = paste0(wd_exports_gc,'co_residual_count_',sampling_type_sel_2,'.pdf'),
             useDingbats = F, width = 5, height = 5)
    }
    
#### GRAND CANYON -- CALCULATE BY SITE ALGORITHM ####
    # GCMRC-LS matched data already have cluster, regressors, log10_ssc_mgL calculated
    # Merge GCMRC data with existing in situ-Landsat paired data
    # Get holdout, gcmrc yes/no columns
    ls_insitu_gc <- getHoldout(ls_insitu_gc)[,agency_cd:='GCMRC']
    
    # Update global ls in situ data.table to include columns in gc data
    ls_insitu_cl <- ls_insitu_cl[,sampling_type:=sample_method]
    # Find common columns in gc and global calibration landsat match sets
    ls_gc_cols <- c(colnames(ls_insitu_gc)[which(colnames(ls_insitu_gc) %in% colnames(ls_insitu_cl))])
    
    # Bind together GCMRC and global training data
    ls_insitu_cl_gc <- na.omit(
      rbind(
        ls_insitu_cl[,..ls_gc_cols],
          # [-(site_no %in% gc_stns[,site_no])], # remove sites common to GCMRC and USGS
        ls_insitu_gc[, # filter by time close to Landsat time to reduce row N
                   ':='(time_lag_hours = abs(as.numeric(hour(time_round) - 
                                                          as.numeric(hms('10:00:00'))/3600)))][ # add date column
                          time_lag_hours < 2][, # filter by time lag
                                              ..ls_gc_cols]), # select common columns with global dataset
      cols = c(regressors_all,'log10_SSC_mgL'))[ # remove rows with NA values in regressor columns (na.omit function)
        cluster_sel %in% unique(ls_insitu_gc[,cluster_sel])
        ]
    
    # [agency_cd == 'GCMRC', site_no:=paste0()]
    
    
    # Get SSC prediction models for global, cluster, and site
    # Use all regression data from calibration and holdout datasets as input, 
        # Remove automated samples for training data
        # This allows script to run without throwing a memory error
    
    # ls_insitu_cl_gc_sample <- ls_insitu_cl_gc[, .SD[sample(x = .N, size = ceiling(.N/4))], by = site_no] # sample a quarter of samples at all sites
    ls_insitu_cl_gc_sample <- ls_insitu_cl_gc[, .SD[sample(x = .N, 
                                                size = min(.N,min(.N, 500)))], 
                                              by = site_no] # sample max 200 samples at all sites
    
    # Run regression with combined global and gc training data
    # Result is a prediction model for each cluster, and each site. 
    # We are only interested in the site model, because we get cluster model predictions from the
    # cluster model generated without GCMRC calibration data
    site_levels_sel <- unique(ls_insitu_cl_gc[,.(site_no, cluster_sel)])
    # gc_models <- getModels_lasso(ls_insitu_cl_gc[ # for all of the data, not subset
    gc_models <- getModels_lasso(ls_insitu_cl_gc_sample, # using a subset of the data
                        # sampling_type != 'Automated'],
              # , site_no:=factor(site_no, levels = site_levels_sel)],
                                 regressors_all)
    
    # Select data for calculating SSC prediction for each station
    gc_clusters <- unique(ls_insitu_gc[,cluster_sel])
    
    gc_st_model <- getSSC_pred_st(ls_insitu_cl_gc[
      sampling_type %chin% c('Automated','Manual') # to subset for only sites in GCMRC
      ],
                regressors_all, gc_models[[3]], 
        levels_all = site_levels_sel
      ) # add level spec.: ,levels_all = ls_insitu_cl_gc[,.(site_no, cluster_sel)]
    
    # Merge site prediction with data.table that has cluster prediction
    
    setkeyv(gcmrc_all_pred, c('site_no', 'B1','B2', 'sample_dt','SSC_mgL'))
    setkeyv(gc_st_model, c('site_no', 'B1','B2', 'sample_dt','SSC_mgL'))
    

    gc_cl_st <- gc_st_model[,.(site_no, B1,B2, sample_dt, pred_st, SSC_mgL, sample_dt)][
            gcmrc_all_pred][
              !is.na(pred_st)
            ]
    
    # Add year, month, station name
    setkey(gc_cl_st, site_no)
    setkey(gc_stns, site_no)
    gc_cl_st <- gc_cl_st[,':='(year = year(sample_dt),
                     month = month(sample_dt))][
                       gc_stns 
                     ][
                       GCMRC_Station_nm == "Colorado River above Little Colorado River near Desert View, AZ" & 
                         month %in% c(9,10) & SSC_mgL == 4, log10_SSC_mgL := NA # SSC = 4 appears to be an error value
                       ]
#### GRAND CANYON -- CALCULATE AND PLOT OVERALL MODEL ERROR AND BIAS ####
    # Calculate model error for the Grand Canyon sites
    gc_model_error <- getErrorBias(getHoldout(gc_cl_st)[,':='(drainage_area_km2 = NA, 
                                                  pred_gl=NA, 
                                                  agency_cd = NA,
                                                  begin_date = NA,
                                                  end_date = NA, 
                                                  p63 = NA)], 'grand_canyon')
    
    gc_stns_melt <- melt(gc_cl_st, 
                         id.vars = c('sample_dt','station_nm','GCMRC_Station_nm','month'),
                         measure.vars = c('pred_cl','pred_st','log10_SSC_mgL'))[
      !(GCMRC_Station_nm %in% c('Colorado River at Potash, UT', 'Colorado River at Lees Ferry, AZ'))][
        ,GCMRC_Station_nm := factor(GCMRC_Station_nm)
      ]
    
    levels(gc_stns_melt$variable)[c(1,2,3)] <- c('Cluster model','Station model','in situ meas.')
    levels(gc_stns_melt$GCMRC_Station_nm)[c(1,2,3,4)] <- c('Colorado R.\n60 Mile Station',
                                                           'Colorado R.\nGrand Canyon, AZ',
                                                           'Colorado R.\n30 Mile Station',
                                                           'Green R.\nMineral Bottom, UT')
    
    gc_stns_melt$GCMRC_Station_nm<- ordered(gc_stns_melt$GCMRC_Station_nm, levels = c('Green R.\nMineral Bottom, UT','Colorado R.\n30 Mile Station',
                                                                               'Colorado R.\n60 Mile Station','Colorado R.\nGrand Canyon, AZ'))
    
    co_labels <- data.frame(month = rep(1,4), value = 1, GCMRC_Station_nm = levels(gc_stns_melt$GCMRC_Station_nm)[c(1,2,3,4)])
    co_pred_boxes <- ggplot(gc_stns_melt, 
                            # aes(x = as.factor(month-month%%2), 
                            aes(x = as.factor(month), 
                                y = 10^value)) + geom_boxplot(outlier.shape = NA, aes(fill = variable), lwd = 0.25) + 
      geom_text(data = co_labels, aes(label = GCMRC_Station_nm, x = month, y = value), hjust = 0, vjust = 0, size = 4) +
      scale_y_log10(limits = c(1, 10000), labels = fancy_scientific) +
      scale_fill_manual(values = c('light blue','#1e8ecf','orange')) +
      facet_wrap(.~GCMRC_Station_nm) + theme_evan_facet +
      theme(
        strip.background = element_blank(),
        # strip.text = element_text(hjust = 0,size = 9),
        strip.text = element_blank(),
        axis.text = element_text(color = 'black'),
        legend.position = c(0.91,0.94),
        legend.background = element_blank()) +
      labs(
        x = 'Month',
        y = 'SSC (mg/L)',
        fill = ""
      )
    
    
    ggsave(co_pred_boxes, filename = paste0(wd_figures,'CO_model_compare.pdf'),useDingbats = F, width = 7, height = 7)
    ggsave(co_pred_boxes, filename = paste0(wd_figures,'CO_model_compare.png'), width = 7, height = 7)




#### GRAND CANYON -- ERROR ANALYSIS BY STATION, GRAIN SIZE ####
# Relationship between grain size and other factors

    
  # Calculate grain size at different stations per month
    gcmrc_stns_monthly <- gcmrc_manual_raw[,.(avg.SSC_mgL_monthly = mean(SSC_mgL, na.rm = T),
                                              avg.p63_monthly = mean(`Cross-section silt&clay concentration (mg/L)`/SSC_mgL, na.rm = T)), 
                                              keyby = .(`USGS Station #`, `Station name`, month)][
                                             , ':='(site_no = `USGS Station #`,
                                                    GCMRC_Station_nm = `Station name`,
                                                    `USGS Station #` = NULL,
                                                    `Station name` = NULL,
                                                    month = as.numeric(month))
                                             ]
  # Calculate model error at Grand Canyon stations per month
    gcmrc_stn_month_error <- gc_cl_st[,.(avg.SSC = 10^mean(log10_SSC_mgL, na.rm = T),
                                         avg.cl.model.error = mean((pred_cl - log10_SSC_mgL), na.rm = T),
                                         avg.st.model.error = mean((pred_st - log10_SSC_mgL), na.rm = T),
                                         Latitude = mean(Latitude),
                                         Longitude = mean(Longitude)),
                                      keyby = .(site_no, GCMRC_Station_nm, month)][
                                        # ,month:=as.numeric(month)
                                      ][GCMRC_Station_nm %in% c('Green River at Mineral Bottom nr Cynlnds Ntl Park',
                                                                'Colorado River near river mile 30',
                                                                'Colorado River above Little Colorado River near Desert View, AZ',
                                                                'Colorado River near Grand Canyon, AZ')]
    
    gc_seasonal_SSC_plot <- ggplot(gcmrc_stn_month_error, aes(x = as.factor(month), y = avg.SSC, 
                                      color = factor(GCMRC_Station_nm,labels = 
                                                       c('Green R.\nMineral Bottom, UT','Colorado R.\n30 Mile Station',
                                                         'Colorado R.\n60 Mile Station','Colorado R.\nGrand Canyon, AZ')),
                                      group = GCMRC_Station_nm)) + 
      geom_line() + 
      geom_point(pch = 21) +
      facet_wrap(.~factor(GCMRC_Station_nm,
                          levels = c('Green River at Mineral Bottom nr Cynlnds Ntl Park',
                                     'Colorado River near river mile 30',
                                     'Colorado River above Little Colorado River near Desert View, AZ',
                                     'Colorado River near Grand Canyon, AZ'),
                          labels = 
                            c('Green R.\nMineral Bottom, UT','Colorado R.\n30 Mile Station',
                              'Colorado R.\n60 Mile Station','Colorado R.\nGrand Canyon, AZ')), nrow = 4) +
      season_facet + 
      scale_color_manual(values = c('#D96725','#245473','#93BBBF','#F2B366')) +
      theme(legend.position = 'top') + 
      labs(x = 'Month',
           y = 'Avg. SSC (mg/L)',
           color = 'Station')
  ggsave(gc_seasonal_SSC_plot, filename = paste0(wd_figures,'gc_seasonal_SSC_plot.pdf'), useDingbats = F, width = 5, height = 5)  
  # Join grain size and model error calculations    
    
    setkeyv(gcmrc_stns_monthly,c('site_no','month'))
    setkeyv(gcmrc_stn_month_error,c('site_no','month'))
    
    gcmrc_stns_grain_error <- gcmrc_stn_month_error[gcmrc_stns_monthly[,GCMRC_Station_nm:=NULL]]

  # Plot Linear relationship between grain size and model error, computed at a station on monthly averages 
    ggplot(gcmrc_stns_grain_error[GCMRC_Station_nm != 'Colorado River at Potash, UT']) + 
      geom_point(aes(x = avg.p63_monthly, y = avg.cl.model.error), color = 'light blue') + 
      geom_point(aes(x = avg.p63_monthly, y = avg.st.model.error), color = '#5e3c99') + 
      facet_wrap(.~GCMRC_Station_nm) +
      geom_hline(yintercept = 0) +
      season_facet
    
  # Plot grain size and model error variations against month
  gc_grain_size_error <- ggplot(melt(gcmrc_stns_grain_error[GCMRC_Station_nm != 'Colorado River at Potash, UT'],
                                     measure.vars = c('avg.cl.model.error', 'avg.p63_monthly'))) + 
    geom_line(aes(x = month, y = value, color = factor(variable, labels = c('Model error','Fraction < sand')))) +
    scale_color_manual(values = c('light blue', '#5e3c99')) + 
    scale_x_continuous(limits = c(1,12), breaks = c(2,4,6,8,10,12)) +
    facet_wrap(.~factor(GCMRC_Station_nm, levels = c('Green River at Mineral Bottom nr Cynlnds Ntl Park',
                                                     'Colorado River near river mile 30',
                                                     'Colorado River above Little Colorado River near Desert View, AZ',
                                                     'Colorado River near Grand Canyon, AZ'),
                                          labels = c('Green R.\nMineral Bottom, UT','Colorado R.\n30 Mile Station',
                                                     'Colorado R.\n60 Mile Station','Colorado R.\nGrand Canyon, AZ'))) +
    season_facet + 
    geom_hline(aes(yintercept = 0), lty = 'dashed', color = 'black') +
    guides(color = guide_legend(keyheight = 0.1)) + 
    theme(legend.position = 'bottom') +
    labs(
      x = 'Month',
      y = 'Avg monthly model error & \nAvg. monthly P63',
      color = ''
    )
  
  ggsave(ggarrange(co_grain_size_bias_plot_auto, gc_grain_size_error, ncol = 1, labels = c('A','B')), 
         filename = paste0(wd_figures,'gcmrc_error_gs_season.pdf'), 
         width = 4, height = 6)
  
#### --- ####
#### AT-A-STATION MODEL IMPROVEMENT/TESTING ANALYSES ####
#### AT-A-STATION -- BOOTSTRAP INCLUDE TEST ####  
  # Generate a list of sites with N > 100
  ssc_n100_matches <- ssc_station_summary[N_samples > 100]
  ssc_n10_matches <- ssc_station_summary[N_samples > 10]
  setkeyv(ssc_station_summary, c('agency_cd','site_no','Latitude','Longitude', 'drainage_area_km2','begin_date','end_date'))
  setkeyv(ssc_model_cl_iterate_pred, c('agency_cd','site_no','Latitude','Longitude', 'drainage_area_km2','begin_date','end_date'))
  
  # Prepare a data.table of all in situ/satellite matches
  # ls_insitu_bootprep <- ls_insitu_cl[lag_days < 3][
  ls_insitu_bootprep <- ssc_model_cl_iterate_pred[
        # abs_lag_days < 3][
        # ,':='(cluster = NULL, cluster_sel = NULL) # may or may not be required
        ][
        site_no %chin% ssc_n100_matches$site_no]
  
  getBootError <- function(site_no_sel){
    # test
    # site_no_sel <- '05AD007'
    # Calculate model for that station's cluster
    
    # Subset data for holdout site (apply cluster if no cluster already)
    # Select data from given site
    ls_insitu_site_all <- ls_insitu_bootprep[site_no == site_no_sel]
    # Count number of rows
    ls_insitu_site_all_N_samples <- nrow(ls_insitu_site_all)
    # Identify cluster of given site
    site_sel_cluster <- ls_insitu_site_all[1,cluster]
    
    # Select only data from cluster of selected site, but without that site's data
    ls_insitu_cl_spec <- ls_insitu_bootprep[cluster == site_sel_cluster & site_no!=site_no_sel]
    # Generate model for those data, without that site's data (~3 seconds)
    ls_insitu_cl_spec_model <- getClusterModels_lasso(ls_insitu_cl_spec, regressors_all)
    site_cl_predict_model <- ls_insitu_cl_spec_model[[2]]
  
    # Apply CLUSTER model to all at-a-site rows (not just sample)
    cl_model_all_apply <- getSSC_pred(ls_insitu_site_all,regressors_all,site_cl_predict_model)
    
    # Sample at N rows, with replacement (N in c(1,2,3,4,5,10,20,50,100))
    # Progressively sample site
    for(i in c(1,2,3,4,5,10,20,30,40,50,100)){
      print(paste0('sample size: ', i))
      # test
      # i <- 25
      if(i < ls_insitu_site_all_N_samples){
        for(j in 1:20){
          sample_rows <- sample(x = c(1:nrow(ls_insitu_site_all)), size = i, replace = T)
          indiv_site_sample <- ls_insitu_site_all[sample_rows][,holdout25:='in']
          # CLUSTER CALCULATIONS  
          # Calculate CLUSTER median model bias (not log-transformed, just median predicted - median actual)
          cl_model_bias <- cl_model_all_apply[sample_rows,.(median_bias = median(pred_cl - log10_SSC_mgL))]
          
          # Add bias to CLUSTER model for ALL at-a-site rows (not just sample)
          cl_model_sample_bias_add <- cl_model_all_apply[,':='(pred_cl = pred_cl - as.numeric(cl_model_bias))]
          cl_error <- getErrorBias_simple(cl_model_sample_bias_add)[,.(mape_cl = mape_cl_ind, bias_cl = bias_cl)]
          if(i>4){
            # AT-A-SITE CALCULATIONS
            # Calculate SITE model using lasso regression and regression variables (must be 1 variable for the first several Ns)
            indiv_site_model <- getClusterModels_lasso(indiv_site_sample, regressors_all)
            # Apply SITE model to all at-a-site rows (not just sample)
            indiv_model_sample_apply <- getSSC_pred(ls_insitu_site_all,regressors_all,indiv_site_model[[2]])
            # Bias for both CLUSTER and SITE models will be zero for the rows included in the sample, 
            st_error <- getErrorBias_simple(indiv_model_sample_apply[,holdout25:='holdout'])[
                    ,.(mape_st = mape_cl_ind, bias_st = bias_cl)]
          }else{
            st_error <- data.table(mape_st = NA, bias_st = NA)
        }
          # The question is are those rows representative?
          # To test, use getErrorBias to calculate error and bias at each sample size N for both SITE and CLUSTER models
          
          # Save row of [N, SITE_error, SITE_bias, CLUSTER_error CLUSTER_bias] for each N
          if(i == 1 & j == 1){
            error_master <- cbind(site_no_sel,i, st_error, cl_error)
          }else{
            error_all <- cbind(site_no_sel,i, st_error, cl_error)
            error_master <- rbind(error_master, error_all)
          }
          
        }
        # print(error_master)
        }
        
    }
    return(error_master)
  }
      
  for(k in 104:length(unique(ssc_n100_matches$site_no))){
    site_no_sel <- unique(ssc_n100_matches$site_no)[k]
    print(paste0('site', k, ' ', site_no_sel))
    boot_error_site <- getBootError(site_no_sel)
    if(!is.null(boot_error_site)){
    if(k == 1){
      boot_error_master <- boot_error_site
    }else{
      boot_error_master <- rbind(boot_error_master, boot_error_site)
    }
    }
  }
  
  
  boot_error_master_melt <- melt(boot_error_master, measure.vars = c('mape_cl','mape_st'))[,
                              ':='(model = factor(variable, levels = c('mape_cl','mape_st'), 
                                            labels = c('Modified-station','Standalone')))] 
    # Plot bias and error (preferably get them on the same plot, change the boxes to points with error bars)
  boot_error_master_plot <- ggplot(boot_error_master_melt, 
           aes(x = as.factor(i), y = value, fill = model)) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.25) + scale_y_continuous(limits = c(0,3)) + 
      season_facet + 
      scale_fill_manual(values = c('#025159','#F28705')) +
      geom_text(data = data.table(i = c(1), y_pos = 3, model = 'Modified-station', 
                                N_stations = paste0('N stations: ', nrow(boot_bias_master_melt)/440,
                                                    '\nSimulations/station: 440\n(20 simulations/Sample N)')),
                aes(y = y_pos, label = N_stations), size = 2.5, hjust = 0, vjust = 1) +
      theme(legend.position = c(0.8, 0.8)) + 
      labs(
        x = 'N samples',
        y = 'Relative absolute error',
        fill = 'Model type'
      )
    
  boot_bias_master_melt <- melt(boot_error_master, measure.vars = c('bias_cl','bias_st'))[,
                                ':='(model = factor(variable, levels = c('bias_cl','bias_st'), 
                                               labels = c('Modified-station','Standalone')))]
  boot_bias_master_plot <- ggplot(boot_bias_master_melt, 
           aes(x = as.factor(i), y = value, fill = model)) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.25) + scale_y_continuous(limits = c(0,3)) + 
      season_facet + 
      scale_fill_manual(values = c('#025159','#F28705')) +
      theme(legend.position = c(0.8, 0.8)) + 
      labs(
        x = 'N samples',
        y = 'Relative absolute bias',
        fill = 'Model type'
      )
    
  boot_error_bias_plot <- ggarrange(boot_error_master_plot, boot_bias_master_plot, nrow = 2, align = 'hv')
  
  ggsave(boot_error_bias_plot, filename = paste0(wd_figures,'boot_error_bias_plot.pdf'), useDingbats = F, height = 6, width = 4)
  
  for(i in 1:length(unique(ssc_n100_matches$site_no))){
    site_no_sel_plot <- unique(ssc_n100_matches$site_no)[i]
    agency_cd_plot <- as.character(ssc_n100_matches[site_no == site_no_sel_plot, agency_cd][1])
    boot_error_master_site_plot <- ggplot(boot_error_master_melt[site_no_sel ==site_no_sel_plot], 
           aes(x = as.factor(i), y = value, fill = model)) + 
      geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(0,3)) + 
      season_facet + 
      scale_fill_manual(values = c('light blue','#1e8ecf')) +
      # facet_wrap(.~site_no_sel) +
      # geom_text(data = data.table(i = 1, y_pos = 2.5, model = 'Cluster model', N_stations = paste0('N stations = ', nrow(boot_bias_master_melt)/440)),
      #           aes(y = y_pos, label = N_stations), size = 2.5, hjust = 0) +
      theme(legend.position = c(0.8, 0.8)) + 
      labs(
        x = 'N samples',
        y = 'Relative absolute error',
        fill = 'Model type',
        title = paste0(agency_cd_plot, ' ', site_no_sel_plot)
      )
    boot_bias_master_site_plot <- ggplot(boot_bias_master_melt[site_no_sel == site_no_sel_plot], 
           aes(x = as.factor(i), y = value, fill = model)) + 
      geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = c(0,3)) + 
      season_facet + 
      scale_fill_manual(values = c('light blue','#1e8ecf')) +
      # facet_wrap(.~site_no_sel) +
      # geom_text(data = data.table(i = 1, y_pos = 2.5, model = 'Cluster model', N_stations = paste0('N stations = ', nrow(boot_bias_master_melt)/440)),
      #           aes(y = y_pos, label = N_stations), size = 2.5, hjust = 0) +
      theme(legend.position = c(0.8, 0.8)) + 
      labs(
        x = 'N samples',
        y = 'Relative absolute bias',
        fill = 'Model type'
        # title = paste0(agency_cd, ' ', site_no_sel)
      )
    
    boot_error_bias_site_plot <- ggarrange(boot_error_master_site_plot, boot_bias_master_site_plot, nrow = 2, align = 'hv')
    
    ggsave(boot_error_bias_site_plot, filename = paste0(wd_station_standalone,site_no_sel_plot, '_boot_error_bias_plot.pdf'), 
           useDingbats = F, height = 6, width = 4)
  }
  
  
      # Then do this for 100 stations, aggregate the results
  
      
      
  
  
#### AT-A-STATION -- STANDALONE MODELS ####
  ls_insitu_indivprep <- ssc_model_cl_iterate_pred[
    # abs_lag_days < 3][
    # ,':='(cluster = NULL, cluster_sel = NULL) # may or may not be required
    ][
      site_no %chin% ssc_n10_matches$site_no]
  getIndivCalib <- function(site_no_sel){
    # test
    # site_no_sel <- '05AD007'
    
    # Select data from given site
    ls_insitu_site_all <- ls_insitu_indivprep[site_no == site_no_sel]
    # Count number of rows for site
    ls_insitu_site_all_N_samples <- nrow(ls_insitu_site_all)
    # Get cluster of site (for consistency with getClusterModels_lasso function, not used)
    cluster_sel_site <- ls_insitu_site_all$cluster_sel[1]
    # Get agency code for site
    agency_cd_plot <- as.character(ls_insitu_site_all$agency_cd[1])
    
          if(ls_insitu_site_all_N_samples>=10){
            # AT-A-SITE CALCULATIONS
            # Calculate SITE model using lasso regression and regression variables
            indiv_site_model <- getClusterModels_lasso(ls_insitu_site_all, regressors_all)
            # Apply SITE model to all at-a-site rows (not just sample)
            indiv_model_apply <- getSSC_pred(ls_insitu_site_all,regressors_all,indiv_site_model[[2]])
            # Calculate relative error and station bias
            st_error <- getErrorBias_simple(indiv_model_apply[,holdout25:='holdout'])[
              ,.(mape_st = mape_cl_ind, bias_st = bias_cl)]
            # Calculate maximum values to set plot limits (minimum at SSC = 1 mg/L)
            plot_limits = c(1, max(max(indiv_model_apply$SSC_mgL,na.rm = T), max(10^indiv_model_apply$pred_cl,na.rm = T)))
            # Plot actual vs. predicted SSC at station
            st_plot <- ggplot(indiv_model_apply[order(-abs(lag_days))], aes(x = SSC_mgL, y = 10^pred_cl)) + 
                          geom_point(aes(color = as.factor(abs(lag_days)))) +
                          geom_abline(intercept = 0, slope = 1) +
                          # Include relative error as text
                          geom_text(data = st_error, aes(x = 1, y = 0.8*plot_limits[2], 
                                                         label = paste0('rel. error = ',round(mape_st,2))), 
                                                         size = 3, hjust = 0) +
                          scale_y_log10(limits = plot_limits, labels = fancy_scientific) +
                          scale_x_log10(limits = plot_limits, labels = fancy_scientific) +
                          scale_color_brewer(palette = 'Greens', direction = -1) +
              season_facet + 
              theme(legend.position = c(0.8, 0.2)) + 
              labs(
                x = expression(paste(italic("in situ"), " SSC (mg/L)")),
                y = "Satellite-estimated SSC (mg/L)",
                color = '|Lead/Lag (d)|',
                title = paste0(agency_cd_plot, ' ', site_no_sel)
              )
            # Save plot to drive
            ggsave(st_plot, filename = paste0(wd_standalone_figures, 'indiv_calib_plot_', site_no_sel, '.pdf'), width = 4, height = 4, useDingbats = F)
            cv.opt <- coef(indiv_site_model[[2]][[cluster_sel_site]], s = "lambda.1se")
            coef_ex <- cbind(rownames(cv.opt),as.numeric(cv.opt))
            colnames(coef_ex) <- c('variable', site_no_sel)
            
            return(list(coef_ex, cbind(site_no_sel, st_error)))
          }
  }
  
  for(k in 1:length(unique(ssc_n10_matches$site_no))){
    site_no_sel <- unique(ssc_n10_matches$site_no)[k]
    print(paste0('site', k, ' ', site_no_sel))
    site_indiv_calib <- getIndivCalib(site_no_sel)
    site_indiv_calib_coeff <- data.table(site_indiv_calib[[1]])
    site_indiv_calib_error <- data.table(site_indiv_calib[[2]])
    if(!is.null(site_indiv_calib)){
    if(k == 1){
      site_indiv_calib_master_coeff <- site_indiv_calib_coeff
      site_indiv_calib_master_error <- site_indiv_calib_error
    }else{
      site_indiv_calib_master_coeff <- cbind(site_indiv_calib_master_coeff, site_indiv_calib_coeff[,..site_no_sel])
      site_indiv_calib_master_error <- rbind(site_indiv_calib_master_error, site_indiv_calib_error)
    }}
  }
  
  fwrite(site_indiv_calib_master_coeff, file = paste0(wd_standalone_models,'site_indiv_calibration_coefficients.csv'))
  fwrite(site_indiv_calib_master_error, file = paste0(wd_standalone_models,'site_indiv_calibration_errors.csv'))
#### AT-A-STATION -- INDIV. BANDS/RATIOS ####
  
  # regressors_primary <- c('B3.B2','B3.B1','B4.B1','B4_median') # only these regressors
  regressors_primary <- c('B4.B2','B4.B1','B4','B3.B2','B3.B1','B3','B2.B1') # only these regressors
  
  for(i in 1:length(unique(ssc_n10_matches$site_no))){
    site_no_sel <- unique(ssc_n10_matches$site_no)[i] 
    ls_insitu_site_all <- ls_insitu_indivprep[site_no == site_no_sel] # select site data
    cluster <- ls_insitu_site_all$cluster[1] # get cluster
    r2_summary_site <- data.frame(cbind(site_no_sel, cluster))
    sigma_summary_site <- data.frame(cbind(site_no_sel, cluster))
    
    for(j in 1:length(regressors_primary)){
      regressor_sel <- regressors_primary[j]
      lm_cols <- c('log10_SSC_mgL',regressor_sel)
      ssc_regress <- ls_insitu_site_all[,..lm_cols]
      lm_sel <- lm(log10_SSC_mgL~., data = ssc_regress)
      ls_insitu_site_all[,c(paste0(regressor_sel, '_pred'))] <- predict(lm_sel)
      r2_summary_site[1,c(paste0(regressor_sel, '_r2'))] <- glance(lm_sel)[2]
      sigma_summary_site[1,c(paste0(regressor_sel, '_sigma'))] <- glance(lm_sel)[3]
    }
    if(i == 1){
      r2_summary_master <- r2_summary_site
      sigma_summary_master <- sigma_summary_site
    }else{
      r2_summary_master <- rbind(r2_summary_master, r2_summary_site)
      sigma_summary_master <- rbind(sigma_summary_master, sigma_summary_site)
    }
  }
  
  which.regressor <- function(x){return(regressors_primary[which.max(x)])}
  
  r2_summary_master$primary_regressor <- apply(r2_summary_master[,c(-1,-2)], 1, which.regressor)
  
  r2_summary_master$primary_regressor_name <- gsub("_median", "",r2_summary_master$primary_regressor)
  r2_summary_master <- setDT(r2_summary_master)[, primary_regressor_name:= gsub("\\.", "/",primary_regressor)]
  
  r2_summary_master_byCluster <- r2_summary_master[,.N, 
                                                   by = .(cluster, primary_regressor_name)][,
                                                     N_cluster := sum(N), by = .(cluster)
                                                   ][,perc_regressor := N/N_cluster]
  primary_regressor_count_plot <- ggplot(r2_summary_master_byCluster, 
                                          aes(x = primary_regressor_name, y = perc_regressor, fill = cluster)) + 
    geom_bar(stat = 'identity') + 
    facet_wrap(.~paste0('Cluster ', cluster), nrow = 1) +
    season_facet +
    scale_fill_brewer(palette = 'Paired') +
    scale_y_continuous(expand = expand_scale(add = c(0,0.1))) +
    labs(
      x = "Strongest individual regressor",
      y = 'Fraction of stations'
    ) + 
    rotate()
  
  # SAVE FIGURE
  ggsave(primary_regressor_count_plot,filename = paste0(wd_figures,'primary_regressor_count_plot.pdf'), 
         width = 7, height = 2.5, useDingbats = F)
  
  
#### AT-A-STATION -- AUTOCORRELATION IN SSC TIMESERIES ####
  
  # Import all in situ SSC data from USGS stations
  USGS_raw_ssc <- read_csv('1984_2017_ssc_usgs_all.dat')
  
  USGS_raw_ssc$site_no <- sapply(USGS_raw_ssc$site_no, getUSGS_site_no)
  
  # Count measurements at each station
  USGS_raw_station_count <- USGS_raw_ssc %>% subset(!is.na(DAILY_SSC)) %>% 
    group_by(SNAME, site_no) %>% 
    summarise_at(c('DAILY_SSC'),list(count = length)) %>% 
    subset(count > 730)
  
  # Find station info at stations with most measurements
  USGS_long_SSC_site_info <- readNWISsite(USGS_raw_station_count$site_no)
  
  # Initialize dataframe for autocorrelation output
  station_acf_master <- data.frame(matrix(ncol = 4))
  colnames(station_acf_master) <- c('acf','lag','station_nm', 'drainage_area_km2')
  
  # For each station, calculate autocorrelation and export plot of correlation at 1:n day lags
  for(i in 1:nrow(USGS_raw_station_count)){
    station_sel <- USGS_raw_station_count$SNAME[i] # select station
    site_no_sel <- USGS_raw_station_count$site_no[i] # selected station number
    # Find drainage area for station
    drainage_area_sel <- USGS_long_SSC_site_info %>% subset(site_no == site_no_sel) %>% select(drain_area_va) %>% as.numeric()
    # Select station from larger dataset of all stations
    USGS_station_ts_sel <- USGS_raw_ssc %>% subset(SNAME == station_sel & !is.na(DAILY_SSC))
    # Calculate number of days between day and previous row day
    day_diffs <- diff(mdy(USGS_station_ts_sel$datetime))
    # Find days with >1 day gap (i.e., gaps in sequential data)
    ts_gap_indices <- c(1,which(day_diffs != 1),nrow(USGS_station_ts_sel))
    # Find row index of longest gap-free sequence of sampling
    ts_gap_max <- which.max(diff(ts_gap_indices))
    # Find start row and end row endex of longest gap-free sequence of sampling
    ts_longest_start_stop <- c(ts_gap_indices[ts_gap_max]:ts_gap_indices[(ts_gap_max + 1)])
    # Only compute autocorrelation for sites with 2+ years of continuous sampling
    if(length(ts_longest_start_stop) > 730){
      # Subset station data to include only longest continuous period
      USGS_station_ts_longest <- USGS_station_ts_sel[ts_longest_start_stop,]
      # Calculate SSC autocorrelation for up to 16 day lag
      station_acf <- acf(log10(USGS_station_ts_longest$DAILY_SSC + 1), 
                         lag.max = 17,plot = F)$acf[1:17] %>% data.frame() %>% 
        cbind(lag = data.frame(0:16)) %>% 
        mutate(station_nm = station_sel,
               drainage_area = drainage_area_sel * 2.58999)
      
      colnames(station_acf) <- c('acf','lag','station_nm', 'drainage_area_km2')
      
      # Combine with previous station observations to assemble summary dataset
      station_acf_master <- rbind(station_acf, station_acf_master)
      
      # Plot correlations
      # Create column for SSC lag of 1 to 16 days
      USGS_station_ts_lag_corrs <- USGS_station_ts_longest %>% select(DAILY_SSC) %>%
        mutate(lag1 = lag(DAILY_SSC,1),
               lag2 = lag(DAILY_SSC,2),
               lag3 = lag(DAILY_SSC,3),
               lag4 = lag(DAILY_SSC,4),
               lag5 = lag(DAILY_SSC,5),
               lag6 = lag(DAILY_SSC,6),
               lag7 = lag(DAILY_SSC,7),
               lag8 = lag(DAILY_SSC,8),
               lag9 = lag(DAILY_SSC,9),
               lag10 = lag(DAILY_SSC,10),
               lag11 = lag(DAILY_SSC,11),
               lag12 = lag(DAILY_SSC,12),
               lag13 = lag(DAILY_SSC,13),
               lag14 = lag(DAILY_SSC,14),
               lag15 = lag(DAILY_SSC,15),
               lag16 = lag(DAILY_SSC,16)
        ) %>% melt(measure.vars = paste0('lag',c(1:16))) %>%
        mutate(lag = factor(paste0('Lag = ', gsub('lag',"",variable), ' d'),
                            levels = paste0('Lag = ', c(1:16), ' d'), ordered = T))
      # Change the lag column to have format Lag = X d
      station_acf <- station_acf %>% mutate(lag = factor(paste0('Lag = ', gsub('lag',"",lag), ' d'),
                                                         levels = paste0('Lag = ', c(0:16), ' d'), ordered = T))
      
      # Plot each lag SSC against that day SSC -- 16 plots in a grid
      USGS_station_ts_lags_plot <- ggplot(data = USGS_station_ts_lag_corrs %>% subset(DAILY_SSC != 0 & value != 0),
                                          aes(x = DAILY_SSC, y = value)) + 
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = '#D9583B') +
        geom_text(data = station_acf %>% subset(lag != 'Lag = 0 d'), 
                  aes(x = min(USGS_station_ts_lag_corrs$DAILY_SSC, na.rm = T), 
                      y = max(USGS_station_ts_lag_corrs$DAILY_SSC, na.rm = T),
                      label = paste0('R = ', round(acf,2))), hjust = 0, vjust = 1, size = 2) +
        scale_x_log10(labels = fancy_scientific) + 
        scale_y_log10(labels = fancy_scientific) +
        facet_wrap(.~lag) + 
        season_facet +
        labs(x = 'SSC (mg/L)',
             y = 'SSC lag (mg/L)',
             title = station_sel)
      
      # Save lag SSC plots
      ggsave(USGS_station_ts_lags_plot, filename = paste0(wd_autocorrelation, 'Lag_correlations_', site_no_sel, '.pdf'),
             width = 5, height = 5, useDingbats = F)
    }
  }
  
  # Break up drainage areas into categories
  station_acf_2 <- station_acf_master %>%
    # mutate(drainage_bins = cut(drainage_area_km2, breaks = c(0,1000,5000,10000,50000,100000,1e8),
    #                            labels = c('0-1,000','1,000-5,000','5000-10,000','10,000-50,000','50,000-100,000','>100,000'))) %>%
    mutate(drainage_bins = cut(drainage_area_km2, breaks = c(0,10000,100000,1e8),
                               labels = c('0-10,000','10,000-100,000','>100,000')))
  # Summarize autocorrelation by drainage area 
  station_acf_master_summary <- station_acf_master_summary %>%
    group_by(drainage_bins, lag) %>%
    summarise_at('acf', mean, na.rm = T) %>% subset(!is.na(drainage_bins))
  
  # Plot summarized autocorrelation vs. lag
  ssc_autocorrelation_drainage_plot <- ggplot(station_acf_master_summary,
                                              aes(x = lag, y = acf, color = drainage_bins, group = drainage_bins)) + 
    geom_line() + 
    scale_color_brewer(palette = 'Paired') +
    scale_y_continuous(limits = c(0,1)) +
    season_facet + 
    theme(legend.position = 'right') + 
    labs(
      x = 'Lag (days)',
      y = 'Autocorrelation in SSC',
      color = parse(text = paste("'Drainage area (km'","^2","*')'", sep = ""))
    )
  ggsave(ssc_autocorrelation_drainage_plot, filename = paste0(wd_autocorrelation,'ssc_autocorrelation_drainage_plot.pdf'),
         width = 5, height = 5, useDingbats = F)
  
  # Plot individual autcorrelation for each station, colored by drainage area category
  ssc_autocorrelation_indiv_drainage_plot <- ggplot(station_acf_2,
                                                    aes(x = lag, y = acf, color = drainage_bins, group = station_nm)) + 
    # geom_line() + 
    geom_smooth(aes(group = drainage_bins)) +
    scale_color_brewer(palette = 'Paired') +
    scale_y_continuous(limits = c(0,1)) +
    season_facet + 
    theme(legend.position = 'right') + 
    labs(
      x = 'Lag (days)',
      y = 'Autocorrelation in SSC',
      color = parse(text = paste("'Drainage area (km'","^2","*')'", sep = ""))
    )
  
  ggsave(ssc_autocorrelation_indiv_drainage_plot, filename = paste0(wd_autocorrelation,'ssc_autocorrelation_indiv_drainage_plot.pdf'),
         width = 5, height = 5, useDingbats = F)
  
  
  
  
  
#### --- ####
  
  
#### MAP AVERAGE SSC ####
    # Generate ssc map plots - average ssc and error
    getSSC_avg_error_plots <- function(data,model_type, region){
      if(model_type == 'bias_st'){
        filename_paste <-  'st'
      } else if(model_type == 'bias_cl'){
        filename_paste <-  'cl'
      }else{
        filename_paste <- 'gl'
      }
    
      station_bias_map_plot <- ggplot(data = ssc_station_summary[N_samples > 12], 
                                      aes(x = Longitude, y = Latitude, 
                                          fill = cut(!!sym(model_type), breaks = c(-1,0.5,1, 2, 1e7),
                                                     labels = error_fractions)
                                      )) + 
        geom_map(data = world_map, map = world_map, aes(map_id=region, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) + # world map
        geom_point(
          size = 2.5,
          # aes(size = cut(N_samples, breaks = c(0,20, 100, 500, 1e5), 
          #                         labels = c('1-20', '21-100','101-500', '> 500'))),
          pch = 21, color = 'black', stroke = 0.2) +
        scale_size_discrete(range = c(1.5,4)) +
        theme_bw() + 
        scale_fill_manual(values = error_colors) + scale_color_brewer(palette = 'PuOr') +
        # coordinates need to extend to full extent of polygon data or it gets cut off
        scale_y_continuous(limits = whemi_limits_y) + scale_x_continuous(limits = whemi_limits_x) + 
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(size = 0.5),
          text = element_text(size=12),
          axis.text = element_text(size = 12)
        ) + 
        labs(
          x = "Longitude",
          y = "Latitude",
          size = 'N samples',
          fill = 'Station relative bias'
        )
      
      tai_bias_subplot <- station_bias_map_plot + 
        scale_y_continuous(limits = tai_limits_y) + 
        scale_x_continuous(limits = tai_limits_x, breaks = c(120, 121, 122)) + 
        theme(legend.position = "none",
              axis.title = element_blank()) + labs(main = 'Taiwan')
      
      anno_bias_tai <- annotation_custom(grob = ggplotGrob(tai_bias_subplot), # joins plot with centroid coordinates
                                    xmin = -140, ymin = -20, 
                                    xmax = -100, ymax = 15)
      SSC_error_inset_map_plot <- station_bias_map_plot + anno_bias_tai
      
      ggsave(SSC_error_inset_map_plot, filename = paste0(wd_exports,'whemi_', filename_paste, '_avg_ssc_map.pdf'), useDingbats = F, width = 8, height = 7)
      ggsave(SSC_error_inset_map_plot, filename = paste0(wd_exports,'whemi_', filename_paste, '_avg_ssc_map.png'),width = 8, height = 7)
      
      spec_limits <- c(0,4) # set limits for gradient
      
      ssc_map_whemi_plot <- ggplot(data = ssc_station_summary[N_samples > 12], 
                                      aes(x = Longitude, y = Latitude, 
                                          fill = log10_SSC_mgL
                                      )) + 
        geom_map(data = world_map, map = world_map, aes(map_id=region, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) + # world map
        geom_point(
          size = 2.5,
          # aes(size = cut(N_samples, breaks = c(0,20, 100, 500, 1e5), 
          #                         labels = c('1-20', '21-100','101-500', '> 500'))),
          pch = 21, color = 'black', stroke = 0.2) +
        scale_size_discrete(range = c(1.5,4)) +
        theme_bw() + 
        # scale_fill_manual(values = error_colors) + scale_color_brewer(palette = 'PuOr') +
        scale_fill_gradientn(colors = c('navy', 'dark green', 'yellow', 'red'),
                             limits = c(0,3), oob = squish, labels = c('0','10','100','1000')) +
        # coordinates need to extend to full extent of polygon data or it gets cut off
        scale_y_continuous(limits = whemi_limits_y) + scale_x_continuous(limits = whemi_limits_x) + 
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(size = 0.5),
          text = element_text(size=12),
          axis.text = element_text(size = 12)
        ) + 
        labs(
          x = "Longitude",
          y = "Latitude",
          fill = 'SSC (mg/L)'
        )
      
      
      
      tai_ssc_subplot <- ssc_map_whemi_plot + 
        scale_y_continuous(limits = tai_limits_y) + 
        scale_x_continuous(limits = tai_limits_x, breaks = c(120, 121, 122)) + 
        theme(legend.position = "none",
              axis.title = element_blank())
      
      anno_ssc_tai <- annotation_custom(grob = ggplotGrob(tai_ssc_subplot), # joins plot with centroid coordinates
                                    xmin = -140, ymin = -20, 
                                    xmax = -100, ymax = 15)
      SSC_map_inset_plot <- ssc_map_whemi_plot + anno_ssc_tai
      
      ggsave(SSC_map_inset_plot, filename = paste0(wd_exports,'whemi_', filename_paste, '_avg_ssc_map.pdf'), useDingbats = F, width = 8, height = 7)
      ggsave(SSC_map_inset_plot, filename = paste0(wd_exports,'whemi_', filename_paste, '_avg_ssc_map.png'), width = 8, height = 7)
      
      
      cluster_map_whemi_plot <- ggplot(data = ssc_station_summary[N_samples > 12], 
                                   aes(x = Longitude, y = Latitude, 
                                       fill = as.factor(cluster_sel)
                                   )) + 
        geom_map(data = world_map, map = world_map, aes(map_id=region, x = long, y = lat), color = "grey30", fill = "grey95", size = 0.25) + # world map
        geom_point(
          size = 2.5,
          # aes(size = cut(N_samples, breaks = c(0,20, 100, 500, 1e5), 
          #                         labels = c('1-20', '21-100','101-500', '> 500'))),
          pch = 21, color = 'black', stroke = 0.2) +
        scale_size_discrete(range = c(1.5,4)) +
        theme_bw() + 
        # scale_fill_manual(values = error_colors) + scale_color_brewer(palette = 'PuOr') +
        scale_fill_brewer(palette = 'Paired') +
        # coordinates need to extend to full extent of polygon data or it gets cut off
        scale_y_continuous(limits = whemi_limits_y) + scale_x_continuous(limits = whemi_limits_x) + 
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(size = 0.5),
          text = element_text(size=12),
          axis.text = element_text(size = 12)
        ) + 
        labs(
          x = "Longitude",
          y = "Latitude",
          fill = 'Cluster group'
        )
      
      
      
      tai_cluster_subplot <- cluster_map_whemi_plot + 
        scale_y_continuous(limits = tai_limits_y) + 
        scale_x_continuous(limits = tai_limits_x, breaks = c(120, 121, 122)) + 
        theme(legend.position = "none",
              axis.title = element_blank())
      
      anno_tai_cluster <- annotation_custom(grob = ggplotGrob(tai_cluster_subplot), # joins plot with centroid coordinates
                                    xmin = -140, ymin = -20, 
                                    xmax = -100, ymax = 15)
      cluster_map_inset_plot <- cluster_map_whemi_plot + anno_tai_cluster
      
      ggsave(cluster_map_inset_plot, filename = paste0(wd_exports,'whemi_', filename_paste, '_cluster_map.pdf'), useDingbats = F, width = 8, height = 7)
      ggsave(cluster_map_inset_plot, filename = paste0(wd_exports,'whemi_', filename_paste, '_cluster_map.png'), width = 8, height = 7)
      
      
      
      SSC_avg_err_comb_plot <- ggarrange(plotlist = list(SSC_map_inset_plot, SSC_error_inset_map_plot), ncol = 2, align = 'hv')
      # SSC_avg_err_comb_plot <- ggarrange(plotlist = list(SSC_cluster_map_plot, SSC_map_plot, SSC_error_map_plot), ncol = 3, align = 'hv')
      ggsave(SSC_avg_err_comb_plot, filename = paste0(wd_exports,'whemi_', filename_paste, 'ssc_avg_err_comb_map.pdf'), useDingbats = F, width = 12, height = 5)
      ggsave(SSC_avg_err_comb_plot, filename = paste0(wd_exports,'whemi_', filename_paste, 'ssc_avg_err_comb_map.png'),width = 12, height = 5)
      
      return(list(cluster_map_inset_plot,SSC_map_inset_plot,SSC_error_inset_map_plot))
    }

  ssc_map_gl_plots <- getSSC_avg_error_plots(data = ssc_station_summary, region = 'na',model_type = 'bias_gl')
  ssc_map_cl_plots <- getSSC_avg_error_plots(data = ssc_station_summary, region = 'na',model_type = 'bias_cl')
  ssc_map_st_plots <- getSSC_avg_error_plots(data = ssc_station_summary, region = 'na',model_type = 'bias_st')
  
    ggsave(ssc_map_cl_plots[[1]], width = 7.5, height = 6, filename = paste0(wd_figures, 'cluster_map_n', cluster_n_best,'.pdf'), useDingbats = F)
    ggsave(ggarrange(ssc_map_cl_plots[[2]], ssc_map_cl_plots[[3]], nrow = 1, align = 'hv', labels = c('A','B')), 
           width = 12, height = 5, filename = paste0(wd_figures, 'cluster_maps_combined.pdf'))

#### MAKE REGRESSION PLOTS ####

    ssc_pred <- ssc_model_cl_iterate_pred
    cl_colors <- brewer.pal(name = 'Paired',n=12)[c(2,4,6,8,11)] # for five clusters
    
    ssc_density_global <- get_sscPlot(ssc_pred,"byGlobal",'yes', 'no')
    # ssc_plot_global_hold <- get_sscPlot(ssc_pred,"byGlobal",'no','yes')
    ssc_plot_global <- get_sscPlot(ssc_pred,"byGlobal",'overlay','no')

    # ggsave(ssc_density_global, filename = paste0(wd_exports,'ssc_density_global.pdf'))
    ggsave(ssc_plot_global, filename = paste0(wd_exports,'ssc_plot_global_lasso.pdf'), width = 6.5, height = 7)

    ssc_density_cluster <- get_sscPlot(ssc_pred,"byCluster",'yes','no')
    # ssc_plot_cluster_hold <- get_sscPlot(ssc_pred,"byCluster",'no','yes')
    ssc_plot_cluster <- get_sscPlot(ssc_pred,"byCluster",'overlay','no')
    
    
    # ssc_plot_cluster_co <- get_sscPlot(co_models[[1]],"byCluster",'no','no')
    # ssc_plot_cluster_hold_am <- get_sscPlot(ssc_pred_am,"byCluster",'no','yes')
    
    ssc_density_cluster_p63 <- get_sscPlot(gs_avg_model[[1]],"byCluster",'yes','no')
    # ssc_plot_cluster_hold_p63 <- get_sscPlot(gs_avg_model[[1]],"byCluster",'no','yes')
    ssc_plot_cluster_p63 <- get_sscPlot(gs_avg_model[[1]],"byCluster",'overlay','no')
    
    # ggsave(ssc_density_cluster, filename = 'paste0(wd_exports,ssc_density_cluster.pdf'))
    # ggsave(ssc_plot_cluster_hold, filename = paste0(wd_exports,'ssc_plot_cluster_hold_lasso.pdf'), width = 6.5, height = 7)
    ggsave(ssc_plot_cluster, filename = paste0(wd_exports,'ssc_plot_cluster_lasso.pdf'), width = 6.5, height = 7)
    # ggsave(ssc_res_cluster, filename = paste0(wd_exports,'ssc_res_cluster.pdf'))
    # ggsave(ssc_plot_cluster_hold_p63, filename = paste0(wd_exports,'ssc_plot_cluster_hold_lasso_p63.pdf'), width = 6.5, height = 7)
    ggsave(ssc_plot_cluster_p63, filename = paste0(wd_exports,'ssc_plot_cluster_lasso_p63.pdf'), width = 6.5, height = 7)
    
    ssc_density_site <- get_sscPlot(ssc_pred,"bySite",'yes','no')
    ssc_plot_site <- get_sscPlot(ssc_pred,"bySite",'overlay','no')
    # ssc_density_site_p63 <- get_sscPlot_p63(ssc_station_pred,"bySite",'yes','no')
    # ssc_plot_site_p63 <- get_sscPlot_p63(ssc_station_pred,"bySite",'no','no')
    
    # ggsave(ssc_density_site, filename = paste0(wd_exports,'ssc_density_site.pdf'))
    ggsave(ssc_plot_site, filename = paste0(wd_exports,'ssc_plot_site_lasso.pdf'), width = 6.5, height = 7)
    # ggsave(ssc_plot_site_p63, filename = paste0(wd_exports,'ssc_plot_site_lasso_p63.pdf'), width = 6.5, height = 7)

    print(ssc_model_all_errorbias)
    all_regression_plots <- ggarrange(plotlist = list(ssc_plot_global,ssc_plot_cluster, ssc_plot_site, ssc_plot_cluster_p63), 
                                        labels = c('Base\nrel. error = 0.97',
                                                   # \nrel. bias = 0.76', 
                                                 'Global\nrel. error = 0.73',
                                                 # \nrel. bias = 0.54', 
                                                 'Station\nrel. error = 0.49',
                                                 # \nrel. bias = 0.08', 
                                                 'Grain-size\nrel. error = 0.65'),
                                                 # \nrel. bias = 0.58'),
                                      ncol = 2, nrow = 2, align = 'hv', common.legend = F, label.y = 1, label.x = 0.22, hjust = 0)
    
    ggsave(all_regression_plots, width = 8, height = 11, filename = paste0(wd_figures,'all_regression_plots.pdf'), useDingbats = F)
    ggsave(all_regression_plots, width = 8, height = 11, filename = paste0(wd_figures,'all_regression_plots.png'))

    all_regression_density_plots <- ggarrange(plotlist = list(ssc_density_global,ssc_density_cluster, ssc_density_site, ssc_density_cluster_p63), 
                                              labels = c('Base\nrel. error = 0.97',
                                                         # \nrel. bias = 0.76', 
                                                         'Global\nrel. error = 0.73',
                                                         # \nrel. bias = 0.54', 
                                                         'Station\nrel. error = 0.49',
                                                         # \nrel. bias = 0.08', 
                                                         'Grain-size\nrel. error = 0.65'),
                                              # \nrel. bias = 0.58'),
                                      ncol = 2, nrow = 2, align = 'hv', common.legend = F, label.y = 1, label.x = 0.22, hjust = 0)
    
    ggsave(all_regression_density_plots, width = 8, height = 11, filename = paste0(wd_figures,'all_regression_density_plots.pdf'), useDingbats = F)
    ggsave(all_regression_density_plots, width = 8, height = 11, filename = paste0(wd_figures,'all_regression_density_plots.png'))
    

#### CALCULATE STATISTICS OF LS-IN SITU DATASET ####
    # Total number of raw Landsat samples
    n_landsat_all <- ls_clean[,.(N_ls_samples = .N)]
    # Total number of in situ samples
    n_insitu_all <- all_acy_insitu_raw[!is.na(log10(SSC_mgL)),.(N_insitu_samples = .N)]
    # Number of in situ daily means
    n_insitu_daily_mean <- all_acy_insitu_daily_mean[!is.na(log10(SSC_mgL)),.(N_insitu_samples = .N)]
    # Total number of samples by absolute value of lead/lag
    n_ls_insitu_lag <- ls_insitu_raw[!is.na(log10_SSC_mgL),.(N_insitu_samples = .N), by = abs_lag_days][order(abs_lag_days)]
    # Number of in situ-satellite pairs, by lag days
    n_match_samples_byLag <- setkey(ls_insitu_raw[!is.na(B1) & !is.na(log10_SSC_mgL),.(N_samples = .N), 
                                                  by = lag_days],lag_days)
    n_match_samples_lag_lt_3 <- n_match_samples_byLag[abs(lag_days) < 3,.(N_insitu_samples = sum(N_samples))]
    # Number of in situ-landsat matched samples, per site
    n_insitu_samples_bySite <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3,.(N_insitu_samples = .N), by = .(site_no, agency_cd)]
    # by agency
    n_insitu_samples_bySite_agency <- n_insitu_samples_bySite[,.('N in situ samples' = sum(N_insitu_samples)), by = agency_cd]
    # Number of sites
    n_insitu_sites <- length(unique(n_insitu_samples_bySite$site_no))
    # by agency
    n_insitu_sites_byAgency <- n_insitu_samples_bySite[,.('N stations' = .N), by = agency_cd]
    # by > 10 samples
    n_insitu_sites_byAgency_n10 <- n_insitu_samples_bySite[N_insitu_samples > 9,.('N stations (>= 10 samples)' = .N), by = agency_cd]
    # Number of landsat-in situ matched samples with P63
    n_ls_insitu_p63 <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3 & !is.na(p63),
                                     .('N P63' = .N)]
    # by agency 
    n_ls_insitu_p63_byAgency <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3 & !is.na(p63),
                                              .('N P63' = .N), by = agency_cd]
    # Number of in situ-landsat stations with POC
    n_ls_insitu_P63_sites_byAgency <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3 & !is.na(p63),
                                                    .('N P63' = .N), by = .(agency_cd, site_no)][
                                                      , .('N P63 stations' = .N), by = agency_cd
                                                      ]
    # Number of in situ-landsat matched samples with POC
    n_ls_insitu_POC <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3 & !is.na(POC_mgL),
                                     .(N_insitu_samples = .N)]
    n_ls_insitu_POC_byAgency <- rbind(ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3 & !is.na(POC_mgL),
                                                    .('N POC' = .N), by = agency_cd],
                                      data.table(agency_cd = 'WRA',
                                                 'N POC' = sum(tawian_carbon_site_summary$POC_N)))
    # Number of in situ-landsat stations with POC
    n_ls_insitu_POC_sites_byAgency <- rbind(ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3 & !is.na(POC_mgL),
                                                          .('N POC' = .N), by = .(agency_cd, site_no)][
                                                            , .('N POC stations' = .N), by = agency_cd
                                                            ],
                                            data.table(agency_cd = 'WRA',
                                                       'N POC stations' = length(unique(tawian_carbon_site_summary$River))))
    
    # Number of in situ-landsat matched samples by sensor
    n_ls_insitu_sensor <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3,
                                        .(N_insitu_samples = .N), by = sensor] 
    # by agency
    n_ls_insitu_sensor_byAgency <- ls_insitu_raw[!is.na(log10_SSC_mgL) & abs_lag_days < 3,
                                                 .('N Landsat 5' = sum(sensor == 'Landsat 5'),
                                                   'N Landsat 7' = sum(sensor == 'Landsat 7')), by = agency_cd]                                
    
    n_samples_master_table <- merge(merge(merge(merge(n_insitu_sites_byAgency[
      n_insitu_sites_byAgency_n10, on = 'agency_cd'][
        n_insitu_samples_bySite_agency, on = 'agency_cd'][
          n_ls_insitu_sensor_byAgency, on = 'agency_cd'],
      n_ls_insitu_P63_sites_byAgency, by = 'agency_cd', all = T),
      n_ls_insitu_p63_byAgency, by = 'agency_cd', all = T),
      n_ls_insitu_POC_sites_byAgency, by = 'agency_cd', all = T),
      n_ls_insitu_POC_byAgency, by = 'agency_cd', all = T)[
        order(-`N in situ samples`)
        ]
    fwrite(n_samples_master_table, file = paste0(wd_exports,'n_samples_master_table.csv'))
    
    ssc_sites_many_n <- n_sat_samples[site_no %chin% all_acy_sts]
    
    n_insitu_samples_byAgency <- setkey(ls_insitu_raw[!is.na(log10_SSC_mgL),.(N_samples_agency = .N), by = agency_cd],agency_cd)
    # holdout_sts <- c(5559600,"Madeira_PortoVelho",9364500,9315000,"Shuangyuanda Bridge")
    holdout_sts <- c("Madeira_PortoVelho", '12340500', '05325000', '06115200')
    
    lag_by_agency_plot <- ggplot(setkey(ls_insitu_raw,agency_cd)[
      n_insitu_samples_byAgency], 
      aes(x = lag_days, 
          group = agency_cd, fill = agency_cd)) + 
      geom_bar() +
      season_facet + 
      facet_wrap(.~reorder(agency_cd,-N_samples_agency), scales = 'free_y') +
      scale_fill_brewer(palette = 'Paired') + 
      theme(legend.position = 'bottom') + 
      labs(x = 'Image-in situ SSC lag (d)',
           y = 'N Paired Samples',
           fill = 'Sampling Agency')
    
    ggsave(lag_by_agency_plot, filename = paste0(wd_exports,'landsat_ssc_lag_by_agency_plot.pdf'), useDingbats = F, width = 5, height = 5)
#### GENERATE EXAMPLE MAPS WITH LANDSAT IMAGES ####
    # Select stations from near Missouri/Mississippi R. confluence
    ex_stns <- c('06935965','06935972','06937000', # Missouri
                      '05587060','05587455','05587500','05587550','05586300', # Upper Mississippi, Illinois
                      '07005500','07010000') # Mississippi
    # Possible dates '2007-08-26'
    # Possible dates '2006-08-23'
    # Possible dates '2006-09-08'
    # Possible dates '2008-04-06'
    # Possible dates '2008-07-27'
    
    ex_stns <- c('07022000', # Mississippi
                 '03612600', '03612500', # Ohio
                 '03609750','03609750','03438220' # Tennessee, Cumberland
                 )
    
    # Possible dates: '2000-09-16' and '2000-09-17' ls 5 and 7
    # Possible dates: '2000-04-25' and '2000-04-26' ls 5 and 7
    # Possible dates: '2000-07-14' and '2000-07-15' ls 5 and 7
    
    ex_stns <- c('06800500', # Elkhorn
                 '06805000', # Salt
                 '06805500','06796000','06774000','06770500','06768000','06765698', # Platte
                 '06793000','06791150', '06785000', '06792000', # Loup
                 '06610000', '06610505','06807000' # Missouri
                 )
    
    # Possible dates: '1987-07-30' ls 5
    # Possible dates: '2003-07-10' ls 5
    # Possible dates: '2003-05-07' ls 5, regional
    # Possible dates: '1993-07-30' ls 5
    
    ex_stns <- c('06329500' # Yellowstone, Sidney
                 )
    
    
    # Select in situ data from near Missouri/Mississippi R. confluence
    ex_insitu_raw <- all_acy_insitu_daily_mean[site_no %chin% ex_stns]
    
    ex_insitu_raw_plot <- ex_insitu_raw[
      sample_dt > ymd('2007-03-01') &
      sample_dt < ymd('2007-11-01')
      ]
    
    # Plot data from selected timeframe
    ex_landsat_lag_range <- data.frame(station_nm = unique(ex_insitu_raw_plot$station_nm), 
                                       site_no = unique(ex_insitu_raw_plot$site_no))
    ex_insitu_timeseries_plot <- ggplot(data = ex_insitu_raw_plot, 
                                             aes(x = sample_dt, y = SSC_mgL, 
                                             color = reorder(station_nm, as.numeric(site_no)), group = station_nm)) + 
    
      geom_line() +
      geom_point() +
    
      season_facet + 
    
      theme(legend.position = 'right') +
      labs(
        x = '',
        y = 'SSC (mg/L)',
        color = 'USGS Station Name'
      )
    print(ggplotly(ex_insitu_timeseries_plot))
    # Subset of points for display
    ex_ls_insitu <- ssc_pred[
      site_no == '06329500' & 
        sensor == 'Landsat 5' &
        # year(sample_dt) == 2007 &
        sample_dt > ymd('2007-06-10') &
        sample_dt < ymd('2007-07-20') &
        sample_dt != ymd('2007-07-06') &
        abs(lag_days) < 3]
    # Plot subset of points regression
    ssc_plot_cluster_ex <- get_sscPlot(ex_ls_insitu[order(-abs(lag_days))],
      "bySite",'no','no') + geom_point(aes(color = abs(lag_days)))
    print(ex_ls_insitu[,.(SSC_mgL, 10^pred_gl, 10^pred_cl, 10^pred_st)])
   fwrite(ex_ls_insitu, file = paste0(wd_exports,'yellowstone_sidney_2007.csv'))

#### --- ####


#### FUNCTIONS ####
  # Generates random holdout set for each cluster
  getHoldout <- function(datatable){
  setDT(datatable)
  rows <- c()
  for(cl_sel in unique(datatable$cluster_sel)){
    cluster_sub <- which(datatable[,cluster_sel] == cl_sel)
    rows <- c(rows,sample(cluster_sub,round(length(cluster_sub)*0.25)))
    # print(length(rows))
  }
  
  datatable[,holdout25:='in'][rows,holdout25:='holdout']
  return(datatable)
  
} 

  # Use lasso regression to calculate base, cluster, station calibration models for each cluster
  getModels_lasso <- function(data, regressors){
  
  # for testing
  # data <- ls_insitu_cl_sample
  
  lm_data <- na.omit(setDT(data)[,
                         ':='(pred_gl = NA,res_gl = NA,pred_cl = NA,
                              res_cl = NA, pred_st = NA, res_st = NA,
                              ssc_subset = cluster_sel
                         )][!is.na(ssc_subset)], cols = c('log10_SSC_mgL',regressors))
  # clusters
  subsets <- unique(lm_data$ssc_subset)
  
  n_clusters <- length(subsets)
  cluster_funs <- list(rep(NA,length(subsets)))
  site_funs <- list(rep(NA,length(subsets)))
  
  
  for(i in subsets){ # for individual cluster models
    regressors_sel <- regressors[-which(regressors == 'site_no')]
    lm_data_lm <- lm_data[ssc_subset == i] # only chooses sites within cluster
    # lm_data_hold <- lm_data_lm[-which(lm_data_lm$site_no %in% holdout_sts),] # for eliminating certain sites
    lm_data_hold <- lm_data_lm[which(lm_data_lm$holdout25 == 'in'),]
    glm_y <- as.matrix(lm_data_hold[,log10_SSC_mgL])
    glm_x <- as.matrix(lm_data_hold[,..regressors_sel])
    
    ssc_lm <- cv.glmnet(x = glm_x, y = glm_y, family = 'gaussian', type.measure = "mse", nfolds = 10)
    cv.opt <- coef(ssc_lm, s = "lambda.1se")
    coef_ex <- cbind(rownames(cv.opt),as.numeric(cv.opt))
    colnames(coef_ex) <- c('variable', 'value')
    
    write.table(coef_ex, sep = ",", file = paste0(wd_exports,'cluster_ncl',n_clusters,'_', i,'_lasso_fit_coeff.csv'), row.names = F)
    glm_x <- NA
    glm_y <- NA
    
    cluster_funs[[i]] <- ssc_lm
    glm_pred <- predict(ssc_lm, newx = as.matrix(lm_data_lm[,..regressors_sel]), s = "lambda.1se")
    lm_data$pred_cl[which(lm_data$ssc_subset == i)] <- glm_pred
    
    # lm_data$res_cl[which(lm_data$ssc_subset == i)] <- resid(ssc_lm)
  }
  
  lm_data$ssc_subset <- lm_data$cluster_sel # sites
  lm_data$site_code <- lm_data$site_no
  subsets <- unique(lm_data$ssc_subset)
  for(i in subsets){ # for individual cluster models
    
    regressors_sel <- regressors
    lm_data_lm <- subset(lm_data, ssc_subset == i) # only chooses sites within cluster
    # lm_data_hold <- lm_data_lm[-which(lm_data_lm$site_no %in% holdout_sts),] # for eliminating certain sites
    # lm_data_hold <- lm_data_lm[which(lm_data_lm$holdout25 == 'in'),]
    lm_data_hold <- lm_data_lm # need to make sure all sites are included in lm_data_hold in order to use holdout
    
    glm_y <- as.matrix(lm_data_hold[,log10_SSC_mgL])
    glm_x <- model.matrix( ~ ., lm_data_hold[,..regressors])
    
    ssc_lm <- cv.glmnet(x = glm_x,y = glm_y, family = 'gaussian', type.measure = "mse", nfolds = 10)
    cv.opt <- coef(ssc_lm, s = "lambda.1se")
    coef_ex <- cbind(rownames(cv.opt),as.numeric(cv.opt))
    colnames(coef_ex) <- c('variable', 'value')
    
    write.table(coef_ex, sep = ",", file = paste0(wd_exports, 'cluster_ncl',n_clusters,'_', i,'site_lasso_fit_coeff.csv'), row.names = F)
    
    site_funs[[i]] <- ssc_lm
    glm_x <- NA
    lm_data_hold <- NA
    
    glm_x_new <- model.matrix( ~ ., lm_data_lm[,..regressors])
    glm_pred <- predict(ssc_lm, newx = glm_x_new, s = "lambda.1se")
    lm_data$pred_st[which(lm_data$ssc_subset == i)] <- glm_pred
    # lm_data$res_cl[which(lm_data$ssc_subset == i)] <- resid(ssc_lm)
  }
  
  lm_data$ssc_subset <- 1
  subsets <- unique(lm_data$ssc_subset) # global - only one value for subset
  for(i in subsets){ # for individual models
    regressors_sel <- regressors[-which(regressors == 'site_no')]
    lm_data_lm <- lm_data[ssc_subset == i] # only chooses sites within cluster
    # lm_data_hold <- lm_data_lm[-which(lm_data_lm$site_no %in% holdout_sts),] # for eliminating certain sites
    lm_data_hold <- lm_data[which(lm_data$holdout25 == 'in'),]
    glm_y <- as.matrix(lm_data_hold[,log10_SSC_mgL])
    glm_x <- as.matrix(lm_data_hold[,..regressors_sel])
    
    ssc_lm <- cv.glmnet(x = glm_x,y = glm_y, family = 'gaussian', type.measure = "mse", nfolds = 10)
    cv.opt <- coef(ssc_lm, s = "lambda.1se")
    coef_ex <- cbind(rownames(cv.opt),as.numeric(cv.opt))
    colnames(coef_ex) <- c('variable', 'value')
    
    write.table(coef_ex, sep = ",", file = paste0(wd_exports,'global_ncl',n_clusters,'_', i, i,'_lasso_fit_coeff.csv'), row.names = F)
    
    glm_pred <- predict(ssc_lm, newx = as.matrix(lm_data_lm[,..regressors_sel]), s = "lambda.1se")
    lm_data$pred_gl[which(lm_data$ssc_subset == i)] <- glm_pred
  }
  return(list(lm_data, cluster_funs, site_funs))
}
  
  # Use lasso regression to calculate only cluster calibration models for each cluster
  getClusterModels_lasso <-  function(data, regressors){
  
  # for testing
  # data <- indiv_site_sample
  
  lm_data <- na.omit(setDT(data)[,
                         ':='(pred_gl = NA,res_gl = NA,pred_cl = NA,
                              res_cl = NA, pred_st = NA, res_st = NA,
                              ssc_subset = cluster_sel
                         )][!is.na(ssc_subset)], cols = c('log10_SSC_mgL',regressors))
  # clusters
  subsets <- unique(lm_data$ssc_subset)
  
  cluster_funs <- list(rep(NA,length(subsets)))
  site_funs <- list(rep(NA,length(subsets)))
  
  
  for(i in subsets){ # for individual cluster models
    regressors_sel <- regressors[-which(regressors == 'site_no')]
    lm_data_lm <- lm_data[ssc_subset == i] # only chooses sites within cluster
    # lm_data_hold <- lm_data_lm[-which(lm_data_lm$site_no %in% holdout_sts),] # for eliminating certain sites
    lm_data_hold <- lm_data_lm[which(lm_data_lm$holdout25 == 'in'),]
    glm_y <- as.matrix(lm_data_hold[,log10_SSC_mgL])
    glm_x <- as.matrix(lm_data_hold[,..regressors_sel])
    
    ssc_lm <- cv.glmnet(x = glm_x, y = glm_y, family = 'gaussian', type.measure = "mse", nfolds = 10)
    cv.opt <- coef(ssc_lm, s = "lambda.1se")
    coef_ex <- cbind(rownames(cv.opt),as.numeric(cv.opt))
    colnames(coef_ex) <- c('variable', 'value')
    
    write.table(coef_ex, sep = ",", file = paste0(wd_exports,'cluster', i,'_lasso_fit_coeff.csv'), row.names = F)
    glm_x <- NA
    glm_y <- NA
    
    cluster_funs[[i]] <- ssc_lm
    glm_pred <- predict(ssc_lm, newx = as.matrix(lm_data_lm[,..regressors_sel]), s = "lambda.1se")
    lm_data$pred_cl[which(lm_data$ssc_subset == i)] <- glm_pred
    
    # lm_data$res_cl[which(lm_data$ssc_subset == i)] <- resid(ssc_lm)
  }
  
  return(list(lm_data, cluster_funs))
}

  # Calculate model relative error and station bias (following Morley, 2018)
  # Version that generates plots
  getErrorBias <- function(dt, subset_name){
  # test
  # dt <- ls_sample_models[[1]]
  # subset_name <- 'ncl5_sample_500perstn'
  # Error computed following Morley et. al, 2018
  # Measures of Model Performance Based On the Log Accuracy Ratio
  
  rel_error <- dt[holdout25 == 'holdout',.(
    mape_gl_ind = (10^median(abs(log10(10^pred_gl/10^log10_SSC_mgL)), na.rm = T)-1),
    mape_cl_ind = (10^median(abs(log10(10^pred_cl/10^log10_SSC_mgL)), na.rm = T)-1),
    mape_st_ind = (10^median(abs(log10(10^pred_st/10^log10_SSC_mgL)), na.rm = T)-1)
  )]
  
  # For plotting: get percentage breaks for vertical dashed lines
  percent_error_dt <- data.table(error_breaks = c(0.1,0.5,1,2), 
                                 error_labels = c('< 10 %', '< 50 %', 
                                                  '< 100 %', '< 200 %'))
  
  # Plot histogram of individual errors
  error_histogram <- ggarrange(
    # ggplot(ls_sample_models[[1]][holdout25 == 'holdout'],
    ggplot(dt[holdout25 == 'holdout'],
           aes(y = (10^abs(log10(10^pred_cl/10^log10_SSC_mgL)))-1)) + 
      geom_boxplot() + 
      scale_y_log10(limits = c(0.001, 200), labels = fancy_scientific) + 
      season_facet +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank()
      ) +
      rotate(), 
    ggplot(dt[holdout25 == 'holdout'], 
           aes(x = (10^abs(log10(10^pred_cl/10^log10_SSC_mgL)))-1)) +
      geom_histogram(binwidth = 0.25) + 
      geom_vline(xintercept = c(0.1,0.5,1,2), color = 'navy', lty = 'dashed') +
      geom_text(data = percent_error_dt, aes(x = error_breaks, y = nrow(dt)/20,
                                             label = error_labels),
                color = 'navy', angle = 90, vjust = -0.2, hjust = 0, size = 2.5) +
      scale_x_log10(limits = c(0.001, 200), labels = fancy_scientific) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.2))) +
      season_facet +
      labs(x = 'Relative error',
           y = 'Count'),
    nrow = 2,
    align = 'v', heights = c(0.25,2)) 
  
  ggsave(error_histogram, filename = paste0(wd_exports,'rel_error_hist_', subset_name, '.pdf'), width = 4, height = 4, useDingbats = F)
  
  # Summarize data at each station
  station_summary_subset <- dt[holdout25 == 'holdout'][,.(
    log10_SSC_mgL = mean(log10_SSC_mgL, na.rm = T),
    N_samples = .N,
    p63 = mean(p63, na.rm = T),
    num_pix = mean(num_pix, na.rm = T),
    pred_gl = mean(pred_gl, na.rm = T),
    pred_cl = mean(pred_cl, na.rm = T),
    pred_st = mean(pred_st, na.rm = T)),
    keyby = .(agency_cd, site_no, cluster_sel, Latitude, Longitude, 
              drainage_area_km2, begin_date, end_date)
    ]
  
  # print(station_summary_subset)
  # Compute relative bias at every station; calculate median of all stations for global, cluster, station models
  # Following Morley et al., 2018
  stn_rel_bias <- station_summary_subset[N_samples > 10][,':='(
    mape_gl_stn = (10^abs(log10(10^pred_gl/10^log10_SSC_mgL))-1),
    mape_cl_stn = (10^abs(log10(10^pred_cl/10^log10_SSC_mgL))-1),
    mape_st_stn = (10^abs(log10(10^pred_st/10^log10_SSC_mgL))-1),
    mape_gl_sign = ifelse(pred_gl > log10_SSC_mgL, 1, -1),
    mape_cl_sign = ifelse(pred_cl > log10_SSC_mgL, 1, -1),
    mape_st_sign = ifelse(pred_st > log10_SSC_mgL, 1, -1))
    ][,':='(
      mape_gl_sign = mape_gl_sign*mape_gl_stn,
      mape_cl_sign = mape_cl_sign*mape_cl_stn,
      mape_st_sign = mape_st_sign*mape_st_stn
    )]
  
  # print(stn_rel_bias)
  # Compute relative bias at every station
  median_rel_bias <- stn_rel_bias[,.(
    bias_gl = median(mape_gl_stn, na.rm = T),
    bias_cl = median(mape_cl_stn, na.rm = T),
    bias_st = median(mape_st_stn, na.rm = T)
  )
  ] 
  median_rel_bias_melt <- melt(median_rel_bias, measure.vars = c('bias_gl','bias_cl','bias_st'))[,
                                                                                                 model := factor(variable, levels = c('bias_gl','bias_cl','bias_st'),
                                                                                                                 labels = c('Base','Cluster','Station'), ordered = T)]
  
  station_bias_melt <- melt(stn_rel_bias, id.vars = c('agency_cd', 'site_no', 'cluster_sel', 'Latitude', 'Longitude', 
                                                      'drainage_area_km2', 'begin_date', 'end_date'),
                            measure.vars = c('mape_gl_sign','mape_cl_sign','mape_st_sign'))[,
                                                                                            model := factor(variable, levels = c('mape_gl_sign','mape_cl_sign','mape_st_sign'),
                                                                                                            labels = c('Base','Cluster','Station'), ordered = T)]
  # print(station_bias_melt)
  # station_bias_cluster_plot <- ggplot(station_bias %>% melt(measure.vars = c('bias_gl','bias_cl','bias_st'))) +
  station_bias_cluster_plot <- ggplot(station_bias_melt) +
    geom_boxplot(aes(x = as.factor(cluster_sel), y = value, 
                     fill = as.factor(cluster_sel)), outlier.shape = NA) + 
    geom_text(data = median_rel_bias_melt,
              aes(x = as.factor(1), y = 2.65, label = paste0('Median bias = ',round(value, 2))),
              hjust = 0,vjust = 0, size = 2.5) +
    scale_fill_brewer(palette = 'Paired') +
    scale_y_continuous(lim = c(-2.75,2.75)) +
    facet_wrap(.~as.factor(model)) +
    season_facet +
    labs(
      x = 'River grouping',
      y = 'Relative bias (at-a-station)'
    )
  # SAVE FIGURE
  ggsave(station_bias_cluster_plot, filename = paste0(wd_exports,'station_bias_', subset_name, '.pdf'), 
         useDingbats = F, width = 5, height = 4)
  
  return(list(cbind(rel_error,median_rel_bias), station_bias_cluster_plot))
}
  # No plots generated
  getErrorBias_simple <- function(dt){
  # test
  # dt <- ls_sample_models[[1]]
  # subset_name <- 'ncl5_sample_500perstn'
  # Error computed following Morley et. al, 2018
  # Measures of Model Performance Based On the Log Accuracy Ratio
  
  rel_error <- dt[holdout25 == 'holdout',.(
    mape_gl_ind = (10^median(abs(log10(10^pred_gl/10^log10_SSC_mgL)), na.rm = T)-1),
    mape_cl_ind = (10^median(abs(log10(10^pred_cl/10^log10_SSC_mgL)), na.rm = T)-1),
    mape_st_ind = (10^median(abs(log10(10^pred_st/10^log10_SSC_mgL)), na.rm = T)-1)
  )]
  
  # For plotting: get percentage breaks for vertical dashed lines
  percent_error_dt <- data.table(error_breaks = c(0.1,0.5,1,2), 
                                 error_labels = c('< 10 %', '< 50 %', 
                                                  '< 100 %', '< 200 %'))
  
  
  # Summarize data at each station
  station_summary_subset <- dt[holdout25 == 'holdout'][,.(
    log10_SSC_mgL = mean(log10_SSC_mgL, na.rm = T),
    N_samples = .N,
    p63 = mean(p63, na.rm = T),
    num_pix = mean(num_pix, na.rm = T),
    pred_gl = mean(pred_gl, na.rm = T),
    pred_cl = mean(pred_cl, na.rm = T),
    pred_st = mean(pred_st, na.rm = T)),
    keyby = .(agency_cd, site_no, cluster_sel, Latitude, Longitude, 
              drainage_area_km2, begin_date, end_date)
    ]
  
  # print(station_summary_subset)
  # Compute relative bias at every station; calculate median of all stations for global, cluster, station models
  # Following Morley et al., 2018
  stn_rel_bias <- station_summary_subset[N_samples > 10][,':='(
    mape_gl_stn = (10^abs(log10(10^pred_gl/10^log10_SSC_mgL))-1),
    mape_cl_stn = (10^abs(log10(10^pred_cl/10^log10_SSC_mgL))-1),
    mape_st_stn = (10^abs(log10(10^pred_st/10^log10_SSC_mgL))-1),
    mape_gl_sign = ifelse(pred_gl > log10_SSC_mgL, 1, -1),
    mape_cl_sign = ifelse(pred_cl > log10_SSC_mgL, 1, -1),
    mape_st_sign = ifelse(pred_st > log10_SSC_mgL, 1, -1))
    ][,':='(
      mape_gl_sign = mape_gl_sign*mape_gl_stn,
      mape_cl_sign = mape_cl_sign*mape_cl_stn,
      mape_st_sign = mape_st_sign*mape_st_stn
    )]
  
  # print(stn_rel_bias)
  # Compute relative bias at every station
  median_rel_bias <- stn_rel_bias[,.(
    bias_gl = median(mape_gl_stn, na.rm = T),
    bias_cl = median(mape_cl_stn, na.rm = T),
    bias_st = median(mape_st_stn, na.rm = T)
  )
  ]
  
  return(cbind(rel_error,median_rel_bias))
}

  # Display log axes labels nicely
  fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- log10(l)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(parse(text = paste("10^",as.character(l),sep = "")))
} 

  # Apply SSC calibration models to make predictions based on new surface reflectance inputs (cluster needed)
  # For base, cluster, and site predictions
  getSSC_pred <- function(lm_data, regressors, cluster_funs){ # Version that includes site specification
  # lm_data <- ls_insitu_site_all
  # cluster_funs <- indiv_site_model[[2]]
  # lm_data <- na.omit(ls_insitu_gc, cols = c(regressors_all, 'cluster_sel', 'log10_SSC_mgL'))
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
  # Just for site-specific predictions
  getSSC_pred_st <- function(lm_data, regressors, cluster_funs, levels_all){ # Version that includes site specification
  # lm_data <- ls_insitu_cl_gc[
  #   sampling_type %chin% c('Automated','Manual') # to subset for only sites in GCMRC
  #   ]
  # cluster_funs <- gc_models[[3]]
  # levels_all = site_levels_sel
  # regressors <- regressors_all
  
  lm_data$pred_st <- NA
  lm_data$ssc_subset <- lm_data$cluster_sel # clusters
  subsets <- unique(lm_data$ssc_subset)
  for(i in subsets){ # for individual cluster models
    print(i)
    levels_sel <- sort(unique(levels_all[cluster_sel == i]$site_no)) # NEEDED for specifying levels
    regressors_sel <- regressors
    lm_data_lm <- subset(lm_data, ssc_subset == i) # only chooses sites within cluster
    # rownames(coef(ssc_lm, s = "lambda.1se"))[c(-1,-2)]
    lm_data_lm$site_no <- factor(lm_data_lm$site_no, levels = levels_sel) # NEEDED for specifying levels
    ssc_lm <- cluster_funs[[i]]
    glm_x_new <- model.matrix( ~ ., lm_data_lm[,..regressors_sel])
    glm_pred <- predict(ssc_lm, newx = glm_x_new, s = "lambda.1se")
    lm_data$pred_st[which(lm_data$ssc_subset == i)] <- glm_pred
    # lm_data$res_cl[which(lm_data$ssc_subset == i)] <- resid(ssc_lm)
    lm_data_lm <- NA
  }
  return(lm_data)
}

  # Plot results of SSC calibration, with several adjustable parameters
  get_sscPlot <- function(ssc_data,ssc_title,density_yn, validation){
  ## Test data
  # ssc_data <- ssc_model_cl_iterate_pred
  # ssc_title <- "byCluster"
  # density_yn <- 'no'
  # validation <- 'yes'
  ## Start function here
  if(ssc_title == "byGlobal"){
    # ssc_data <- ssc_data[,pred := pred_gl]
    cl_colors <- rep('#1f78b4',5)
    ssc_plot <- ggplot(data = ssc_data,   
                       aes(x = SSC_mgL, y = 10^(pred_gl)))
  } else if(ssc_title == "byCluster"){
    ssc_plot <- ggplot(data = ssc_data,   
                       aes(x = SSC_mgL, y = 10^(pred_cl)))
  } else {
    ssc_plot <- ggplot(data = ssc_data,   
                       aes(x = SSC_mgL, y = 10^(pred_st)))
  }
  n_samples <- nrow(ssc_data)
  if(density_yn == "yes"){
    ssc_plot <- ssc_plot + 
      stat_density_2d(aes(fill = stat(nlevel)), 
                      # h = c(1,1),
                      geom = "polygon") +
      # facet_wrap(. ~ as.factor(X5_clusters_weighted)) +
      scale_fill_distiller(palette = 'GnBu', type = 'seq', direction = 1)
  } else if(density_yn == 'overlay'){
    ssc_plot <- ssc_plot +
      geom_point(aes(
                    # fill = as.factor(cluster_sel)
                     # ,color = as.factor(cluster_sel)), pch = 16,
                     ),color = 'black', pch = 16,
                 # color = "black", pch = 21, alpha = 0.8, # comment to remove outline
                 size = 1, alpha = min(1,1/(sqrt(n_samples*0.001)))) + 
      stat_density_2d(aes(fill = stat(nlevel)), 
                      # h = c(1,1),
                      geom = "polygon", alpha = 0.3) +
      # facet_wrap(. ~ as.factor(X5_clusters_weighted)) +
      scale_fill_distiller(palette = 'GnBu', type = 'seq', direction = 1)
      # scale_color_brewer('Paired', type = 'qual') + scale_fill_brewer('Paired', type = 'qual')
      # scale_color_manual(values = cl_colors) + scale_fill_manual(values = cl_colors)
  }else{
    ssc_plot <- ssc_plot +
      geom_point(aes(
        # fill = as.factor(cluster_sel)
        # ,color = as.factor(cluster_sel)), pch = 16,
      ),color = 'black', pch = 16,
      # color = "black", pch = 21, alpha = 0.8, # comment to remove outline
      size = 1, alpha = min(1,1/(sqrt(n_samples*0.001))))
  }
  if(validation == 'yes'){
    ssc_plot <- ssc_plot + facet_wrap(.~factor(holdout25, levels = c('in','holdout'), ordered = T))
  }
  
  ssc_plot <- ssc_plot + 
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, color = 'orange') +
    scale_y_log10(limits = c(1,50000), labels = fancy_scientific, breaks = breaks, minor_breaks = minor_breaks,
                  expand = expand_scale(add = c(0.7,0.3))) + 
    scale_x_log10(limits = c(1,50000), labels = fancy_scientific, breaks = breaks, minor_breaks = minor_breaks, 
                  expand = expand_scale(add = c(0.7,0.3))) +
    theme(
      legend.title = element_blank(),
      legend.position = 'none',
      text = element_text(size=15),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      strip.text = element_blank()
    ) + 
    labs(
      x = expression(paste(italic("in situ"), " SSC (mg/L)")),
      y = "Satellite-estimated SSC (mg/L)"
      # fill = "Cluster", 
      # color = "Cluster"
    )
  return(ssc_plot)
}

  # Modify USGS station number to have required 8 digits
  getUSGS_site_no <- function(x){
    if(!is.na(x)){
      if(nchar(x) == 7 & !is.na(as.numeric(x))){
        return(as.character(paste0('0',x)))}
      else{return(x)}}else{return(NA)}}

#### THEMES AND PLOTTING PARAMETERS ####

breaks <- 10^(-10:10)
minor_breaks <- rep(5, 21)*(10^rep(-10:10, each=9))

theme_evan <- theme_bw() +
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

theme_evan_facet <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    # panel.grid = element_blank(),
    # legend.position = 'none',
    panel.border = element_rect(size = 0.5),
    strip.background = element_rect(fill = 'white'),
    text = element_text(size=8),
    axis.text = element_text(size = 8), 
    plot.title = element_text(size = 9)
  )

season_facet <- theme_evan_facet + theme(
  legend.position = 'none', 
  strip.background = element_blank(),
  strip.text = element_text(hjust = 0, margin = margin(0,0,0,0, unit = 'pt'))
)




