# landsat-57 #
## Estimating suspended sediment concentration (SSC) using a globally applicable calibration for Landat 5 and 7 surface reflectance satellite images ##

Calibrations have been completed for **Landsat 5 and 7** LT1 Surface Reflectance (SR) Product (USGS/NASA)

The calibrations described here were made using 134,697 *in situ* measurements at 727 stations in the United States (United States Geological Survey, USGS), Canada (Water Survey of Canada, WSC), South America (HYdro-geochemistry of the AMazonian Basin, HYBAM), and Taiwan (Water Resources Agency, WRA).

Based on Landsat surface reflectance data obtained using this Google Earth Engine [Landsat image sampling  protocol](https://github.com/evandethier/earthengine), the code included here is used to 

  -[Download data from online repositories](insitu-data-download.R)
  
  -[Generate calibrations](landsat-57-calibration.R)
  
  -[Apply calibrations to predict SSC for Landsat 5 and 7 SR data](landsat-57-prediction.R)
  
  
### **Standalone calibrations for rivers in the calibration dataset** ###

In addition to the globally applicable calibration methods, standalone calibrations for 138 rivers are available.

These can be found in the [standalone calibration subfolder](landsat-57-standalone-calibrations)

For example, this is a standalone calibration from **WSC station 05GG001**, 

North Saskatchewan River at Prince Albert, Saskatchewan, CA

[53.20344,-105.77214]

<p align="center">
  <img src="/landsat-calibration/ls57-readme-figs/indiv_calib_plot_05GG001.jpg" width="60%" />
</p>
