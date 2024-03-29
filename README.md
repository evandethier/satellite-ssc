# satellite-ssc #
## Estimating suspended sediment concentration (SSC) using a globally applicable calibration for surface reflectance satellite images ##

Calibration has been completed for Landsat 5 and 7 LT1 Surface Reflectance Product (USGS/NASA)

The calibrations described here were made using 134,697 *in situ* measurements at 727 stations in the United States (United States Geological Survey, USGS), Canada (Water Survey of Canada, WSC), South America (HYdro-geochemistry of the AMazonian Basin, HYBAM), and Taiwan (Water Resources Agency, WRA). 

For information on how to implement these techniques for calibration or prediction, refer to the documentation provided in the [**satellite-ssc**](landsat-calibration) subfolders.

For our application of this work to global rivers, [published in Science](https://doi.org/10.1126/science.abn7980), refer to [**outlet-rivers**](outlet-rivers) subfolders.

# Introduction to calibration techniques and results #

**Calibration stations and satellite imagery**

Because the river water and the sediment/organic matter it transports has variable characteristics that can affect its reflectance of light, we use unsupervised K-Means clustering to categorize rivers as one of six groups. We then made a category-specific calibration for each river grouping. The basics of this approach are described below.

<p align="center">
  <img src="/Readme_figures/fig1_cluster_map_n6.jpg" width="85%" >
</p>

The Landsat 5 satellite operated from 1984-2012, and the Landsat 7 satellites has been in operation since 1999. Each satellite revists every location on earth every 16 days, though there is some overlap in the paths they take, meaning some areas are revisited twice every 16 days. Clouds obscure the earth from the view of satellites, so not every image is useful for calibration development or SSC estimation. However, on average the stations in the calibration dataset had 672 usable images of their river location.

<p align="center">
  <img src="/Readme_figures/fig2_n_sat_samples_histogram.jpg" width="60%" >
</p>

We matched these satellite images with *in situ* measurements from the same day or, if no measurements were taken that day, from as much as +/- 2 days from the Landsat image date. Image-*in situ* measurement pairs were distributed throughout the Landsat record, with the most images during the time when both Landsat satellites were operational (1999-2012).

<p align="center">
  <img src="/Readme_figures/fig3_sampling_date_histogram.jpg" width="60%" >
</p>

We classified pixels in each Landsat image as water or not water (using the USGS/NASA LT1 Surface Reflectance product). Then, for each image of a given *in situ* station, we sampled the surface reflectance in the visible, near-infrared, and shortwave infrared wavelengths within 200 m of the station. We then used these image-*in situ* measurement pairs to develop calibrations for estimating SSC based on Landsat imagery.

<p align="center">
  <img src="/Readme_figures/figs1-workflow-yellowstone-sidney-example.jpg" width="70%" >
</p>


**Clustering of optically similar rivers improves calibration**

We tested different numbers of cluster groups, selecting a six-group approach by minimizing relative calibration model error and at-a-station bias.  Average true color appearance of river pixels for a gradient of suspended sediment concentrations shows variability based on cluster. This color variability is related to typical percent organic carbon content. 

<p align="center">
  <img src="/Readme_figures/fig5_cluster_combined_fig.jpg" width="70%" >
</p>


**Influence of variations in sediment grain size and organic matter on calibration success**

We found that variations in the grain-size of suspended sediment affected the success of calibration models. When more than 50% of the suspended sediment was of sand size or coarser, the calibration tended to underpredict the SSC. This is because sand is less reflective per mass than finer grained sediment. Fortunately, most rivers large enough to be used for these methods (~90 meters wide) carry sediment without much sand. 

<p align="center">
  <img src="/Readme_figures/fig8_p63_vs_error_ncl5.jpg" width="60%" >
  <img src="/Readme_figures/fig9_usa_grain_size_map.jpg" width="95%" >
</p>

Percent organic carbon (POC) of the suspended sediment also contributed to some degree of uncertainty. When POC made up a high fraction of suspended sediment concentrations, the calibration tended to overpredict SSC. This challenge mostly applies to rivers that don't have high suspended sediment concentration, as percent POC tends to decrease with increased SSC.

<p align="center">
  <img src="/Readme_figures/fig10_POC_fraction_vs_error_ncl5.jpg" width="70%" >
</p>




