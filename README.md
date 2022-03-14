# Scripts in study [*Forest fires and climate-induced tree range shifts in the western US*](https://doi.org/10.1038/s41467-021-26838-z)

## Abstract:
Due to climate change, plant populations experience environmental conditions to which they are not adapted. Our understanding of the next century’s vegetation geography depends on the distance, direction, and rate at which plant distributions shift in response to a changing climate. In this study we test the sensitivity of tree range shifts (measured as the difference between seedling and mature tree ranges in climate space) to wildfire occurrence, using 74,069 Forest Inventory Analysis plots across nine states in the western United States. Wildfire significantly increased the seedling-only range displacement for 2 of the 8 tree species in which seedling-only plots were displaced from tree-plus-seedling plots in the same direction with and without recent fire. The direction of climatic displacement was consistent with that expected for warmer and drier conditions. The greater seedling-only range displacement observed across burned plots suggests that fire can accelerate climate-related range shifts and that fire and fire management will play a role in the rate of vegetation redistribution in response to climate change.

## Using these scripts
### Climate Data
Download CMIP5 climate data for [1961-1990 here](http://www.cacpd.org/climate_normals/NA_NORM_6190_Bioclim_netCDF.7z) and [1981-2010 here](http://www.cacpd.org/climate_normals/NA_NORM_8110_Bioclim_netCDF.7z) and place folders in YourProjectDirectory/Data/Environmental/Original/

### Species Occurrence Data
Download [statewide FIA information](https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html) (this study used CA, OR, ID, WA, MT) and place folders in YourProjectDirectory/Data/Occurrence/Original/FIA_by_state/

### Running the Scripts
Execute both data preparation scripts in YourProjectDirectory/Scripts/0-Data_Prep/–– climate_raster_prep.R then occurrence_data_prep.R

Once complete, then executing YourProjectDirectory/Scripts/0-primary_analyis.R should perform the majority of analyses presente in the study.
