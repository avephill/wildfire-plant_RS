# Scripts in study *How Competition and Wildfire Affect Tree Migration in the American West* (in review)
## Abstract:
Due to climate change, plant populations experience environmental conditions to which they are not adapted. Our understanding of the next century’s vegetation geography depends on the distance, direction, and rate at which plants redistribute in response to a changing climate. Although plant migration in response to contemporary climate change is widely observed, our understanding of its mechanics is nascent. In this study we test the response of plant migration rates to wildfire occurrence using 33,838 Forest Inventory Analysis plots across five states in the western United States. Wildfire increases the rate of migration towards suitable climates for seven tree species (by more than 22% on average), suggesting that incumbent vegetation can act as a barrier to plant migration and that fire management may be an important facet of the management of vegetation transitions in response to climate change.

## Using these scripts
### Climate Data
Download CMIP5 climate data for [1961-1990 here](http://www.cacpd.org/climate_normals/NA_NORM_6190_Bioclim_netCDF.7z) and [1981-2010 here](http://www.cacpd.org/climate_normals/NA_NORM_8110_Bioclim_netCDF.7z) and place folders in YourProjectDirectory/Data/Environmental/Original/

### Species Occurrence Data
Download [statewide FIA information](https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html) (this study used CA, OR, ID, WA, MT) and place folders in YourProjectDirectory/Data/Occurrence/Original/FIA_by_state/

### Running the Scripts
Execute both data preparation scripts in YourProjectDirectory/Scripts/0-Data_Prep/–– climate_raster_prep.R then occurrence_data_prep.R

Once complete, then executing YourProjectDirectory/Scripts/0-primary_analyis.R should perform the majority of analyses presente in the study.
