# Here we trim spatial data (occurrence and environmental) to parameters specified in 0-primary_analysis
# This is called from 0-primary_analysis

source("Scripts/1-Analysis/0.1.1-Load_Spatial_Data.R", echo = T)
setwd("Data/Environmental/")



files <- list.files(dir_climate)
dirs <- paste(dir_climate, files, sep="/")

# Isolate parameter names
parnames <- sapply(files, function(x) strsplit(x, "[.]")[[1]], USE.NAMES=FALSE)[1,]


raster_results_list <- process.env.raster(layer_directory = dirs)
env.lonlat <- raster_results_list$env.lonlat %>% data.frame()
env.stack <- raster_results_list$env.stack
resolution_raster <- raster_results_list$resolution_raster


# This is for scaling results (climate niches, etc) after primary analysis. 
scaling_vector <- apply(env.lonlat[,!names(env.lonlat) %in% c('x','y')], 2, sd)
scaling_vector <- t(data.frame(scaling_vector))
write.csv(scaling_vector,'Processed/scaling_vector.csv', row.names = F)

layer_names <- names(env.lonlat[,3:length(env.lonlat)])

setwd(project_directory)
