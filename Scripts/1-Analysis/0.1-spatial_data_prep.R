# Here we trim spatial data (occurrence and environmental) to parameters specified in 0-primary_analysis
# This is called from 0-primary_analysis

# political boundaries for visualization
setwd("Data/Environmental")
na_outline <- readOGR('../Occurrence/Original/NA_Outline','NA')
state_mask <- na_outline[na_outline@data$NAME %in% c("California","Oregon","Washington","Montana","Idaho"),]
state_mask <- spTransform(state_mask,CRS(crs_wgs84))
west_us <- state_mask

# Creates ecoregion mask for climate variables and extent of study and
# applies it to plots considered in study
shp_ecoregion <- readOGR("Original/NA_CEC_Eco_Level1/NA_CEC_Eco_Level1.shp")

shp_nfm <- shp_ecoregion[shp_ecoregion@data$NA_L1NAME=="NORTHWESTERN FORESTED MOUNTAINS",]
shp_mwcf <- shp_ecoregion[shp_ecoregion@data$NA_L1NAME=="MARINE WEST COAST FOREST",]
shp_eco <- raster::aggregate(rbind(shp_nfm,shp_mwcf))
ecomask <- spTransform(shp_eco,crs_wgs84)
ecomask_coords <- data.frame(Longitude=prep_plts$LON,Latitude=prep_plts$LAT,ID=prep_plts$CN)
ecomask_coords <- ecomask_coords[complete.cases(ecomask_coords),]
coordinates(ecomask_coords) <- ~ Longitude + Latitude
proj4string(ecomask_coords) <- proj4string(ecomask)
plmasked <- ecomask_coords[ecomask,] 
masked_PLT_CN <- plmasked$ID
plts <- prep_plts[prep_plts$CN %in% masked_PLT_CN,]
plt_difs <- prep_plt_difs[prep_plt_difs$PLT_CN %in% masked_PLT_CN,]
plt_difs_JG <- prep_plt_difs_JG[prep_plt_difs_JG$PLT_CN %in% masked_PLT_CN,]


if(bool_plot_climate_change==T){
  files<-list.files(dir_cc)
  cc_dirs<-paste(dir_cc,files,sep="/")
}
files <- list.files(dir_climate)
dirs <- paste(dir_climate,files,sep="/")

# Isolate parameter names
parnames <- sapply(files, function(x) strsplit(x, "[.]")[[1]], USE.NAMES=FALSE)[1,]

# Converts parameter acronyms to full name
bioclim_dict <- read.csv("bioclim_dict.csv",stringsAsFactors = F)
 
process.env.raster <- function(layer_directory){
  parameter <- raster(layer_directory[1])
  parameter <- mask(parameter,ecomask)
  
  pts_raster <- rasterToPoints(parameter, spatial = TRUE)
  pts_lon_lat <- spTransform(pts_raster, CRS(crs_wgs84))
  pts_df <- as.data.frame(pts_lon_lat)
  
  # These are the lonlat coordinates, determined from the resolution 
  # of the first layer, that will be used to sample all other rasters
  env.lonlat <- pts_df[,2:3] 
  
  col_nam <- as.character(bioclim_dict[bioclim_dict$abr==names(parameter),"fullname"])
  env.lonlat[,col_nam]<-pts_df[,1]
  resolution_raster <- res(parameter)[1] #resolution of parameter, used in calculation of Schoener's D
  
  # Loop through other environmental layers and process them similarly
  for(j in 2:length(layer_directory)){
    parameter<-raster(layer_directory[j])
    x<-raster::extract(parameter,env.lonlat[,1:2])
    col_nam<-as.character(bioclim_dict[bioclim_dict$abr==names(parameter),"fullname"])
    env.lonlat[,col_nam]<-x
  }
  
  env.lonlat <- env.lonlat[complete.cases(env.lonlat),] #gets rid of NAs if any (found some in aevpet)
  #change all the names
  names(env.lonlat) <- sapply(names(env.lonlat), function(x){
    if(x %in% bioclim_dict$abr){
      bioclim_dict[bioclim_dict$abr == x,"fullname"]
    } else{
      return(x)
    }
  })
  return(list(env.lonlat=env.lonlat,resolution_raster=resolution_raster))
}

raster_results_list <- process.env.raster(layer_directory = dirs)
env.lonlat <- raster_results_list$env.lonlat
resolution_raster <- raster_results_list$resolution_raster


# This is for scaling results (climate niches, etc) after primary analysis. 
scaling_vector <- apply(env.lonlat[,!names(env.lonlat) %in% c('x','y')], 2, sd)
scaling_vector <- t(data.frame(scaling_vector))
write.csv(scaling_vector,'Processed/scaling_vector.csv',row.names=F)

layer_names <- names(env.lonlat[,3:length(env.lonlat)])

# If we're considering the change in climate over the century, make climate space for change in climate:
if(bool_plot_climate_change==T){
  cc_env.lonlat <- process.env.raster(layer_directory = cc_dirs)$env.lonlat
  
  cc_scaling_vector <- apply(cc_env.lonlat[,!names(cc_env.lonlat) %in% c('x','y')], 2, sd)
  cc_scaling_vector <- t(data.frame(cc_scaling_vector))
  write.csv(cc_scaling_vector,'Processed/scaling_vector.csv',row.names=F)
}
