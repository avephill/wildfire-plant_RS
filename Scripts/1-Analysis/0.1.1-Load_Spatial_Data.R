# This is for loading in some spatial data that is needed for both primary analysis and 
# data presentation
# Dataframe with all plots and fire information


## Load in occurrence Data
setwd("Data/Occurrence/Processed")
# This dataframe indicates the species compisition of the follow categories for each <- plot:
# Only Adult Trees, Only Saplings, or Both
prep_plt_difs<-read.csv("plot_differences.csv")

prep_plts <- read.csv("pl_fire.csv", stringsAsFactors = F) 

# This is a dataset that indicates what plots have, of each species:
# Only big mature Trees, Only saplings, or Both
prep_plt_difs_JG <- read.csv("plot_differences_sapling_mature.csv", stringsAsFactors = F) 

# We don't want fire in the sapling/mature analysis
prep_plt_difs_JG <- subset(prep_plt_difs_JG, is.na(prep_plt_difs_JG$FIRE)) 

spcd <- read.csv("../Original/species_codes.csv") #this translates species ID's to species names
tr_lite <- read.csv("tr_lite.csv", stringsAsFactors = F)
sd_lite <- read.csv("sd_lite.csv", stringsAsFactors = F)


## Load in Environmental Data
crs_wgs84. <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
setwd(project_directory)
setwd("Data/Environmental")
# This is the directory in which I'll find all 
# the environmental variables used in the study
dir_climate <- "Processed/CMIP8110"



# Here we trim spatial data (occurrence and environmental) to parameters specified
# This R script uses crs_wgs84, dir_climate, dir_cc, prep_plts, prep_plt_difs_JG, prep_plt_difs
# produces plts, plt_difs, plt_difs_JG, env.lonlat, cc_env.lonlat, scale_vector
setwd(project_directory)

crs_wgs84<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"


renameClimColumns <- function(df){
  names(df) <- sapply(names(df), function(x){
    if(x %in% bioclim_dict$abr){
      bioclim_dict[bioclim_dict$abr == x, "fullname"]
    } else{
      return(x)
    }
  })
  return(df)
}

process.env.raster <- function(layer_directory){
  
  study_mask <- raster::intersect(state_mask, ecomask)
  
  parameters <- raster::stack(layer_directory) %>% 
    raster::mask(study_mask)
  
  env.lonlat <- rasterToPoints(parameters) %>% na.omit
  resolution_raster <- res(parameters)[1] #resolution of parameter, used in calculation of Schoener's D
  
  #change all the names
  # env.lonlat <- renameClimColumns(env.lonlat)
  
  return(list(env.lonlat = env.lonlat, resolution_raster = resolution_raster,
              env.stack = parameters))
}


na_outline <- readOGR('Data/Occurrence/Original/NA_Outline','NA')
state_mask <- na_outline[na_outline@data$NAME %in% c("California","Oregon",
                                                     "Washington",
                                                     "Montana","Idaho",
                                                     "Wyoming","Colorado", 
                                                     "New Mexico", "Utah"),] %>% 
  spTransform(CRS(crs_wgs84))
west_us <- state_mask

# Creates ecoregion mask for climate variables and extent of study and
# applies it to plots considered in study
shp_ecoregion <- readOGR("Data/Environmental/Original/NA_CEC_Eco_Level1/")
shp_nfm <- shp_ecoregion[shp_ecoregion@data$NA_L1NAME=="NORTHWESTERN FORESTED MOUNTAINS",]
shp_mwcf <- shp_ecoregion[shp_ecoregion@data$NA_L1NAME=="MARINE WEST COAST FOREST",]
shp_eco <- raster::aggregate(rbind(shp_nfm, shp_mwcf))
ecomask <- spTransform(shp_eco, crs_wgs84)

ecomask_coords <- data.frame(Longitude=prep_plts$LON, Latitude=prep_plts$LAT, ID=prep_plts$CN)
ecomask_coords <- ecomask_coords[complete.cases(ecomask_coords),]
coordinates(ecomask_coords) <- ~ Longitude + Latitude
proj4string(ecomask_coords) <- proj4string(ecomask)
plmasked <- ecomask_coords[ecomask,] 


masked_PLT_CN <- plmasked$ID
plts <- prep_plts[prep_plts$CN %in% masked_PLT_CN,]

plt_difs_JG <- prep_plt_difs_JG[prep_plt_difs_JG$PLT_CN %in% masked_PLT_CN,]

plt_difs <- prep_plt_difs[prep_plt_difs$PLT_CN %in% masked_PLT_CN,]



