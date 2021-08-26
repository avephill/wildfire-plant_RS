# This script takes environmental raster data from CMIP5 and produces rasters
# in the format, size, and locality that we use in this study
# Download 
# http://www.cacpd.org/climate_normals/NA_NORM_6190_Bioclim_netCDF.7z
# and
# http://www.cacpd.org/climate_normals/NA_NORM_8110_Bioclim_netCDF.7z
# and place folders in projectdirectory/Data/Environmental/Original/

library(raster)
if(!exists("project_directory")) project_directory <- "" # Enter location of project directory within quotation marks
setwd(project_directory)
na_outline <- readOGR('Data/Occurrence/Original/NA_Outline','NA')
state_mask <- na_outline[na_outline@data$NAME %in% c("California","Oregon",
                                                     "Washington","Montana",
                                                     "Idaho", "Wyoming",
                                                     "Utah", "New Mexico",
                                                     "Colorado"),]
state_mask <- spTransform(state_mask, CRS(crs_wgs84))
crs_lcc <- "+proj=lcc +lat_1=49.0 +lat_2=77.0 +lat_0=0.0 +lon_0=-95.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m +no_defs"
crs_wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
dir.create("Data/Environmental/Processed/CMIP6190")
dir.create("Data/Environmental/Processed/CMIP8110")
dir.create("Data/Environmental/Processed/CMIPdif")
ncPar <- c("PPT_wt",'MWMT','PPT_sm', 'MCMT', "Eref")

for(i in 1:length(ncPar)){
  
  print(i/length(ncPar))
  par <- ncPar[i]
  par.nc<-paste(par,"nc",sep=".")
  print(par.nc)
  
  if(par=="aevpet.nc"){
    par8110 <- raster(paste("Data/Environmental/Original/NA_NORM_8110_Bioclim_netCDF",par.nc,sep="/"))
    par8110 <- projectRaster(par8110, crs = crs_wgs84) %>% crop(extent(state_mask))
    writeRaster(par8110,paste("Data/Environmental/Processed/CMIP8110",par.nc,sep="/"),overwrite=TRUE)
  } else{
    
    par8110<- raster(paste("Data/Environmental/Original/NA_NORM_8110_Bioclim_netCDF",par.nc,sep="/"))
    crs(par8110)<-crs_lcc
    par8110<-projectRaster(par8110, crs=crs_wgs84) %>% crop(extent(state_mask))
    writeRaster(par8110,paste("Data/Environmental/Processed/CMIP8110",par.nc,sep="/"),overwrite=TRUE)
    
    par6190<- raster(paste("Data/Environmental/Original/NA_NORM_6190_Bioclim_netCDF",par.nc,sep="/"))
    crs(par6190)<-crs_lcc
    par6190<-projectRaster(par6190 ,crs=crs_wgs84) %>% crop(extent(state_mask))
    writeRaster(par6190,paste("Data/Environmental/Processed/CMIP6190",par.nc,sep="/"),overwrite=TRUE)
    
    #sp::plot(par6190)
    parDif <- par8110 - par6190
    names(parDif) <- paste("d_",par,sep="")
    writeRaster(parDif,paste("Data/Environmental/Processed/CMIPdif/","d_",par.nc,sep=""), overwrite=TRUE)
  }
}
