### Quantify consistency of migration direction through climate space 
### between burned/unburned sites and between different life stages

# To compare species in different life stages (off- spring vs. adult), we followed the FIA definition, 
# dividing the data into two types of subgroups, (i) seedling (dbh < 2.5 cm) vs. tree (dbh > 2.5 cm),
# and (ii) sapling (2.5 cm < dbh < 12.7 cm) vs. large tree (dbh > 12.7 cm).

# Between burned and unburned sites:
boot.resamples <- 5000


burnagreement_bootstrap.df <- rbind(clim_scaled_BT %>% mutate(ID = "BT"),
                      clim_scaled_BS %>% mutate(ID = "BS"),
                      clim_scaled_NT %>% mutate(ID = "NT"),
                      clim_scaled_NS %>% mutate(ID = "NS")) %>% 
  mutate(ID = factor(ID))

burn.agreement <- function(data, indices){
  data_samp <- data[indices,]
  
  clim_centroid_scaled_BT <- colMeans(data_samp %>% dplyr::filter(ID == "BT") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  clim_centroid_scaled_BS <- colMeans(data_samp %>% dplyr::filter(ID == "BS") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  clim_centroid_scaled_NT <- colMeans(data_samp %>% dplyr::filter(ID == "NT") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  clim_centroid_scaled_NS <- colMeans(data_samp %>% dplyr::filter(ID == "NS") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  
  dist_BN <- dist(rbind(clim_dif_N, clim_dif_B), method="euclidean") # Distance between vector termini points
  dist_B  <- dist(rbind(clim_centroid_scaled_BS, clim_centroid_scaled_BT), method="euclidean")# This is the magnitude of Fire Vector
  dist_N  <- dist(rbind(clim_centroid_scaled_NS, clim_centroid_scaled_NT), method="euclidean")# This is magnitude of Non-Fire vector
  
  # arccos(A/H) = angle
  # use the law of cosines to find the angle A between the vectors: 
  # dist_BN^2 = dist_B^2 + dist_N^2 - 2*dist_B*dist_N*cos(A)
  A <- as.numeric(acos((dist_B^2 + dist_N^2 - dist_BN^2) / (2*dist_B*dist_N)))
  # print(A)
  if(is.na(A)){
    fire_agreement <- NA
    
  } else if(A > pi/2){
    fire_agreement <- 0
    
  }  else {
    # this here calculates the proportion of fire vector that lies along outside-fire vector, inspired by http://www.nabla.hr/VT-Vectors2Dand3DB5.htm
    # dist_B*cos(A) is the magnitude of dist_B that lies along dist_N
    # then we divide it by the entirety of dist_B and multiply it by 100
    # fire_agreement <- (dist_B*cos(A) / dist_B) * 100
    fire_agreement <- cos(A) * 100
  }
  
  return(fire_agreement)
}

burnagreement_boot <- boot(burnagreement_bootstrap.df, burn.agreement, R = boot.resamples,
                           stype = "i",
                           parallel = "multicore",
                           strata = burnagreement_bootstrap.df$ID)

# gets the standard error of the mean
burnagreement <- burnagreement_boot$t0
burnagreement_se <- summary(burnagreement_boot)$bootSE




# Between Different Life Stages:

lonlat_J <- plts[plts$CN %in% plt_difs_JG[plt_difs_JG$NVL_SEED == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_J)<-c('x','y')

lonlat_G <- plts[plts$CN %in% plt_difs_JG[plt_difs_JG$BOTH == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_G)<-c('x','y')
clim12 <- env.lonlat

occ.sp <- rbind(lonlat_G,lonlat_J)

occ.G <- raster::extract(env.stack, lonlat_G, df = T)  %>% 
  dplyr::select(-ID) %>%  mutate(x = lonlat_G$x, y = lonlat_G$y) %>% na.omit 
occ.G <- occ.G %>% 
  dplyr::select(x, y, names(occ.G %>% dplyr::select(-x,-y)))

occ.J <- raster::extract(env.stack, lonlat_J, df = T) %>% 
  dplyr::select(-ID) %>%  mutate(x = lonlat_J$x, y = lonlat_J$y) %>% na.omit 
occ.J <- occ.J %>% 
  dplyr::select(x, y, names(occ.J %>% dplyr::select(-x,-y)))

env.lonlat_G <- occ.G
env.lonlat_J <- occ.J
env.lonlat_G$ID <- "Both"
env.lonlat_J$ID <- "Novel"
env.lonlat_JG <- rbind(env.lonlat_G, env.lonlat_J)
env.lonlat_JG <- env.lonlat_JG[,!names(env.lonlat_JG) %in% c('x','y')]
env.lonlat_JG$Ver <- "JG"
env.lonlat_scaled_JG <- data.frame(scale(env.lonlat_JG[,!names(env.lonlat_JG) %in% c("ID","Ver")],
                                         center = F, scale=scaling_vector))
env.lonlat_scaled_JG$ID <- env.lonlat_JG$ID




lifeagreement_bootstrap.df <- rbind(env.lonlat_scaled_JG %>% 
                                      dplyr::filter(ID == "Both") %>% mutate(ID = "NG"),
                                    env.lonlat_scaled_JG %>% 
                                      dplyr::filter(ID == "Novel") %>% mutate(ID = "NJ"),
                                    clim_scaled_NT %>% mutate(ID = "NT"),
                                    clim_scaled_NS %>% mutate(ID = "NS")) %>% 
  mutate(ID = factor(ID))
  

life.agreement <- function(data, indices){
  
  # stratified sample
  data_samp <- data[indices,]
  
  clim_centroid_scaled_NG <- colMeans(data_samp %>% dplyr::filter(ID == "NG") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  clim_centroid_scaled_NJ <- colMeans(data_samp %>% dplyr::filter(ID == "NJ") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  clim_centroid_scaled_NT <- colMeans(data_samp %>% dplyr::filter(ID == "NT") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  clim_centroid_scaled_NS <- colMeans(data_samp %>% dplyr::filter(ID == "NS") %>% 
                                        dplyr::select(-ID) %>% na.omit)
  
  
  clim_dif_JG <- clim_centroid_scaled_NJ - clim_centroid_scaled_NG
  
  clim_dif_ST <- clim_centroid_scaled_NS - clim_centroid_scaled_NT # Calculated in 0-primary_analysis (climatic dif of seedling/tree outside of fire)
  dist_STJG <- dist(rbind(clim_dif_JG, clim_dif_ST),method="euclidean") # distance between vector termini
  dist_ST  <- dist(rbind(clim_centroid_scaled_NS, clim_centroid_scaled_NT), method="euclidean")# This is magnitude of Non-Fire vector
  dist_JG <- dist(rbind(clim_centroid_scaled_NJ, clim_centroid_scaled_NG), method="euclidean")# this is magnitude of sapling, big tree vector
  
  #arccos(A/H)= angle
  #use the law of cosines to find the angle A between the vectors: 
  A <- as.numeric(acos((dist_ST^2 + dist_JG^2 - dist_STJG^2)/(2*dist_ST*dist_JG)))
  if(is.na(A)){
    lifestage_agreement <- NA
    
  }else if(A > pi/2){
    lifestage_agreement <- 0
  } else{
    #this here calculates the proportion of dist_ST that lies along dist_JG, inspired by http://www.nabla.hr/VT-Vectors2Dand3DB5.htm
    # dist_ST*cos(A) is the magnitude of dist_ST that lies along dist_JG
    # then we divide it by the entirety of dist_ST and multiply it by 100
    lifestage_agreement <- (dist_ST*cos(A) / dist_ST) * 100
  }
  
  return(lifestage_agreement)
}

lifeagreement_boot <- boot(lifeagreement_bootstrap.df, life.agreement, 
                           R = boot.resamples, parallel = "multicore",
                           stype = "i",
                           strata = lifeagreement_bootstrap.df$ID)
lifeagreement <- lifeagreement_boot$t0
lifeagreement_se <- summary(lifeagreement_boot)$bootSE


