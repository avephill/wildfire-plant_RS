### Quantify consistency of migration direction through climate space 
### between burned/unburned sites and between different life stages

# To compare species in different life stages (off- spring vs. adult), we followed the FIA definition, 
# dividing the data into two types of subgroups, (i) seedling (dbh < 2.5 cm) vs. tree (dbh > 2.5 cm),
# and (ii) sapling (2.5 cm < dbh < 12.7 cm) vs. large tree (dbh > 12.7 cm).

# Between burned and unburned sites:

clim_centroid_scaled_BT <-  colMeans(clim_scaled_BT)
clim_centroid_scaled_BS <- colMeans(clim_scaled_BS)
clim_centroid_scaled_NT <- colMeans(clim_scaled_NT)
clim_centroid_scaled_NS <- colMeans(clim_scaled_NS)

dist_BN <- dist(rbind(clim_dif_N, clim_dif_B), method="euclidean") # Distance between vector termini points
dist_B  <- dist(rbind(clim_centroid_scaled_BS, clim_centroid_scaled_BT), method="euclidean")# This is the magnitude of Fire Vector
dist_N  <- dist(rbind(clim_centroid_scaled_NS, clim_centroid_scaled_NT), method="euclidean")# This is magnitude of Non-Fire vector

#arccos(A/H) = angle
#use the law of cosines to find the angle A between the vectors: 
#dist_BN^2 = dist_B^2 + dist_N^2 - 2*dist_B*dist_N*cos(A)
A <- as.numeric(acos((dist_B^2 + dist_N^2 - dist_BN^2)/(2*dist_B*dist_N)))

if(A > pi/2){
  fire_agreement <- 0
} else{
  #this here calculates the proportion of fire vector that lies along outside-fire vector, inspired by http://www.nabla.hr/VT-Vectors2Dand3DB5.htm
  # dist_B*cos(A) is the magnitude of dist_B that lies along dist_N
  # then we divide it by the entirety of dist_B and multiply it by 100
  fire_agreement <- (dist_B*cos(A) / dist_B) * 100
}


# Between Different Life Stages:

lonlat_J <- plts[plts$CN %in% plt_difs_JG[plt_difs_JG$NVL_SEED == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_J)<-c('x','y')

lonlat_G <- plts[plts$CN %in% plt_difs_JG[plt_difs_JG$BOTH == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_G)<-c('x','y')
clim12 <- env.lonlat

occ.sp <- rbind(lonlat_G,lonlat_J)
occ.G <- na.exclude(sample.sp.globvar(dfsp=lonlat_G,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=resolution_raster))
occ.J <- na.exclude(sample.sp.globvar(dfsp=lonlat_J,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=resolution_raster))

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

clim_centroid_scaled_NG <- colMeans(
  filter(env.lonlat_scaled_JG,ID=="Both") %>%
  select(-ID)
)
clim_centroid_scaled_NJ <- colMeans(
  filter(env.lonlat_scaled_JG,ID=="Novel") %>%
  select(-ID)
)

clim_dif_JG <- clim_centroid_scaled_NJ - clim_centroid_scaled_NG



clim_dif_ST <- clim_dif_N # Calculated in 0-primary_analysis (climatic dif of seedling/tree outside of fire)
dist_STJG <- dist(rbind(clim_dif_JG, clim_dif_ST),method="euclidean") # distance between vector termini
dist_ST <- dist_N # Calculated line 12 of this script
dist_JG <- dist(rbind(clim_centroid_scaled_NJ, clim_centroid_scaled_NG),method="euclidean")# this is magnitude of sapling, big tree vector

#arccos(A/H)= angle
#use the law of cosines to find the angle A between the vectors: 
A <- as.numeric(acos((dist_ST^2 + dist_JG^2 - dist_STJG^2)/(2*dist_ST*dist_JG)))
if(A > pi/2){
  lifestage_agreement <- 0
} else{
  #this here calculates the proportion of dist_ST that lies along dist_JG, inspired by http://www.nabla.hr/VT-Vectors2Dand3DB5.htm
  # dist_ST*cos(A) is the magnitude of dist_ST that lies along dist_JG
  # then we divide it by the entirety of dist_ST and multiply it by 100
  lifestage_agreement <- (dist_ST*cos(A) / dist_ST) * 100
}
