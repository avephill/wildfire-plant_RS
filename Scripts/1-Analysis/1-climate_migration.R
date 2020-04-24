## This script is called by 0-primary_analysis and is used to
## calculate climatic niche difference in two ways: 
##     Schoener's D and Euclidean centroid distance


setwd(project_directory)

### Inside Fire #######################################
phase_name <- "In Fire"

lonlat_S <- plts[plts$CN %in% plt_difs_B[plt_difs_B$NVL_SEED == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_S)<-c('x','y')

lonlat_T <- plts[plts$CN %in% plt_difs_B[plt_difs_B$BOTH == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_T)<-c('x','y')

source("Scripts/1-Analysis/1.1-broennimann2012_adapted.R",echo = T)
a_B <- a # this is the result of niche equivalency test
occ_BT <- occ.T
occ_BS <- occ.S

source("Scripts/1-Analysis/1.2-centroid_analysis.R") # uses occ.T and occ.S from 1.1-broennimann_adapted.R
cent_d_B <- centroid_results$centroid_distance # centroid distance
cent_B <- centroid_results$centroids # centroid coords
pca_B <- centroid_results$pca_list

env.lonlat_S$fire <- phase_name
env.lonlat_T$fire <- phase_name

env.lonlat_BT <- env.lonlat_T
env.lonlat_BS <- env.lonlat_S

# this goes into the master results
df_env_T <- occ.T
df_env_T$ID <- "Both"
df_env_S <- occ.S
df_env_S$ID <- "Novel"
df_env_ST <- bind_rows(df_env_T,df_env_S)
df_env_ST$FIRE <- "yes"
df_env_BST <- df_env_ST



### Outside Fire #######################################
phase_name<- "Outside Fire"

lonlat_S<- plts[plts$CN %in% plt_difs_N[plt_difs_N$NVL_SEED == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_S)<-c('x','y')

lonlat_T<-plts[plts$CN %in% plt_difs_N[plt_difs_N$BOTH == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_T)<-c('x','y')


source("Scripts/1-Analysis/1.1-broennimann2012_adapted.R")
a_N <- a # this is the result of niche equivalency test
occ_NT <- occ.T
occ_NS <- occ.S

source("Scripts/1-Analysis/1.2-centroid_analysis.R") # uses occ.T and occ.S from 1.1-broennimann_adapted.R
cent_d_N <- centroid_results$centroid_distance # centroid distance
cent_N <- centroid_results$centroids # centroid coords
pca_N <- centroid_results$pca_list

env.lonlat_S$fire <- phase_name
env.lonlat_T$fire <- phase_name

env.lonlat_NT <- env.lonlat_T
env.lonlat_NS <- env.lonlat_S


# this goes into the master results
df_env_T <- occ.T
df_env_T$ID <- "Both"
df_env_S <- occ.S
df_env_S$ID <- "Novel"
df_env_ST <- bind_rows(df_env_T,df_env_S)
df_env_ST$FIRE <- "no"
df_env_NST <- df_env_ST



df_env_master_i <- bind_rows(df_env_BST,df_env_NST)
df_env_master_i$species <- species
df_env_master_i$scaled <- "no"
