# This is the primary script for data analysis, where Schoener's D and 
# climate-space centroid distance are calculated for each species, and averaged across them.
# Main results can be found in the results_list produced at the end of this script, and are 
# further distilled into results_df
# This script relies upon both climate_raster_prep.R and occurrence_data_prep.R being run


if(!exists("project_directory")) project_directory <- "" # Enter location of project directory within quotation marks
setwd(project_directory)

library(rgdal)
library(raster)
library(ade4)
library(adehabitatHS)
library(tidyverse)
library(ggfortify)
library(knitr)
library(markdown)
library(pander)
library(ggrepel)
library(cowplot)
library(MVN) # For testing multivariate normality
library(ICSNP)
library(parallel)
library(ecospat)


rep.master <- 250  # How many simulations runs to test statistical significance of Schoener's D
N_threshold  <-  4 # Only uses species with more than 'threshold' number of occurrences in each data class
bool_override_species <- F # Just perform analysis on specified species (skip vetting)
#override_species <- c()
resultsDir <- "Results3Clims"

# Determines the extent to which the migration vectors must agree in order to be included in analysis (in %)
vec_agreement_cutoff_fire <- 50 
vec_agreement_cutoff_lifestage <- 50

source("Scripts/1-Analysis/0.1-spatial_data_prep.R", echo = T)

dir.create(resultsDir)
write.csv(env.lonlat, paste0(resultsDir,"/env_lonlat.csv"), row.names = F)

#####
source("Scripts/1-Analysis/1.1.2-ecospat_adapted.R", echo = T)

species_list <- as.character(unique(na.omit(plt_difs$BOTH)))
species_list <-  species_list[order(species_list,decreasing = F)] 

if(bool_override_species == T){
  species_list <- override_species
}
results_list <- list()
for (k in 1:length(species_list)){
  species <- species_list[k]
  species_ID <- spcd[spcd$SCI_NAME==species,]$SPCD
  print(species)
  
  # N = not burned, B = burned, T = trees and seedlings, S = seedlings only
  
  # Isolate all the plots without fire
  plt_difs_N <- subset(plt_difs,is.na(plt_difs$FIRE))
  n_NS <- length(unique(na.omit(plt_difs_N[plt_difs_N$NVL_SEED == species,]$PLT_CN)))
  n_NT <- length(unique(na.omit(plt_difs_N[plt_difs_N$BOTH == species,]$PLT_CN)))
  
  #all the plots with fire
  plt_difs_B <- subset(plt_difs,plt_difs$FIRE=="yes")
  n_BS <- length(unique(na.omit(plt_difs_B[plt_difs_B$NVL_SEED == species,]$PLT_CN)))
  n_BT <- length(unique(na.omit(plt_difs_B[plt_difs_B$BOTH == species,]$PLT_CN)))
  N <- list(n_NS = n_NS,
            n_NT = n_NT,
            n_BS = n_BS,
            n_BT = n_BT)
  
  if(min(unlist(N)) > N_threshold){
    # If the species meets the set sample size threshold, then Proceed with analysis
    
    migration="Immigration"
    
    # The following Script calculates Schoener's D and centroid distance
    source("Scripts/1-Analysis/1-climate_migration.R")
    
    results_list_i <- list(SchoenerD_N = niche_N.results$a, SchoenerD_B = niche_B.results$a, 
                           centroid_dist_N = cent_d_N, centroid_dist_N_ci = cent_d_N_ci,
                           centroid_dist_B = cent_d_B, centroid_dist_B_ci = cent_d_B_ci,
                           centroid_coords_N = cent_N, centroid_coords_B = cent_B,
                           cent_difdif = centDD_boot, 
                           pca_N = pca_N, pca_B = pca_B, h2.test_B = h2.test_B, 
                           h2.test_N = h2.test_N, N = N, difdifD = difdifD.results)
    
    
    
    # Scale the species climate niche with global climate 
    df_env_master_i_scaled <- as.data.frame(
      scale(df_env_master_i %>% dplyr::select(-x, -y, -ID, -FIRE, 
                                              -scaled, -species),
            center=F, scale=scaling_vector))
    df_env_master_i_scaled$x <- df_env_master_i$x
    df_env_master_i_scaled$y <- df_env_master_i$y
    df_env_master_i_scaled$ID <- df_env_master_i$ID
    df_env_master_i_scaled$FIRE <- df_env_master_i$FIRE
    
    clim_scaled_BT <- 
      df_env_master_i_scaled %>%
      filter(ID == "Both" & FIRE == "yes") %>% 
      dplyr::select(-c("ID","FIRE","x","y"))
    
    clim_scaled_BS <- 
      df_env_master_i_scaled %>%
      filter(ID == "Novel" & FIRE == "yes") %>% 
      dplyr::select(-c("ID","FIRE","x","y"))

    clim_scaled_NT <- 
      df_env_master_i_scaled %>%
      filter(ID == "Both" & FIRE == "no") %>% 
      dplyr::select(-c("ID","FIRE","x","y"))
    
    clim_scaled_NS <- 
      df_env_master_i_scaled %>%
      filter(ID == "Novel" & FIRE == "no") %>% 
      dplyr::select(-c("ID","FIRE","x","y"))
  
    # Calculate centroid difference between scaled climatic niche of S and T for both B and N 
    clim_dif_N <- colMeans(clim_scaled_NS) - colMeans(clim_scaled_NT)
    clim_dif_B <- colMeans(clim_scaled_BS) - colMeans(clim_scaled_BT)
    results_list_i$clim_centroid_dif_N <- clim_dif_N
    results_list_i$clim_centroid_dif_B <- clim_dif_B
    
    df_env_master_i_scaled$species <- species
    df_env_master_i_scaled$scaled <- "yes"
    if(exists("df_env_master")){
      df_env_master <- bind_rows(df_env_master, df_env_master_i, df_env_master_i_scaled)
    } else {
      df_env_master <- bind_rows(df_env_master_i,df_env_master_i_scaled)
    }
    
    # Calculate the degree to which migration direction is consistent 
    # between life stages and between burned/unburned areas
    source("Scripts/1-Analysis/2-migration_vec_agreement.R")
    results_list_i$lifestage_agreement <- as.numeric(lifeagreement)
    results_list_i$lifestage_agreement_SE <- as.numeric(lifeagreement_se)
    results_list_i$fire_agreement <- as.numeric(burnagreement)
    results_list_i$fire_agreement_SE <- as.numeric(burnagreement_se)
    results_list$new_spec <- results_list_i
    names(results_list)[names(results_list)=="new_spec"] <- species
    
  } else{
    print("Minimum sample size not met")
  }
}

# Computes approximate altitudinal spatial distance from temperature 
# and adds it to resultsList

source("Scripts/1-Analysis/3-spatial_distance_conversion.R")

dir.create(resultsDir)

save(results_list, file=paste(resultsDir,"/results_list.RData", sep=""))
write.csv(df_env_master,paste(resultsDir, "/df_env_master.csv", sep=""))

# Distills results_list into Dataframe 
results_vec <- unlist(results_list)
results_df <- data.frame(
  species = names(results_list),
  fire_agreement = as.numeric(unlist(lapply(results_list,`[`,c("fire_agreement")))),
  fire_agreement_SE = as.numeric(unlist(lapply(results_list,`[`,c("fire_agreement_SE")))),
  lifestage_agreement = as.numeric(unlist(lapply(results_list,`[`,c("lifestage_agreement")))),
  lifestage_agreement_SE = as.numeric(unlist(lapply(results_list,`[`,c("lifestage_agreement_SE")))),
  Schoener_D_N = as.numeric(results_vec[grepl("N.obs.D",names(results_vec))]),
  Schoener_Dp_N = as.numeric(results_vec[grepl("N.p.D",names(results_vec))]),
  Schoener_D_B = as.numeric(results_vec[grepl("B.obs.D",names(results_vec))]),
  Schoener_Dp_B = as.numeric(results_vec[grepl("B.p.D",names(results_vec))]),
  centroid_dist_N = as.numeric(results_vec[grepl("centroid_dist_N$",names(results_vec))]),
  h2.test_N = as.numeric(results_vec[grepl("h2.test_N.p.value",names(results_vec))]),
  centroid_dist_B = as.numeric(results_vec[grepl("centroid_dist_B$",names(results_vec))]),
  h2.test_B = as.numeric(results_vec[grepl("h2.test_B.p.value",names(results_vec))]),
  n_NS = as.numeric(results_vec[grepl("n_NS",names(results_vec))]),
  n_NT = as.numeric(results_vec[grepl("n_NT",names(results_vec))]),
  n_BS = as.numeric(results_vec[grepl("n_BS",names(results_vec))]),
  n_BT = as.numeric(results_vec[grepl("n_BT",names(results_vec))]),
  alt.dist_B = as.numeric(results_vec[grepl("alt.dist_B",names(results_vec))]),
  alt.dist_N = as.numeric(results_vec[grepl("alt.dist_N",names(results_vec))]),
  difdifD = as.numeric(results_vec[grepl("difdifD.rawObs",names(results_vec))]),
  difdifDp = as.numeric(results_vec[grepl("difdifD.p.D",names(results_vec))])
)

# Calculate CI for centroids difference in Difference
results_df$cent_difdif <- lapply(results_list, function(x){
  x$cent_difdif <- x$cent_difdif$t0}) %>% 
  unlist()

confType <- 'basic'
results_df$cent_difdifLower <- lapply(results_list, function(x){
  x$cent_difdifLower <- boot.ci(x$cent_difdif, type = confType)[[confType]][4]}) %>% 
  unlist()
results_df$cent_difdifUpper <- lapply(results_list, function(x){
  x$cent_difdifUpper <- boot.ci(x$cent_difdif, type = confType)[[confType]][5]}) %>% 
  unlist()

results_df$cent_difdif90Lower <- lapply(results_list, function(x){
  x$cent_difdif90Lower <- boot.ci(x$cent_difdif, type = confType, conf = .90)[[confType]][4]}) %>% 
  unlist()

results_df$cent_difdif90Upper <- lapply(results_list, function(x){
  x$cent_difdif90Upper <- boot.ci(x$cent_difdif, type = confType, conf = .90)[[confType]][5]}) %>% 
  unlist()


write.csv(results_df, paste0(resultsDir, "/results_df.csv"), row.names = F)
View(results_df)


