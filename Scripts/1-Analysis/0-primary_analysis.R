# This is the primary script for data analysis, where Schoener's D and 
# climate-space centroid distance are calculated for each species, and averaged across them.
# Main results can be found in the results_list produced at the end of this script, and are 
# further distilled into results_df
# This script relies upon both climate_raster_prep.R and occurrence_data_prep.R being run


if(!exists("project_directory")) project_directory <- "~/OfficePortal/ClimateMismatch/Publishable Repo" # Enter location of project directory within quotation marks
setwd(project_directory)

source("Scripts/1-Analysis/1.1.2-broennimann_functions.R")

library(rgdal)
library(raster)
library(ade4)
library(adehabitatHS)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(knitr)
library(markdown)
library(reshape2)
library(pander)
library(ggrepel)
library(cowplot)


rep.master <- 100 # How many simulations runs to test statistical significance of Schoener's D
N_threshold  <-  4 # Only uses species with more than 'threshold' number of occurrences in each data class
bool_plot_climate_change <- T #this makes a plot of the average climatic changes in the last 30 years
bool_override_species <- F
override_species <- c('Pinus ponderosa',"Pseudotsuga menziesii")

# Determines the extent to which the migration vectors much agree in order to be included in analysis (in %)
vec_agreement_cutoff_fire <- 50 
vec_agreement_cutoff_lifestage <- 50


## Load in occurrence Data
setwd("Data/Occurrence/Processed")
# This dataframe indicates the species compisition of the follow categories for each plot:
# Only Adult Trees, Only Saplings, or Both
prep_plt_difs<-read.csv("plot_differences.csv")

# Dataframe with all plots and fire information
prep_plts <- read.csv("pl_fire.csv",stringsAsFactors = F) 

# This is a dataset that indicates what plots have, of each species:
# Only big mature Trees, Only saplings, or Both
prep_plt_difs_JG <- read.csv("plot_differences_sapling_mature.csv",stringsAsFactors = F) 

# We don't want fire in the sapling/mature analysis
prep_plt_difs_JG <- subset(prep_plt_difs_JG,is.na(prep_plt_difs_JG$FIRE)) 

spcd<-read.csv("../Original/species_codes.csv") #this translates species ID's to species names
tr_lite<-read.csv("tr_lite.csv",stringsAsFactors = F)
sd_lite<-read.csv("sd_lite.csv",stringsAsFactors = F)



## Load in Environmental Data
crs_wgs84<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
setwd(project_directory)
setwd("Data/Environmental")
# This is the directory in which I'll find all 
# the environmental variables used in the study
dir_climate <- "Processed/CMIP8110"

# This is the directory in which I'll find all the 
# DIFFERENCE in environmental variables used in the study
dir_cc <- "Processed/CMIPdif"



# Here we trim spatial data (occurrence and environmental) to parameters specified
# This R script uses crs_wgs84, dir_climate, dir_cc, prep_plts, prep_plt_difs_JG, prep_plt_difs
# produces plts, plt_difs, plt_difs_JG, env.lonlat, cc_env.lonlat, scale_vector
setwd(project_directory)
source("Scripts/1-Analysis/0.1-spatial_data_prep.R")

#####
setwd(project_directory)

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
  
  if(min(c(n_BS,n_BT,n_NS,n_NT)) > N_threshold){
    # If the species meets the set sample size threshold, then Proceed with analysis
    
    migration="Immigration"
    
    # The following Script calculates Schoener's D and centroid distance
    source("Scripts/1-Analysis/1-climate_migration.R")
    
    results_list_i <- list(SchoenerD_N = a_N, SchoenerD_B = a_B, 
                           centroid_dist_N=cent_d_N, centroid_dist_B=cent_d_B,
                           centroid_coords_N = cent_N, centroid_coords_B = cent_B, 
                           pca_N = pca_N, pca_B = pca_B, N = N)
    
    
    
    # Scale the species climate niche with global climate 
    df_env_master_i_scaled <- as.data.frame(scale(df_env_master_i[,!names(df_env_master_i) %in% c('x','y','ID','FIRE',"scaled","species")],center=F,scale=scaling_vector))
    df_env_master_i_scaled$x <- df_env_master_i$x
    df_env_master_i_scaled$y <- df_env_master_i$y
    df_env_master_i_scaled$ID <- df_env_master_i$ID
    df_env_master_i_scaled$FIRE <- df_env_master_i$FIRE
    
    clim_scaled_BT <- 
      df_env_master_i_scaled %>%
      filter(ID == "Both" & FIRE == "yes") %>% 
      select(-c("ID","FIRE","x","y"))
    
    clim_scaled_BS <- 
      df_env_master_i_scaled %>%
      filter(ID == "Novel" & FIRE == "yes") %>% 
      select(-c("ID","FIRE","x","y"))

    clim_scaled_NT <- 
      df_env_master_i_scaled %>%
      filter(ID == "Both" & FIRE == "no") %>% 
      select(-c("ID","FIRE","x","y"))
    
    clim_scaled_NS <- 
      df_env_master_i_scaled %>%
      filter(ID == "Novel" & FIRE == "no") %>% 
      select(-c("ID","FIRE","x","y"))
  
    # Calculate centroid difference between scaled climatic niche of S and T for both B and N 
    clim_dif_N <- colMeans(clim_scaled_NS) - colMeans(clim_scaled_NT)
    clim_dif_B <- colMeans(clim_scaled_BS) - colMeans(clim_scaled_BT)
    results_list_i$clim_centroid_dif_N <- clim_dif_N
    results_list_i$clim_centroid_dif_B <- clim_dif_B
    
    df_env_master_i_scaled$species <- species
    df_env_master_i_scaled$scaled <- "yes"
    if(exists("df_env_master")){
      df_env_master <- bind_rows(df_env_master,df_env_master_i,df_env_master_i_scaled)
    } else {
      df_env_master <- bind_rows(df_env_master_i,df_env_master_i_scaled)
    }
    
    # Calculate the degree to which migration direction is consistent 
    # between life stages and between burned/unburned areas
    source("Scripts/1-Analysis/2-migration_vec_agreement.R")
    results_list_i$lifestage_agreement <- as.numeric(lifestage_agreement)
    results_list_i$fire_agreement <- as.numeric(fire_agreement)
    results_list$new_spec <- results_list_i
    names(results_list)[names(results_list)=="new_spec"] <- species
    
  } else{
    print("Minimum sample size not met")
  }
}

dir.create("Results")
save(results_list, file="Results/results_list.RData")
write.csv(df_env_master,"Results/df_env_master.csv")

# Distills results_list into Dataframe 
results_vec <- unlist(results_list)
results_df <- data.frame(
  species = names(results_list),
  fire_agreement = as.numeric(unlist(lapply(results_list,`[`,c("fire_agreement")))),
  lifestage_agreement = as.numeric(unlist(lapply(results_list,`[`,c("lifestage_agreement")))),
  Schoener_D_N = as.numeric(results_vec[grepl("N.obs.D",names(results_vec))]),
  Schoener_Dp_N = as.numeric(results_vec[grepl("N.p.D",names(results_vec))]),
  Schoener_D_B = as.numeric(results_vec[grepl("B.obs.D",names(results_vec))]),
  Schoener_Dp_B = as.numeric(results_vec[grepl("B.p.D",names(results_vec))]),
  centroid_dist_N = as.numeric(results_vec[grepl("centroid_dist_N",names(results_vec))]),
  centroid_dist_B = as.numeric(results_vec[grepl("centroid_dist_B",names(results_vec))]),
  n_NS = as.numeric(results_vec[grepl("n_NS",names(results_vec))]),
  n_NT = as.numeric(results_vec[grepl("n_NT",names(results_vec))]),
  n_BS = as.numeric(results_vec[grepl("n_BS",names(results_vec))]),
  n_BT = as.numeric(results_vec[grepl("n_BT",names(results_vec))])
)
