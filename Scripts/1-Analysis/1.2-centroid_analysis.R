# This file makes the PCA plots for each species, and 
# calculates centroid distance from scaled climatic values
# called from 0-primary_analysis

env.lonlat_T <- occ.T
env.lonlat_S <- occ.S
env.lonlat_T$ID <- "Both"
env.lonlat_S$ID <- "Novel"
pca_df_prep <- rbind(env.lonlat_T,env.lonlat_S)
pca_df <- pca_df_prep[,!names(pca_df_prep) %in% c('x','y')]
last_row_T <- nrow(pca_df[pca_df$ID=="Both",])

pca_df_scaled <- data.frame(scale(pca_df[,!names(pca_df) == "ID"],center=F,scale=scaling_vector))
pca_df_scaled$ID <- pca_df$ID
cent_T <-  colMeans(pca_df_scaled[(pca_df_scaled$ID=="Both") ,!names(pca_df_scaled)=="ID"])
cent_S <- colMeans(pca_df_scaled[(pca_df_scaled$ID=="Novel") ,!names(pca_df_scaled)=="ID"])
centroids <- rbind(cent_T,cent_S)
row.names(centroids) <- c("Both","Novel")

cent_d <- as.numeric(dist(centroids,method="euclidean"))
centroids <- t(data.frame(centroids))

## Only for visualization
pca <- prcomp(pca_df[,unlist(lapply(pca_df, FUN=class))=="numeric"], center=T, scale.=scaling_vector)
maintext<- paste("Immigration",phase_name)
##
pca_list <- list(pca_results=pca,pca_source=pca_df,pca_text=maintext)
centroid_results <- list(centroids=centroids, centroid_distance=cent_d, pca_list=pca_list)

