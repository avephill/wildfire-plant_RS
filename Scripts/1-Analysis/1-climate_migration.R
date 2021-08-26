## This script is called by 0-primary_analysis and is used to
## calculate climatic niche difference in two ways: 
##     Schoener's D and Euclidean centroid distance


setwd(project_directory)

### Inside Fire #######################################
print("Working on Burned Plots")
phase_name <- "In Fire"

lonlat_BS <- plts[plts$CN %in% plt_difs_B[plt_difs_B$NVL_SEED == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_BS)<-c('x','y')

lonlat_BT <- plts[plts$CN %in% plt_difs_B[plt_difs_B$BOTH == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_BT)<-c('x','y')

# source("Scripts/1-Analysis/1.1-broennimann2012_adapted.R")
niche_B.results <- broennimann.niche.analysis(env.lonlat, env.stack, lonlat_BS, lonlat_BT, rep.master)

centroid_B.results <- centroid.analysis(clim.T = niche_B.results$clim.T,
                                        clim.S = niche_B.results$clim.S,
                                        scaling_vector = scaling_vector)
cent_d_B <- centroid_B.results$centroid_distance # centroid distance
cent_d_B_ci <- centroid_B.results$cent_d_ci # confidence interval
cent_B <- centroid_B.results$centroids # centroid coords
pca_B <- centroid_B.results$pca_list
h2.test_B <- centroid_B.results$h2.test_result
# for centroid difference difference
Bboth.df_scaled <- centroid_B.results$both.df_scaled 
Bnovel.df_scaled <- centroid_B.results$novel.df_scaled

env.lonlat_BT <- centroid_B.results$env.lonlat_T %>% mutate(fire = phase_name)
env.lonlat_BS <- centroid_B.results$env.lonlat_S %>% mutate(fire = phase_name)




### Outside Fire #######################################
print("Working on Unburned Plots")
phase_name<- "Outside Fire"

lonlat_NS<- plts[plts$CN %in% plt_difs_N[plt_difs_N$NVL_SEED == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_NS)<-c('x','y')

lonlat_NT<-plts[plts$CN %in% plt_difs_N[plt_difs_N$BOTH == species,]$PLT_CN,c('LON',"LAT")]
colnames(lonlat_NT)<-c('x','y')


niche_N.results <- broennimann.niche.analysis(env.lonlat, env.stack, lonlat_NS, lonlat_NT, rep.master)


centroid_N.results <- centroid.analysis(clim.T = niche_N.results$clim.T,
                                        clim.S = niche_N.results$clim.S,
                                        scaling_vector = scaling_vector)
cent_d_N <- centroid_N.results$centroid_distance # centroid distance
cent_d_N_ci <- centroid_N.results$cent_d_ci
cent_N <- centroid_N.results$centroids # centroid coords
pca_N <- centroid_N.results$pca_list
h2.test_N <- centroid_N.results$h2.test_result

# for centroid difference difference
Nboth.df_scaled <- centroid_N.results$both.df_scaled 
Nnovel.df_scaled <- centroid_N.results$novel.df_scaled

env.lonlat_NT <- niche_N.results$clim.T %>% mutate(fire = phase_name)
env.lonlat_NS <- niche_N.results$clim.S %>% mutate(fire = phase_name)



# this goes into the master results
df_env_T <- niche_B.results$occ.T %>% mutate(ID = "Both")
df_env_S <- niche_B.results$occ.S %>% mutate(ID = "Novel")
df_env_BST <- bind_rows(df_env_T, df_env_S) %>% mutate(FIRE = "yes")

df_env_T <- niche_N.results$occ.T %>% mutate(ID = "Both")
df_env_S <- niche_N.results$occ.S %>% mutate(ID = "Novel")
df_env_NST <- bind_rows(df_env_T, df_env_S) %>% mutate(FIRE = "no")


df_env_master_i <- bind_rows(df_env_BST, df_env_NST)
df_env_master_i$species <- species
df_env_master_i$scaled <- "no"



print("Working on Dif in Dif stats")
##### Niche Difference Difference Test


difdif.overlap.eq.gen <- function (repi, Bz1, Bz2, Nz1, Nz2, kernel.method)
{

  if (!is.null(Bz1$y)) {
    occ.pool <- rbind(Bz1$sp, Bz2$sp, Nz1$sp, Nz2$sp) %>% data.frame() %>% 
      mutate(ID = 1:nrow(.))
    row.names(occ.pool) <- c()
    spB1.sim <- occ.pool %>% sample_n(nrow(Bz1$sp))
    spB2.sim <- occ.pool %>% setdiff(spB1.sim) %>% sample_n(nrow(Bz2$sp))
    spN1.sim <- occ.pool %>% setdiff(rbind(spB1.sim, spB2.sim)) %>% sample_n(nrow(Nz1$sp))
    spN2.sim <- occ.pool %>% setdiff(rbind(spB1.sim, spB2.sim, spN1.sim)) %>% sample_n(nrow(Nz2$sp))
    
    spB1.sim <- spB1.sim %>% dplyr::select(-ID)
    spB2.sim <- spB2.sim %>% dplyr::select(-ID)
    spN1.sim <- spN1.sim %>% dplyr::select(-ID)
    spN2.sim <- spN2.sim %>% dplyr::select(-ID)
  }
  Bz1.sim <- ecospat.grid.clim.dyn_alt(Bz1$glob, Bz1$glob1, data.frame(spB1.sim),
                                   R = length(Bz1$x), kernel.method = kernel.method,
                                   th.env = .05)
  Bz2.sim <- ecospat.grid.clim.dyn_alt(Bz2$glob, Bz2$glob1, data.frame(spB2.sim),
                                   R = length(Bz2$x), kernel.method = kernel.method,
                                   th.env = .05)
  Nz1.sim <- ecospat.grid.clim.dyn_alt(Nz1$glob, Nz1$glob1, data.frame(spN1.sim),
                                   R = length(Nz1$x), kernel.method = kernel.method,
                                   th.env = .05)
  Nz2.sim <- ecospat.grid.clim.dyn_alt(Nz2$glob, Nz2$glob1, data.frame(spN2.sim),
                                   R = length(Nz2$x), kernel.method = kernel.method,
                                   th.env = .05)
  sim.Ddif_abs <- abs(ecospat.niche.overlap(Bz1.sim, Bz2.sim, cor = TRUE)$D - 
    ecospat.niche.overlap(Nz1.sim, Nz2.sim, cor = TRUE)$D)
  sim.Ddif <- ecospat.niche.overlap(Bz1.sim, Bz2.sim, cor = TRUE)$D - 
                        ecospat.niche.overlap(Nz1.sim, Nz2.sim, cor = TRUE)$D
  return(data.frame(absSim = sim.Ddif_abs, sim = sim.Ddif))
}


Bz1 <- niche_B.results$z1
Bz2 <- niche_B.results$z2
Nz1 <- niche_N.results$z1
Nz2 <- niche_N.results$z2
R <- length(Bz1$x)
difdifD.results <- list()
obs.Ddif <- ecospat.niche.overlap(Bz1, Bz2, cor = TRUE)$D  -
  ecospat.niche.overlap(Nz1, Nz2, cor = TRUE)$D 
  

obs.Ddif_abs <- obs.Ddif %>% abs()

cl <- makeCluster(detectCores() - 1)
invisible(clusterEvalQ(cl, library("ecospat")))
invisible(clusterEvalQ(cl, library("dplyr")))
invisible(clusterEvalQ(cl, source("Scripts/1-Analysis/1.1.2-ecospat_adapted.R")))
sim.o <- matrix(unlist(parLapply(cl, 1:rep.master,
                                 difdif.overlap.eq.gen, Bz1, Bz2, Nz1, Nz2, "ks")), 
                byrow = TRUE, ncol = 2) %>% 
  data.frame()
stopCluster(cl)

colnames(sim.o) <- c("D_abs", "D")
difdifD.results$sim <- sim.o
difdifD.results$rawObs <- obs.Ddif
difdifD.results$obs <- list(D = obs.Ddif_abs)

alternative <- "greater"
if (alternative == "greater") {
  difdifD.results$p.D <- (sum(sim.o$D_abs >= obs.Ddif_abs) + 1)/(length(sim.o$D_abs) + 1)
}

ecospat.plot.overlap.test(difdifD.results, "D", paste(species, "Ndif - Bdif"))
difdif.plot.overlap(difdifD.results$sim$D, difdifD.results$rawObs, 
                    p = 1, title = "B_D - N_D")
difdif.plot.overlap(difdifD.results$sim$D_abs, difdifD.results$obs$D, 
                    p = difdifD.results$p.D, title = "B_D - N_D Absolute")

####




### Centroid Difference Difference Test ####

centDD_bootstrap.df <- rbind(Bboth.df_scaled  %>% mutate(ID = "BT" ),
                      Bnovel.df_scaled %>% mutate(ID = "BS"),
                      Nboth.df_scaled  %>% mutate(ID = "NT" ),
                      Nnovel.df_scaled %>% mutate(ID = "NS")) %>% 
  mutate(ID = factor(ID))

centroid.DD <- function(data, indices){
  data_samp <- data[indices,]
  
  BT.centroid <- data_samp %>% dplyr::filter(ID == "BT") %>% 
    dplyr::select(-ID) %>% colMeans()
  BS.centroid <- data_samp %>% dplyr::filter(ID == "BS") %>% 
    dplyr::select(-ID) %>% colMeans()
  Bcentroids <- rbind(BT.centroid, BS.centroid)
  
  Bcent_d <- as.numeric(dist(Bcentroids, method="euclidean"))
  
  NT.centroid <- data_samp %>% dplyr::filter(ID == "NT") %>% 
    dplyr::select(-ID) %>% colMeans()
  NS.centroid <- data_samp %>% dplyr::filter(ID == "NS") %>% 
    dplyr::select(-ID) %>% colMeans()
  Ncentroids <- rbind(NT.centroid, NS.centroid)
  
  Ncent_d <- as.numeric(dist(Ncentroids, method="euclidean"))
  cent_dD <- Bcent_d - Ncent_d
  
  return(cent_dD)
}


centDD_boot <- boot(centDD_bootstrap.df, centroid.DD, R = 20000, parallel = "multicore",
                    strata = centDD_bootstrap.df$ID)



