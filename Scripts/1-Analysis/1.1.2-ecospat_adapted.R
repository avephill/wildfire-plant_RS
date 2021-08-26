
broennimann.niche.analysis <- function(env.lonlat, env.stack, lonlat_S, lonlat_T, rep.master){
  #resolution of the gridding of the climate space
  R=100
  
  clim.bkg <- env.lonlat %>% data.frame() %>% dplyr::select(-x,-y)
  occ.T <- raster::extract(env.stack, lonlat_T, df = T)  %>% 
    dplyr::select(-ID) %>%  mutate(x = lonlat_T$x, y = lonlat_T$y) %>% na.omit 
  occ.T <- occ.T %>% 
    dplyr::select(x, y, names(occ.T %>% dplyr::select(-x,-y)))
  
  occ.S <- raster::extract(env.stack, lonlat_S, df = T) %>% 
    dplyr::select(-ID) %>%  mutate(x = lonlat_S$x, y = lonlat_S$y) %>% na.omit 
  occ.S <- occ.S %>% 
    dplyr::select(x, y, names(occ.S %>% dplyr::select(-x,-y)))
  
  clim.T <- occ.T %>% dplyr::select(-x, -y)
  clim.S <- occ.S %>% dplyr::select(-x, -y)
  
  pca.env <- dudi.pca(clim.bkg, center = T, scale = T, scannf = F, nf = 2)
  
  scores.bkg <- pca.env$li	#scores for global climate
  scores.T <- suprow(pca.env, clim.T)$lisup				#scores for spa
  scores.S <- suprow(pca.env, clim.S)$lisup				#scores for spb
  
  z1 <- ecospat.grid.clim.dyn_alt(glob = scores.bkg, 
                                  glob1 = scores.bkg,
                                  sp = scores.T, R = R, kernel.method = "ks",
                                  th.env = 0.05)
  
  
  z2 <- ecospat.grid.clim.dyn_alt(glob = scores.bkg, 
                                  glob1 = scores.bkg,
                                  sp = scores.S, R = R, kernel.method = "ks",
                                  th.env = 0.05)
  
  
  a <- ecospat.niche.equivalency.test_alt(z1, z2, rep = rep.master,
                                          alternative = "lower",
                                          ncores = detectCores() - 1, kernel.method = "ks",
                                          th.env = 0.05)
  
  return(list(a = a, occ.T = occ.T, occ.S = occ.S, z1 = z1, z2 = z2, 
              clim.T = clim.T, clim.S = clim.S))
  
}


ecospat.niche.equivalency.test_alt <- function (z1, z2, rep, alternative = "greater", ncores = 1,
                                                kernel.method = "adehabitat", th.sp = 0, th.env = 0) 
{
  R <- length(z1$x)
  l <- list()
  obs.o <- ecospat.niche.overlap(z1, z2, cor = TRUE)
  if (ncores == 1) {
    sim.o <- as.data.frame(matrix(unlist(lapply(1:rep, overlap.eq.gen_alt, 
                                                z1, z2, kernel.method, th.sp, th.env)), byrow = TRUE, ncol = 2))
  }
  else {
    cl <- makeCluster(ncores)
    invisible(clusterEvalQ(cl, library("ecospat")))
    invisible(clusterEvalQ(cl, source("Scripts/1-Analysis/1.1.2-ecospat_adapted.R")))
    sim.o <- as.data.frame(matrix(unlist(parLapply(cl, 1:rep, 
                                                   overlap.eq.gen_alt, z1, z2, kernel.method, th.sp, th.env)), byrow = TRUE, ncol = 2))
    stopCluster(cl)
  }
  colnames(sim.o) <- c("D", "I")
  l$sim <- sim.o
  l$obs <- obs.o
  if (alternative == "greater") {
    l$p.D <- (sum(sim.o$D >= obs.o$D) + 1)/(length(sim.o$D) + 
                                              1)
    l$p.I <- (sum(sim.o$I >= obs.o$I) + 1)/(length(sim.o$I) + 
                                              1)
  }
  if (alternative == "lower") {
    l$p.D <- (sum(sim.o$D <= obs.o$D) + 1)/(length(sim.o$D) + 
                                              1)
    l$p.I <- (sum(sim.o$I <= obs.o$I) + 1)/(length(sim.o$I) + 
                                              1)
  }
  return(l)
}

overlap.eq.gen_alt <- function (repi, z1, z2, kernel.method, th.sp, th.env) 
{
  if (is.null(z1$y)) {
    occ.pool <- c(z1$sp, z2$sp)
    rand.row <- sample(1:length(occ.pool), length(z1$sp))
    sp1.sim <- occ.pool[rand.row]
    sp2.sim <- occ.pool[-rand.row]
  }
  if (!is.null(z1$y)) {
    occ.pool <- rbind(z1$sp, z2$sp)
    row.names(occ.pool) <- c()
    rand.row <- sample(1:nrow(occ.pool), nrow(z1$sp))
    sp1.sim <- occ.pool[rand.row, ]
    sp2.sim <- occ.pool[-rand.row, ]
  }
  
  z1.sim <- ecospat.grid.clim.dyn_alt(z1$glob, z1$glob1, data.frame(sp1.sim), 
                                  R = length(z1$x), kernel.method = kernel.method, 
                                  th.sp = th.sp, th.env = th.env)
  z2.sim <- ecospat.grid.clim.dyn_alt(z2$glob, z2$glob1, data.frame(sp2.sim), 
                                  R = length(z2$x), kernel.method = kernel.method, 
                                  th.sp = th.sp, th.env = th.env)
  par(mfrow=c(1,2))
  plot(z1.sim$z.cor)
  plot(z2.sim$z.cor)
  par(mfrow=c(1,1))
  o.i <- ecospat.niche.overlap(z1.sim, z2.sim, cor = TRUE)
  sim.o.D <- o.i$D
  sim.o.I <- o.i$I
  return(c(sim.o.D, sim.o.I))
}



ecospat.grid.clim.dyn_alt <- function (glob, glob1, sp, R = 100, th.sp = 0, th.env = 0, geomask = NULL, 
          kernel.method = "adehabitat", extend.extent = c(0, 0, 0, 
                                                          0)) 
{
  if (is.null(kernel.method) | (kernel.method != "ks" & kernel.method != 
                                "adehabitat")) {
    stop("supply a kernel method ('adehabitat' or 'ks')")
  }
  glob <- as.matrix(glob)
  glob1 <- as.matrix(glob1)
  sp <- as.matrix(sp)
  l <- list()
  if (ncol(glob) > 2) {
    stop("cannot calculate overlap with more than two axes")
  }
  if (ncol(glob) == 1) {
    xmin <- min(glob[, 1]) + extend.extent[1]
    xmax <- max(glob[, 1]) + extend.extent[2]
    glob1.dens <- ecospat.kd_alt(x = glob1, ext = c(xmin, xmax), 
                             method = kernel.method, th = 0, R = R)
    sp.dens <- ecospat.kd_alt(x = sp, ext = c(xmin, xmax), method = kernel.method, 
                          th = 0, env.mask = glob1.dens$y > 0, R = R)
    x <- sp.dens$x
    y <- sp.dens$y
    z <- sp.dens$y * nrow(sp)/sum(sp.dens$y)
    Z <- glob1.dens$y * nrow(glob)/sum(glob1.dens$y)
    z.uncor <- z/max(z)
    z.cor <- z/Z
    z.cor[is.na(z.cor)] <- 0
    z.cor[z.cor == "Inf"] <- 0
    z.cor <- z.cor/max(z.cor)
  }
  if (ncol(glob) == 2) {
    xmin <- apply(glob, 2, min, na.rm = T)
    xmax <- apply(glob, 2, max, na.rm = T)
    ext <- c(xmin[1], xmax[1], xmin[2], xmax[2]) + extend.extent
    glob1.dens <- ecospat.kd_alt(x = glob1, ext = ext, method = kernel.method, 
                             th = th.env, R = R)
    if (!is.null(geomask)) {
      crs(geomask) <- NA
      glob1.dens <- raster::mask(glob1.dens, geomask, 
                                 updatevalue = 0)
    }
    sp.dens <- ecospat.kd_alt(x = sp, ext = ext, method = kernel.method, 
                          th = th.sp, env.mask = glob1.dens > 0, R = R)
    x <- seq(from = ext[1], to = ext[2], length.out = R)
    y <- seq(from = ext[3], to = ext[4], length.out = R)
    l$y <- y
    Z <- glob1.dens * nrow(glob1)/raster::cellStats(glob1.dens, 
                                                    "sum")
    z <- sp.dens * nrow(sp)/raster::cellStats(sp.dens, "sum")
    z.uncor <- z/raster::cellStats(z, "max")
    z.cor <- z/Z
    z.cor[is.na(z.cor)] <- 0
    z.cor <- z.cor/raster::cellStats(z.cor, "max")
  }
  w <- z.uncor
  w[w > 0] <- 1
  l$x <- x
  l$z <- z
  l$z.uncor <- z.uncor
  l$z.cor <- z.cor
  l$Z <- Z
  l$glob <- glob
  l$glob1 <- glob1
  l$sp <- sp
  l$w <- w
  return(l)
}

ecospat.kd_alt <- function (x, ext, R = 100, th = 0, env.mask = c(), method = "adehabitat") 
{
  if (method == "adehabitat") {
    if (ncol(x) == 2) {
      xr <- data.frame(cbind((x[, 1] - ext[1])/abs(ext[2] - 
                                                     ext[1]), (x[, 2] - ext[3])/abs(ext[4] - ext[3])))
      mask <- adehabitatMA::ascgen(sp::SpatialPoints(cbind((0:(R))/R, 
                                                           (0:(R)/R))), nrcol = R - 2, count = FALSE)
      x.dens <- adehabitatHR::kernelUD(sp::SpatialPoints(xr[, 
                                                            1:2]), h = "href", grid = mask, kern = "bivnorm")
      x.dens <- raster::raster(xmn = ext[1], xmx = ext[2], 
                               ymn = ext[3], ymx = ext[4], matrix(x.dens$ud, 
                                                                  nrow = R))
      if (!is.null(th)) {
        th.value <- quantile(raster::extract(x.dens, 
                                             x), th)
        x.dens[x.dens < th.value] <- 0
      }
      if (!is.null(env.mask)) {
        x.dens <- x.dens * env.mask
      }
    }
    else if (ncol(x) == 1) {
      xr <- seq(from = min(ext), to = max(ext), length.out = R)
      x.dens <- density(x[, 1], kernel = "gaussian", from = min(xr), 
                        to = max(xr), n = R, cut = 0)
      if (!is.null(env.mask)) {
        x.dens$y <- x.dens$y * env.mask
      }
      if (!is.null(th)) {
        xr <- sapply(x, findInterval, x.dens$x)
        th.value <- quantile(x.dens$y[xr], th)
        sprm <- which(x.dens$y < th.value)
        x.dens$y[sprm] <- 0
      }
    }
  }
  if (method == "ks") {
    if (ncol(x) == 2) {
      x.dens <- ks::kde(x, xmin = ext[c(1, 3)], xmax = ext[c(2, 
                                                             4)], gridsize = c(R, R))
      x.dens <- raster::flip(raster::t(raster::raster(x.dens$estimate)), 
                             direction = "y")
      raster::extent(x.dens) <- c(xmn = ext[1], xmx = ext[2], 
                                  ymn = ext[3], ymx = ext[4])
      if (!is.null(th)) {
        th.value <- quantile(raster::extract(x.dens, 
                                             x), th)
        x.dens[x.dens < th.value] <- 0
      }
      if (!is.null(env.mask)) {
        x.dens <- x.dens * env.mask
      }
    }
    else if (ncol(x) == 1) {
      x.dens <- ks::kde(x, xmin = min(ext), xmax = max(ext), 
                        gridsize = c(R, R))
      x.dens$y <- x.dens$estimate
      x.dens$x <- x.dens$eval.points
      if (!is.null(env.mask)) {
        x.dens$y <- x.dens$y * env.mask
      }
      if (!is.null(th)) {
        xr <- sapply(x, findInterval, x.dens$x)
        th.value <- quantile(x.dens$y[xr], th)
        sprm <- which(x.dens$y < th.value)
        x.dens$y[sprm] <- 0
      }
    }
  }
  return(x.dens)
}



difdif.plot.overlap <- function (sim, obs, p, title) 
{

  r0 <- c(sim, obs)
  l0 <- max(sim) - min(sim)
  w0 <- l0/(log(length(sim), base = 2) + 1)
  xlim0 <- range(r0) + c(-w0, w0)
  h0 <- hist(sim, plot = FALSE, nclass = 10)
  y0 <- max(h0$counts)
  hist(sim, plot = TRUE, nclass = 10, xlim = xlim0, col = grey(0.8), 
       main = title, xlab = "Difference in D", sub = paste("p.value = ", 
                                              round(p, 5)))
  lines(c(obs, obs), c(y0/2, 0), col = "red")
  points(obs, y0/2, pch = 18, cex = 2, col = "red")
  invisible()
}


# This file makes the PCA plots for each species, and 
# calculates centroid distance from scaled climatic values
# called from 0-primary_analysis
centroid.distance <- function(data, indices){
  
  data_samp <- data[indices,]
  both.centroid <- data_samp %>% dplyr::filter(ID == "Both") %>% 
    dplyr::select(-ID) %>% colMeans()
  novel.centroid <- data_samp %>% dplyr::filter(ID == "Novel") %>% 
    dplyr::select(-ID) %>% colMeans()
  centroids <- rbind(both.centroid, novel.centroid)
  
  cent_d <- as.numeric(dist(centroids, method="euclidean"))
  return(cent_d)
}


centroid.analysis <- function(clim.T, clim.S, scaling_vector){
  env.lonlat_T <- clim.T %>% mutate(ID = "Both")
  env.lonlat_S <- clim.S %>% mutate(ID = "Novel")
  
  pca_df_prep <- rbind(env.lonlat_T, env.lonlat_S)
  pca_df <- pca_df_prep[,!names(pca_df_prep) %in% c('x','y')]
  last_row_T <- nrow(pca_df[pca_df$ID=="Both",])
  
  pca_df_scaled <- data.frame(scale(pca_df[,!names(pca_df) == "ID"], center=F, scale=scaling_vector))
  pca_df_scaled$ID <- pca_df$ID
  
  both.df_scaled <- pca_df_scaled[(pca_df_scaled$ID=="Both") ,!names(pca_df_scaled)=="ID"]
  novel.df_scaled <- pca_df_scaled[(pca_df_scaled$ID=="Novel") ,!names(pca_df_scaled)=="ID"]
  
  bootstrap.df <- rbind(both.df_scaled %>% mutate(ID = "Both"),
                        novel.df_scaled %>% mutate(ID = "Novel")) %>% 
    mutate(ID = factor(ID))
  
  
  cent_boot <- boot(bootstrap.df, centroid.distance, R = 5000, 
                    strata = bootstrap.df$ID,
                    parallel = "multicore")
  
  
  # Let's just use 'chi' for now
  h2.result <- HotellingsT2(both.df_scaled, novel.df_scaled, test = 'chi')
  
  cent_T <-  colMeans(both.df_scaled)
  cent_S <- colMeans(novel.df_scaled)
  centroids <- rbind(cent_T, cent_S)
  row.names(centroids) <- c("Both", "Novel")
  
  cent_d <- as.numeric(dist(centroids, method="euclidean"))
  centroids <- t(data.frame(centroids))
  
  ## Only for visualization
  pca <- prcomp(pca_df[,unlist(lapply(pca_df, FUN=class))=="numeric"], center=T, scale.=scaling_vector)
  maintext<- paste("Immigration", phase_name)
  ##
  pca_list <- list(pca_results = pca, pca_source = pca_df, pca_text = maintext)
  centroid_results <- list(centroids = centroids, centroid_distance = cent_d, 
                           pca_list = pca_list, h2.test_result = h2.result, cent_d_boot = cent_boot,
                           both.df_scaled = both.df_scaled, novel.df_scaled = novel.df_scaled,
                           env.lonlat_T = env.lonlat_T, env.lonlat_S = env.lonlat_S)
  
  return(centroid_results)
}
