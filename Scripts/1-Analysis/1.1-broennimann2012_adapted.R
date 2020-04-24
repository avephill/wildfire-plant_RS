##
## This script is adapted from a script originally written by Olivier Broennimann. 
## Departement of Ecology and Evolution (DEE). 
## 2nd of November 2010. University of Lausanne. Department of ecology and evolution, Switzerland
##
## DESCRIPTION
## This script investigates ecological niche overlap from occurrence and spatial environmental data
## 


#################################################################################################
############################## preparation of datasets ##########################################
#################################################################################################

# total climate space
clim12 <- env.lonlat 

# Occurrences of both S and T
occ.sp <- rbind(lonlat_T,lonlat_S)

# create sp occurrence dataset by adding climate variables from the global climate datasets
# resolution should be the resolution of the climate data grid
occ.T <- na.exclude(sample.sp.globvar(dfsp=lonlat_T,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=resolution_raster))
occ.S <- na.exclude(sample.sp.globvar(dfsp=lonlat_S,colspxy=1:2,colspkept=NULL,dfvar=clim12,colvarxy=1:2,colvar="all",resolution=resolution_raster))


#################################################################################################
############################## ANALYSIS - selection of parameters ###############################
#################################################################################################


# selection of variables to include in the analyses
nvar <- length(layer_names)

#number of interation for the tests of equivalency and similarity defined in 0-primary_analysis.R


#resolution of the gridding of the climate space
R=100

#################################################################################################
################### row weigthing and grouping factors for ade4 functions  ######################
#################################################################################################

row.w.1.occ<-1-(nrow(occ.T)/nrow(rbind(occ.T,occ.S))) # relative prevalence of occ1 to occ2
row.w.2.occ<-1-(nrow(occ.S)/nrow(rbind(occ.T,occ.S))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim12)),rep(row.w.1.occ, nrow(occ.T)),rep(row.w.2.occ, nrow(occ.S)))

#######################
# global dataset for the analysis and rows for each sub dataset

data.env.occ <- rbind(clim12, occ.T,occ.S) %>%
  select(-x,-y)

row.clim12 <- 1:nrow(clim12)
row.T<- (nrow(clim12)+1):(nrow(clim12)+nrow(occ.T))
row.S<- (nrow(clim12)+nrow(occ.T) + 1): (nrow(clim12)+nrow(occ.T) + nrow(occ.S))


#################################### PCA-ENV ####################################################
#################################################################################################

# measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas


pca.cal <- dudi.pca(data.env.occ, center = T, scale = T, scannf = F, nf = 2)

# predict the scores on the axes
scores.clim12<- pca.cal$li[row.clim12,] # the row coordinates/principle components of each
scores.T<- pca.cal$li[row.T,]
scores.S<- pca.cal$li[row.S,]

# calculation of occurence density and test of niche equivalency and similarity 
glob<-scores.clim12
#this is imported from grid_clim, may speed things up:
xmin<-min(glob[,1]);xmax<-max(glob[,1]);ymin<-min(glob[,2]);ymax<-max(glob[,2])			# The max and mins of the last PCA plot
globr<-data.frame(cbind((glob[,1]-xmin)/abs(xmax-xmin),(glob[,2]-ymin)/abs(ymax-ymin)))	# makes the subset climate between 0 and 1, relative to the min and max of global climate	
xyR<- cbind((1:R)/R,(1:R)/R)   #this is for devoping a grid
coordsR<-SpatialPoints(xyR)   #this is for devoping a grid
mask<-ascgen(coordsR,nrcol=98,count=F)   #this is for devoping a grid	
glob.dens<-kernelUD(SpatialPoints(globr[,1:2]) ,grid=mask,kern="bivnorm")
glob.dens<-glob.dens[1]

#I think the following makes axes that span the entirety of the global climate PCA
#probably for setting plot boundaries later
x<-seq(from=min(glob[,1]),to=max(glob[,1]),length.out=R)				# breaks on score gradient 1 
y<-seq(from=min(glob[,2]),to=max(glob[,2]),length.out=R)				# breaks on score gradient 2


z1 <- grid.clim(glob=scores.clim12,sp=scores.T,R)

z2 <- grid.clim(glob=scores.clim12,sp=scores.S,R)

a <- niche.equivalency.test(z1,z2,rep=rep.master)# test of niche equivalency and similarity according to Warren et al. 2008
