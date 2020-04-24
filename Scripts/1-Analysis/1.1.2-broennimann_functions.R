
##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
## remove occurences in a dataframe that are closer to each other than a specified distance threshold
##
##ARGUMENTS
##df: dataframe with x, y, and variables
##colxy: the range of columns for x and y in df
##colvar: the range of columns for variables in df	
##min.dist: minimun distance threshold in the sub-dataframe
##

occ.desaggragation <-function(df,colxy,colvar=NULL,min.dist,plot=T){

initial<-df
train<-initial
xx<-colxy[1]
yy<-colxy[2]
kept<-0 ;out<-0; keep<-c()
x11(2,2,pointsize = 12); par(mar=c(0,0,0,0)); plot.new()

while(nrow(train)>0){
		
	i<-sample(1:nrow(train),1)
		
	if(sum(sqrt((train[,xx]-train[i,xx])^2 + (train[,yy]-train[i,yy])^2)<=min.dist)>1) {
		out<-out+1
		plot.new(); text(0.5,0.8,paste("# initial:",nrow(initial))); text(0.5,0.5,paste("# kept: ",kept)); text(0.5,0.2,paste("# out: ",out))
	}
	else {
		keep<-c(keep,row.names(train[i,]))
		kept<-kept+1
		plot.new(); text(0.5,0.8,paste("# initial:",nrow(initial))); text(0.5,0.5,paste("# kept: ",kept)); text(0.5,0.2,paste("# out: ",out))
	}

	train<-train[-i,]
}
keep.row<-rep(F,nrow(initial))

for(k in 1:nrow(initial)){
	if( sum(row.names(initial)[k]==keep)==1) keep.row[k]<-T
}
dev.off()

if(is.null(colvar))final<-initial[keep.row,colxy]
if(ncol(df)==2)final<-initial[keep.row,colxy]
if(!is.null(colvar)&ncol(df)>2)final<-initial[keep.row,c(colxy,colvar)]

if(plot==T){
	x11()
	plot(initial[,colxy],main="distribution of occurences",sub=paste("# initial (black):",nrow(initial)," | # kept (red): ",kept),pch=19,col="black",cex=0.2)
	points(final[,colxy],pch=19,col="red",cex=0.2)
}
return(final)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##
## add environmental values to a species dataframe.
## the xy (lat/long) coordinates of the species occurrences are compared to those of the environment dataframe
## and the value of the closest pixel is added to the species dataframe. 
## when the closest environment pixel is more distant than resolution, NA is added instead of the value.
## (similar to sample() in ArcGIS)

##ARGUMENTS
##dfsp: species dataframe with x, y and optional other variables
##colspxy: the range of columns for x and y in dfsp
##colspkept: the columns of dfsp that should be kept in the final dataframe (by default: xy )
##dfvar: environmental dataframe with x, y and environmental variables
##colvarxy: the range of columns for x and y in dfvar
##colvar: the range of enviromental variables columns in dfvar. (by default: all exept xy )
##resolution: distance between x,y of species and environmental datafreme after which values shouldn't be added 
##(typically, the resolution of the data in dfvar)

sample.sp.globvar <-function(dfsp,colspxy,colspkept="xy",dfvar,colvarxy,colvar="all",resolution){

if(sum(colspkept=="xy")==1)colspkept<-colspxy
if(sum(colvar=="all")==1) {
	if(!is.null(colspkept)) colvar<-(1:ncol(dfvar))[-colvarxy]
	if(is.null(colspkept))	colvar<-(1:ncol(dfvar))
}
colspx<-colspxy[1];colspy<-colspxy[2];colvarx<-colvarxy[1];colvary<-colvarxy[2]

x<-dfsp[,colspx]
X<-dfvar[,colvarx]
y<-dfsp[,colspy]
Y<-dfvar[,colvary]

train<-data.frame(matrix(nrow=nrow(dfsp),ncol=length(colvar)))
names(train)<-names(dfvar)[colvar]

#x11(2,2,pointsize = 12); par(mar=c(0,0,0,0));
for (i in 1:nrow(dfsp)){
	dist<-sqrt((X-x[i])^2 + (Y-y[i])^2)
	min<-min(dist)
	if(min<=resolution){
		if(length(colvar)>1)train[i,]<-dfvar[dist==min,colvar][1,]
		if(length(colvar)==1) train[i,]<-dfvar[dist==min,colvar][1]
	}
	#plot.new(); text(0.5,0.5,paste(paste("sampling:","\n","runs to go: ",nrow(dfsp)-i))); 
}
#dev.off()

if(!is.null(colspkept))final<-cbind(dfsp[,colspkept],train)
if(is.null(colspkept))final<-train

return(final)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##Investigate spatial autocorrelation by drawing a mantel Correlogram (autocorrelation vs distance)
##
##ARGUMENTS
##df: dataframe with x, y, and variables
##colxy: the range of columns for x and y in df
##colvar: the range of columns for variables in df
##n: number of random occurences used for the test (computation time increase tremendiously when using more than 500occ.) 	
##max: maximum distance to be computed in the correlogram
##nclass: number of class of distance to be computed in the correlogram
##nperm: number of permutation in the randomization process


mantel.correlogram <- function(df,colxy,n,colvar,max,nclass,nperm){

library(ecodist)

envnorm<-data.frame(t((t(df[,colvar])-mean(df[,colvar]))/sd(df[,colvar])))
row.rand<-sample(1:nrow(df),n,replace=T)
envdist<-dist(envnorm[row.rand,])
geodist<-dist(df[row.rand,colxy])
b<- seq(from = min(geodist), to = max, length.out = nclass)
crlg<-mgram(envdist,geodist,breaks=b,nperm=nperm)
plot(crlg)
abline(h=0)
}

##################################################################################################
##written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
##October 09. University of Lausanne. Switzerland
##
##DESCRIPTION
##randomly sample pseudoabsences from an environmental dataframe covering the study area
##A minimum distance from presences can be set.
##ARGUMENTS
##nbabsences: number of pseudoabsences desired 
##glob: environmental dataframe covering the study area to sample, with x,y 
##colxyglob: the range of columns for x and y in glob
##colvar: the range of columns for x and y in glob. colvar="all" keeps all the variables in glob in the final dataframe. colvar=NULL keeps only x and y
##presence: occurence dataframe 
##colxypresence: the range of columns for x and y in presence
##mindist: minimum distance from prensences closer to wich pseudoabsences shouldn't be drawn (buffer distance around presences)

rand.pseudoabsences<-function(nbabsences, glob, colxyglob,colvar="all", presence, colxypresence, mindist){

colxglob<-colxyglob[1]
colyglob<-colxyglob[2]
colxpresence<-colxypresence[1]
colypresence<-colxypresence[2]

keep<-c()

no.i<-1
while(no.i <= nbabsences){
	ki<-sample(1:nrow(glob),1)
	if(sum(((glob[ki,colxglob]- presence[,colxpresence])^2 + (glob[ki,colyglob]- presence[,colypresence])^2) <= mindist^2)==0) {
		keep[no.i]<-ki
 		no.i<-no.i+1
	}
}
if(sum(colvar=="all")==1) colvar<-(1:ncol(glob))[-colvarxy]
if(!is.null(colvar))pseudoabs<-glob[keep,c(colxyglob,colvar)]
if(is.null(colvar))pseudoabs<-glob[keep,colxyglob]

return(pseudoabs)
}




niche.equivalency.test<-function(z1,z2,rep){

  R<-length(z1$x)
  l<-list()
  
  obs.o<-niche.overlap(z1,z2,cor=T)									#observed niche overlap
  sim.o<-data.frame(matrix(nrow=rep.master,ncol=2))							#empty list of random niche overlap
  names(sim.o)<-c("D","I")
  total<-rep.master

  for (i in 1:total){
    if(is.null(z1$y)){ #overlap on one axis
      
      occ1.sim<-sample(z1$x,size=nrow(z1$sp),replace=T,prob=z1$z.cor) 		#random sampling of occurrences following the corrected densities distribution
      occ2.sim<-sample(z2$x,size=nrow(z2$sp),replace=T,prob=z2$z.cor)
      
      occ.pool<-c(occ1.sim,occ2.sim)   # pool of random occurrences
      rand.row<-sample(1:length(occ.pool),length(occ1.sim),replace=T) 		# random reallocation of occurrences to datasets
      sp1.sim<-occ.pool[rand.row]
      sp2.sim<-occ.pool[-rand.row]

      glob<-z1$glob
      sp<-data.frame(sp1.sim)
      z1.sim <- grid.clim(glob,sp,R)

      glob <- z2$glob
      sp <- data.frame(sp2.sim)
      z2.sim <- grid.clim(glob,sp,R)
    }
    
    if(!is.null(z1$y)){ 										#overlap on two axes
      coordinates<-which(z1$z.cor>0,arr.ind=T)						# array of cell coordinates ((1,1),(1,2)...)
      weight<-z1$z.cor[z1$z.cor>0] #densities in the same format as cells
      #print(sample(1:nrow(coordinates),size=nrow(z1$sp),replace=T))
      coordinates.sim1<-coordinates[sample(1:nrow(coordinates),size=nrow(z1$sp),replace=T,prob=weight),] # #random sampling of coordinates following z1$z.cor distribution
      occ1.sim<-cbind(z1$x[coordinates.sim1[,1]],z1$y[coordinates.sim1[,2]]) 	# random occurrences following the corrected densities distribution
      
      coordinates<-which(z2$z.cor>0,arr.ind=T)						# array of cell coordinates ((1,1),(1,2)...)
      weight<-z2$z.cor[z2$z.cor>0] 								#densities in the same format as cells
      coordinates.sim2<-coordinates[sample(1:nrow(coordinates),size=nrow(z2$sp),replace=T,prob=weight),] #random sampling of coordinates following z1$z.cor distribution
      occ2.sim<-cbind(z2$x[coordinates.sim2[,1]],z2$y[coordinates.sim2[,2]]) 	# random occurrences following the corrected densities distribution
      
      occ.pool<-rbind(occ1.sim,occ2.sim) # pool of random occurrences
      rand.row<-sample(1:nrow(occ.pool),nrow(occ1.sim),replace=T) 			# random reallocation of occurrences to datasets
      sp1.sim<-occ.pool[rand.row,]
      sp2.sim<-occ.pool[-rand.row,]
      
      glob<-z1$glob
      sp<-data.frame(sp1.sim)
      z1.sim <- grid.clim(glob,sp,R)
      
      glob<-z2$glob
      sp<-data.frame(sp2.sim)
      z2.sim <- grid.clim(glob,sp,R)
    }
    
    o.i<-niche.overlap(z1.sim,z2.sim,cor=F)							# overlap between random and observed niches
    sim.o$D[i]<-o.i$D											# storage of overlaps
    sim.o$I[i]<-o.i$I
  
  }

  sim.o<-na.omit(sim.o) 
  obs.o<-na.omit(obs.o)
  l$sim<-sim.o	# storage
  l$obs<-obs.o	# storage
  l$p.D<- min((sum(sim.o$D <= obs.o$D ) + 1),(sum(sim.o$D >= obs.o$D ) + 1))  *2/(length(sim.o$D) + 1) 	# storage of p-values
  if((sum(sim.o$D <= obs.o$D ) + 1) <  (sum(sim.o$D >= obs.o$D ) + 1)){ # (-) means that the obs D is significantly less than expected, positive means significantly more than expected
    l$p.D <- l$p.D * -1
  }
  
  l$p.I<- min((sum(sim.o$I <= obs.o$I ) + 1),(sum(sim.o$I >= obs.o$I ) + 1))  *2/(length(sim.o$I) + 1)	# storage of p-values
  
  return(l)
}





grid.clim<-function(glob,sp,R){

  # glob: global background dataset for the whole study area, 
  # sp: occurrence dataset
  # R: resolution of the grid
  l<-list()
  
  
  if(ncol(glob)==2){ #if scores in two dimensions (e.g. PCA)					
    
    spr<-data.frame(cbind((sp[,1]-xmin)/abs(xmax-xmin),(sp[,2]-ymin)/abs(ymax-ymin))) 			# makes the subset spcies occurrence climates between 0 and 1, relative to min and max of global climate

    # using a gaussian kernel density function, with RxR bins.
    spr<-spr[complete.cases(spr[,1:2]),]
    while(nrow(spr) < 5){
      spr<-rbind(spr,spr)
      print("we're multiplying!")
    }
    
    sp.dens<-kernelUD(SpatialPoints(spr[,1:2]),grid=mask,kern="bivnorm")[1]					# calculate the density of occurrences in a grid of RxR pixels along the score gradients
    # using a gaussian kernel density function, with RxR bins.
    
    sp.dens<-sp.dens[1]


    z<-matrix(sp.dens$ud*nrow(sp)/sum(sp.dens$ud),nrow=R,ncol=R,byrow=F) 			#rescale density to the number of occurrences in sp
    Z<-matrix(glob.dens$ud*nrow(glob)/sum(glob.dens$ud),nrow=R,ncol=R,byrow=F) 	#rescale density to the number of sites in glob
    
    z[z<max(z)/1000]<-0 											# remove infinitesimally small number generated by kernel density function
    Z[Z<max(Z)/1000]<-0 											# remove infinitesimally small number generated by kernel density function
    
    z.uncor<-z/max(z)											# rescale between [0:1] for comparison with other species	
    z<-z/Z												# correct for environment prevalence
    z[is.na(z)]<-0											# remove n/0 situations
    z[z=="Inf"]<-0											# remove n/0 situations
    z.cor<-z/max(z)											# rescale between [0:1] for comparison with other species	
    l$x<-x;l$y<-y;l$z.uncor<-z.uncor;l$z.cor<-z.cor;l$Z<-Z;l$glob<-glob;l$sp<-sp;
  }
  return(l)
}

niche.overlap<-function(z1,z2,cor){ 
  
  # z1 = species 1 occurrence density grid created by grid.clim
  # z2 = species 2 occurrence density grid created by grid.clim
  # cor=T correct occurrence densities of each species by the prevalence of the environments in their range
  
  l<-list()
  
  if(cor==F){		
    p1<-z1$z.uncor/sum(z1$z.uncor) # rescale occurence densities so that the sum of densities is the same for both species
    p2<-z2$z.uncor/sum(z2$z.uncor) # rescale occurence densities so that the sum of densities is the same for both species
  }
  
  if(cor==T){
    p1<-z1$z.cor/sum(z1$z.cor)	# rescale occurence densities so that the sum of densities is the same for both species
    p2<-z2$z.cor/sum(z2$z.cor)	# rescale occurence densities so that the sum of densities is the same for both species
  }
  
  D <-1-(0.5*(sum(abs(p1-p2))))				# overlap metric D
  I <-1-(0.5*(sqrt(sum((sqrt(p1)-sqrt(p2))^2))))	# overlap metric I
  l$D<-D
  l$I<-I
  return(l)
}