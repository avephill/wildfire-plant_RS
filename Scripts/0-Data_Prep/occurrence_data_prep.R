### This preps FIA occurrence data for analysis 0-primary_analysis.R, by cropping occurrence
### points to area of study and producing occurrence datasets in easily workable formats

### First download statewide FIA information from 
###     https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html
### and place in Data/Occurrence/Original/StateswFullInfo/

# All trees (standing live and dead) with a diameter at breast height (dbh) of at least 12.7 cm, 
# are inventoried on forested subplots. Within each subplot, a 2.07 m radius microplot offset 3.66 m 
# from subplot center is established where only live trees with a dbh 
# between 2.5 and 12.7 cm are inventoried. 
# Within each microplot, all live tree seedlings are tallied according to species.

# page 41 of core_ver7 info pdf for FIA for more info
# DSTR code 30= fire (crown or ground), 31=ground, 32=crown. Must have been within 5 years but there is date associated
# using disturbance requires significant threshold damage-- mortality/damage to 25% of all trees or 50% of an ind. species scount
# Saplings: Live trees 1.0 to 4.9 inches (2.5-12.5 cm) in diameter (DBH/DRC).
# Seedlings: Live trees smaller than 1.0 inch (2.5 cm) DBH/DRC that are at least 6 inches (15.2 cm) 
# in height for softwoods and 12-inches (30.5 cm) in height for hardwoods.

# To compare species in different life stages (off- spring vs. adult), we followed the FIA definition, 
# dividing the data into two types of subgroups, (i) seedling (dbh < 2.5 cm) vs. tree (dbh > 2.5 cm),
# and (ii) sapling (2.5 cm < dbh < 12.7 cm) vs. large tree (dbh > 12.7 cm).
library(tidyverse)
library(raster)
library(gridExtra)
library(tcltk2)
library(rgdal)
if(!exists("project_directory")) project_directory <- "" # Enter location of project directory within quotation marks
setwd(project_directory)

setwd("Data/Occurrence/Original")
spcd<-read.csv("species_codes.csv")

states <- c("CA", "OR", "WA", "ID", "MT", "WY", "UT", "NM", "SD", "CO")
# states <- "SD"
pl <- data.frame()
sd <- data.frame()
tr <- data.frame()
cond <- data.frame()

for(state in states){
  state_pl <-   read.csv(paste0("States/",state,"/",state,"_PLOT.csv"))
  state_sd <-   read.csv(paste0("States/",state,"/",state,"_SEEDLING.csv"))
  state_tr <-   read.csv(paste0("States/",state,"/",state,"_TREE.csv"))
  state_cond <- read.csv(paste0("States/",state,"/",state,"_COND.csv")) %>% 
    dplyr::select(PLT_CN, DSTRBCD1, DSTRBCD2, DSTRBCD3)
  
  pl  <- rbind(pl, state_pl)
  sd  <- rbind(sd, state_sd)
  tr  <- rbind(tr, state_tr)
  cond <- rbind(cond, state_cond)
}



# Only keep plots that have Level of Sampling Detail that effectively makes this plot a full tally of trees and seedlings

    # LOD = 1 : Data collected for vegetation structure only; total aerial canopy cover and canopy cover by     layer for tally tree species (all sizes), non-tally tree species (all size), shrubs/subshrubs/woody         vines, forbs, and graminoids
    
    # LOD = 2 : Vegetation structure data (LOD = 1) plus understory species composition data collected          including up to four most abundant species per GROWTH_HABIT_CD per subplot of: seedlings and saplings of     any tree species (tally or non-tally) <5 inches d.b.h. (d.r.c. for woodland species), shrubs/subshrubs      /woody vines, forbs, and graminoids.
    
    # LOD = 3 : Vegetation structure data, understory species composition data (LOD = 2), plus up to four       most abundant tree species (tally or non-tally)  5 inches d.b.h. (d.r.c for woodland species) per           GROWTH_HABIT_CD per subplot.

pl.filt<-pl[pl$P2VEG_SAMPLING_LEVEL_DETAIL_CD %in% c(2,3),]

# Here we determine whether a site burned by looking at the 'Disturbance' column of the
# Condition table, and indicate the result in our pl.filt data.frame

pl.postfilt <- pl.filt %>% 
  rename(PLT_CN = CN) %>% # rename pl ID column so it's same as condition table
  left_join(cond %>% filter(DSTRBCD1 %in% c(30,31,32) | # join condition table after making a fire column
                            DSTRBCD2 %in% c(30,31,32) |
                            DSTRBCD3 %in% c(30,31,32)) %>% 
              mutate(FIRE = "yes"), 
            by = "PLT_CN") %>% 
  rename(CN = PLT_CN)


setwd("../Processed")
write.csv(pl.postfilt,"pl_fire.csv")



# From here we decrease the size of the FIA Tree and Seed tables, and replace species code with binomial

sp_lite <- data.frame(SPCD=spcd$SPCD,SCI_NAME=spcd$SCI_NAME)
tr_sp <- merge(tr,sp_lite,by="SPCD",all=T) # translate spcd into binomial
tr_sp <- tr_sp[!is.na(tr_sp$PLT_CN),] # remove entries without plt number
tr_lite <- data.frame(CN=tr_sp$CN,PLT_CN=tr_sp$PLT_CN,STATECD=tr_sp$STATECD,
                   BINOM=tr_sp$SCI_NAME,DIA=tr_sp$DIA,TOTAGE=tr_sp$TOTAGE,
                   INVYR=tr_sp$INVYR,SUBPLOT_CN=tr_sp$SUBP)

move_to_seed <- tr_lite[tr_lite$TOTAGE <= 5 & !is.na(tr_lite$TOTAGE),c("CN","PLT_CN","STATECD","BINOM","DIA","INVYR","SUBPLOT_CN")]# this catches seedlings in tree table
# removes rows from tr_lite that are actually seedlings and have been placed in move_to_seed
tr_lite<-tr_lite[-which(tr_lite$TOTAGE <= 5 & !is.na(tr_lite$TOTAGE)),c("CN","PLT_CN","STATECD","BINOM","DIA","INVYR","SUBPLOT_CN")]
tr_lite<-tr_lite[order(tr_lite$CN),] #now here is tree dataframe with names
# CN PLT_CN STATECD                 BINOM
# 707281  1     64      41           Alnus rubra
# 706915  2     64      41           Alnus rubra
# 525473  3     64      41 Pseudotsuga menziesii
write.csv(tr_lite,"tr_lite.csv")


sp_lite<-data.frame(SPCD=spcd$SPCD,SCI_NAME=spcd$SCI_NAME)
sdsp<-merge(sd,sp_lite,by="SPCD",all=T) 
sdsp<-sdsp[!is.na(sdsp$PLT_CN),]

pre_sd_lite <- data.frame(CN=sdsp$CN, PLT_CN=sdsp$PLT_CN, STATECD=sdsp$STATECD, 
                          BINOM=sdsp$SCI_NAME, INVYR=sdsp$INVYR, SUBPLOT_CN=sdsp$SUBP)
sd_lite <- rbind(pre_sd_lite,move_to_seed[names(move_to_seed) != "DIA"])
sd_lite <- sd_lite[order(sd_lite$CN),] #now here is seed dataframe with names
# CN PLT_CN STATECD                 BINOM
# 707281  1     64      41           Alnus rubra
# 706915  2     64      41           Alnus rubra
# 525473  3     64      41 Pseudotsuga menziesii
write.csv(sd_lite,"sd_lite.csv")




# Here we make a table of all plots and the difference between seedlings and mature tree species composition
setwd(paste0(project_directory,"Data/Occurrence/Processed"))
pl.filt<-read.csv("pl_fire.csv",stringsAsFactors = F)
tr_lite<-read.csv("tr_lite.csv",stringsAsFactors = F)
sd_lite<-read.csv("sd_lite.csv",stringsAsFactors = F)


seeds_novel <- c()
trees_novel <- c()
plts_dif <- data.frame(matrix(ncol=5))
names(plts_dif) <- c("PLT_CN","NVL_TREE","NVL_SEED","BOTH","FIRE")
pl_rows <- nrow(pl.filt)
pb <- tkProgressBar(title = "progress bar", min = 0,max = pl_rows, width = 300)

for(i in 1:pl_rows){ # Iterate through each plot
  # Get unique trees and seedling of each plot
  t <- as.character(na.omit(unique(tr_lite[tr_lite$PLT_CN==pl.filt[i,]$CN,]$BINOM)))
  s <- as.character(na.omit(unique(sd_lite[sd_lite$PLT_CN==pl.filt[i,]$CN,]$BINOM)))
  
  t <- t[order(t)]
  s <- s[order(s)]
  
  if(length(s)<1){
    
  }else{
    if(!identical(s,t)){ # There is a difference in seedling species and tree species
      
      inter <- intersect(s,t)
      
      s_novel_i <- s[!(s %in% inter) ]# seedlings not shared between the two
      t_novel_i <- t[!(t %in% inter)]
      
      seeds_novel<-c(seeds_novel,s_novel_i)
      trees_novel<-c(trees_novel,t_novel_i)
      
      # Make sure there're NAs where s_novel_i and t_novel_i aren't same length, and plt is repeated
      max.len <- max(length(s_novel_i), length(t_novel_i), length(inter)) 
      s_novel_i = c(s_novel_i, rep(NA, max.len - length(s_novel_i)))
      t_novel_i = c(t_novel_i, rep(NA, max.len - length(t_novel_i)))
      inter= c(inter, rep(NA, max.len - length(inter)))
      
      plts_to_add <- c()
      fire_vals <- c() # Does this plot say 'yes' in fire?
      
      for(j in 0:(max.len-1)){
        plts_to_add <- c(plts_to_add, pl.filt[i,]$CN)
        fire_vals <- c(fire_vals, pl.filt[i,]$FIRE)
      }
      
      plts_add <- data.frame(PLT_CN=plts_to_add, NVL_SEED=s_novel_i, 
                             NVL_TREE=t_novel_i, BOTH=inter, FIRE=fire_vals)
      plts_dif <- rbind(plts_dif, plts_add)
      
    }
    
    
  }
  Sys.sleep(0.1) #for progress bar
  setTkProgressBar(pb, i, label=paste( round(i/pl_rows*100, 0),"% done")) #for progress bar
}
close(pb)
plts_dif<- plts_dif[2:nrow(plts_dif),] #first row is NA
write.csv(plts_dif,"plot_differences.csv")







# Here we make a table of all plots and the difference 
# between young trees and old tree species composition. This is the same analysis as the section
# immediately prior, but at different life stages


sapling <- tr_lite[tr_lite$DIA <= 12.7,] # centimeters
mature <- tr_lite[tr_lite$DIA > 12.7,]

CN_list <- c(sapling$PLT_CN,mature$PLT_CN)
pl.filt_JG <- subset(pl.filt, pl.filt$CN %in% CN_list)

#this does all the work and separation
sapling_novel<-c()
mature_novel<-c()
plt_difs_JG<-data.frame(matrix(ncol=5))
names(plt_difs_JG)<-c("PLT_CN","NVL_TREE","NVL_SEED","BOTH","FIRE")
total = nrow(pl.filt_JG)
pb <- tkProgressBar(title = "progress bar", min = 0, max = total, width = 300) # for progress bar
for(i in 1:total){
  t<-as.character(na.omit(unique(mature[mature$PLT_CN==pl.filt_JG[i,]$CN,]$BINOM)))
  s<-as.character(na.omit(unique(sapling[sapling$PLT_CN==pl.filt_JG[i,]$CN,]$BINOM)))
  
  t<-t[order(t)]
  s<-s[order(s)]
  
  if(length(s)<1){ #there must be a seedling
    
  }else{
    if(!identical(s,t)){ #wow there is a difference in seedling species and tree species!
      
      inter<-intersect(s,t)
      
      novs<-s[!(s %in% inter) ]# seedlings not shared between the two
      novt<-t[!(t %in% inter)]
      
      sapling_novel<-c(sapling_novel,novs)
      mature_novel<-c(mature_novel,novt)
      
      # all this makes sure there're NAs where novs and novt aren't same length, and plt is repeated
      max.len = max(length(novs), length(novt), length(inter)) 
      novs = c(novs, rep(NA, max.len - length(novs)))
      novt = c(novt, rep(NA, max.len - length(novt)))
      inter= c(inter, rep(NA, max.len - length(inter)))
      plts_to_add_JG<-c()
      firevals<-c()#does this plot say 'yes' in fire?
      for(j in 0:(max.len-1)){
        plts_to_add_JG<-c(plts_to_add_JG,pl.filt_JG[i,]$CN)
        firevals<-c(firevals,pl.filt_JG[i,]$FIRE)
      }
      
      plts_add_JG<- data.frame(PLT_CN=plts_to_add_JG,NVL_SEED=novs,NVL_TREE=novt,BOTH=inter,FIRE=firevals)
      plt_difs_JG<-rbind(plt_difs_JG,plts_add_JG)
      
    }
    
    
  }
  Sys.sleep(0.1) #for progress bar
  setTkProgressBar(pb, i, label=paste( round(i/total*100, 0),"% done")) #for progress bar
}
close(pb)
plt_difs_JG <- plt_difs_JG[2:nrow(plt_difs_JG),] #first row is NA
write.csv(plt_difs_JG,"plot_differences_sapling_mature.csv")