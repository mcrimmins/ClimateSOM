# PCA of PRISM precip to define regions
# MAC 06/3020

library(raster)
library(RStoolbox)

# set rasteroptions
rasterOptions(progress = 'text')

# load PRISM for SW Region
prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2019_PRISM_daily_prcp.grd") 

# crop to SW/NMX region
#e <- extent(-115,-102,25.5, 37.5)
# AZ and western NM
e <- extent(-115.5,-106,31.3, 37.5)
# AZ only
#e <- extent(subset(states, NAME_1=="Arizona"))
#e<-extent(aznm)
#prcp <- (crop(prcp, e)) # can also add in rotate to convert lon to neg

# dates - find and remove leap days
startYr<-1981
  dates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
  colnames(dates)<-"date"
  dates$month<-as.numeric(format(dates$date, "%m"))
  dates$day<-as.numeric(format(dates$date, "%d"))
  dates$year<-as.numeric(format(dates$date, "%Y"))
  dates$doy<-as.numeric(format(dates$date, "%j"))
  #dates$doy_ly<-leap_every_year(dates$date) # this day of year without leap day shifts

# subset layers to months of interest
#mos<-c(6,7,8,9)
mos<-c(7,8,9)
  subDates<-dates[which(dates$month %in% mos),]
  subLayers<-prcp[[which(dates$month %in% mos)]]
  # crop to region
  subLayers<-crop(subLayers,e)

#pcs<-rasterPCA(subLayers)  

# percent rank of all JAS days for each grid cell  
  perc.rank<-function(x) trunc(rank(x,ties.method = "average"))/length(x)
  percRankPrecip <- calc(subLayers, fun=perc.rank)
  percRankPrecip <-(percRankPrecip[[nlayers(percRankPrecip)]])*100
  
  
  
  writeRaster(percRankPrecip, filename = "/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2019_PRISM_daily_prcp.grd",
              overwrite=TRUE)
  
  # create NA MASK
  # load PRISM for SW Region
  prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2019_PRISM_daily_prcp.grd") 
  perc<- stack("/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2019_PRISM_daily_prcp.grd") 
    prcp<-prcp[[1]]
    perc<-perc[[1]]
    e <- extent(-115.5,-106,31.3, 37.5)
      prcp<-crop(prcp,e)
      perc<-crop(perc,e)
      prcp[prcp >= 0] <- 1
      # test of mask
      mr <- mask(perc, prcp)
      # write out mask
      writeRaster(prcp, filename = "~/RProjects/SOMs/monsoonPrecip/SWUS_PRISM_MASK.grd",
                  overwrite=TRUE)
      
  
  
# plot all days in a year
  library(rasterVis)
  subLayers[subLayers == 0] <- NA  
  year<-2019
  at<-c(seq(0,50,2),80)
  mapTheme <- rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
  levelplot(subLayers[[which(subDates$year==year)]], contour=FALSE, 
            margin=FALSE, par.settings=mapTheme, 
            main=paste0(year," Precip  - PRISM-daily 1981-2019"))+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  
  