# PCA of PRISM precip to define regions
# MAC 06/3020

library(raster)
library(RStoolbox)
library(lubridate)

# set rasteroptions
rasterOptions(progress = 'text')

# functions
leap_every_year <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}

# load PRISM for SW Region
prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2020_PRISM_daily_prcp.grd") 

# crop to SW/NMX region
#e <- extent(-115,-102,25.5, 37.5)
# AZ and western NM
e <- extent(-115.5,-106,31.3, 37.5)
# AZ only
#e <- extent(subset(states, NAME_1=="Arizona"))
#e<-extent(aznm)
prcp <- (crop(prcp, e)) # can also add in rotate to convert lon to neg

# dates - find and remove leap days
startYr<-1981
  dates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
  colnames(dates)<-"date"
  dates$month<-as.numeric(format(dates$date, "%m"))
  dates$day<-as.numeric(format(dates$date, "%d"))
  dates$year<-as.numeric(format(dates$date, "%Y"))
  dates$doy<-as.numeric(format(dates$date, "%j"))
  dates$doy_ly<-leap_every_year(dates$date) # this day of year without leap day shifts

# subset layers to months of interest
#mos<-c(6,7,8,9)
mos<-c(7,8,9)
  subDates<-dates[which(dates$month %in% mos),]
  subLayers<-prcp[[which(dates$month %in% mos)]]
  # crop to region
  subLayers<-crop(subLayers,e)

# PC analysis  
pcs<-rasterPCA(subLayers, nSamples=10000, nComp=6)  

# plot all days in a year
library(rasterVis)
# map layers
states <- getData('GADM', country='United States', level=1)
az<-subset(states, NAME_1=="Arizona")
nm<-subset(states, NAME_1=="New Mexico")
aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico")
at<-c(seq(-350,350,25))
mapTheme <- rasterTheme(region = c("blue", "white","red"))
levelplot(pcs$map, contour=FALSE, at=at, 
          margin=FALSE, par.settings=mapTheme, 
          main=paste0("PCA Precip, JAS PRISM-daily 1981-2019"))+
  layer(sp.polygons(aznm, col = 'gray40', lwd=1))

plot(pcs$model$center)
plot(pcs$model$sdev[1:100])
plot(pcs$model$loadings[1,1:20])
screeplot(pcs$model, npcs=6)


# percent rank of all JAS days for each grid cell  
  perc.rank<-function(x) trunc(rank(x,ties.method = "average"))/length(x)
  percRankPrecip <- calc(subLayers, fun=perc.rank)
  percRankPrecip <-(percRankPrecip[[nlayers(percRankPrecip)]])*100
  
  writeRaster(percRankPrecip, filename = "/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2020_PRISM_daily_prcp.grd",
              overwrite=TRUE)
  
  # create NA MASK
  # load PRISM for SW Region
  prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2019_PRISM_daily_prcp.grd") 
  perc<- stack("/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2019_PRISM_daily_prcp.grd") 
    prcp<-prcp[[1]]
    perc<-perc[[1]]
   # e <- extent(-115.5,-106,31.3, 37.5) # AZ and wNM
    e <- extent(-114.85,-103,31.3, 37) # all AZ and NM
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
  
  
  # pentad median
  # # calculate pentads
  #  pentadMean<-function(x){
  #    movingFun(x, 5, mean, type = "around", na.rm = TRUE)
  #  }
    # pentadMedian<-function(x){
    #   movingFun(x, 5, median, type = "around", na.rm = TRUE)
    # }

    # beginCluster(6)
    #   prcp <- clusterR(prcp, calc, args=list(fun=pentadMedian))
    # endCluster()
  
    beginCluster(6)
       dailyClimo <- raster::clusterR(subLayers, stackApply,
                                          args=list(indices = subDates$doy_ly,
                                                    fun = mean, na.rm = TRUE))
    endCluster()
    doyDates<-seq.Date(as.Date("2016-07-01"),as.Date("2016-09-30"),1)
    names(dailyClimo)<-format(doyDates, "%b-%d")
    dailyClimo[dailyClimo == 0] <- NA  
    
    library(rasterVis)
    # map layers
    states <- getData('GADM', country='United States', level=1)
    az<-subset(states, NAME_1=="Arizona")
    nm<-subset(states, NAME_1=="New Mexico")
    aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico")
    
    at<-c(seq(0,20,1))
    #mapTheme <- rasterTheme(region=rev(terrain_hcl(12)))
    mapTheme <- rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
    pPrecip<-levelplot(dailyClimo[[which(doyDates>="2016-07-01" & doyDates<="2016-09-30")]],
                       contour=FALSE, margin=FALSE, par.settings=mapTheme, at=at,
                       main="Daily Median Precipitation PRISM, 1981-2019")+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    
# get monthly PRISM from RCC ACIS
    library(RCurl)
    library(jsonlite)
    library(raster)    
    # create current date
    #dateRangeStart=paste0(year,"-06-15")
    #dateRangeEnd= paste0(year,"-09-30")
    dateRangeStart="1900-01-01"
    dateRangeEnd= "2021-12-31"
  
    perc.rank<-function(x) trunc(rank(x,ties.method = "average"))/length(x)
    
    # generate dates -- keep with PRISM date
    allDates<-seq(as.Date(dateRangeStart), as.Date(dateRangeEnd),by="month")
    
    # AZ/NM bbox -115.004883,31.184609,-102.524414,37.387617
    ACISbbox<-"-115.5,31.3,-106,37.5" 
    
    # ACIS query
    jsonQuery=paste0('{"bbox":"',ACISbbox,'","sdate":"',dateRangeStart,'","edate":"',dateRangeEnd,'","grid":"21","elems":"mly_pcpn","meta":"ll,elev","output":"json"}') # or uid
    #jsonQuery=paste0('{"bbox":"',ACISbbox,'","sdate":"',dateRangeStart,'","edate":"',dateRangeEnd,'","grid":"2","elems":"pcpn","meta":"ll","output":"json"}') # or uid
    
    out<-postForm("http://data.rcc-acis.org/GridData",
                  .opts = list(postfields = jsonQuery,
                               httpheader = c('Content-Type' = 'application/json', Accept = 'application/json')))
    out<-fromJSON(out)
    
    # convert to list of matrices, flipud with PRISM
    matrixList <- vector("list",length(out$data))
    for(i in 1:length(out$data)){
      matrixList[[i]]<-apply(t(out$data[[i]][[2]]),1,rev)
    }
    
    # read into raster stack
    rasterList<-lapply(matrixList, raster)
    gridStack<-stack(rasterList)
    gridExtent<-extent(min(out$meta$lon), max(out$meta$lon), min(out$meta$lat), max(out$meta$lat))
    gridStack<-setExtent(gridStack, gridExtent, keepres=FALSE, snap=FALSE)
    names(gridStack)<-allDates
    # set 0 and neg to NA
    gridStack[gridStack < 0] <- NA
    ##
    allDates<-as.data.frame(allDates)
       allDates$month<-as.numeric(format(allDates$allDates, "%m"))
       allDates$year<-as.numeric(format(allDates$allDates, "%Y"))
    idx<-which(allDates$month %in% c(7,8,9)) # grab only summer months
    #idx<-which(allDates$month %in% c(7)) # grab only summer months
    allDates<-allDates[idx,]
    gridStack<-gridStack[[idx]]
    sumSeas<-stackApply(gridStack, allDates$year, fun = sum)
      seasAvgPrecip<-cellStats(sumSeas, 'mean')
      seasAvgPrecip<-cbind.data.frame(unique(allDates$year),seasAvgPrecip)
      seasAvgPrecip$percRank<-perc.rank(seasAvgPrecip$seasAvgPrecip) 
      colnames(seasAvgPrecip)<-c("year","avgPrecip","percRank")
      # names
      seasAvgPrecip$anomName<-"normal"
      seasAvgPrecip$anomName[seasAvgPrecip$percRank<=0.33] <- "dry"
      seasAvgPrecip$anomName[seasAvgPrecip$percRank>=0.66] <- "wet"
      
      library(cowplot)   
      ggplot(seasAvgPrecip, aes(year,avgPrecip, fill=as.factor(seasAvgPrecip$anomName)) )+
        geom_bar(stat = 'identity')+
        ggtitle("Regional Average Total Precip (July-Aug-Sept)")+
        geom_hline(yintercept=mean(seasAvgPrecip$avgPrecip), color="black")+
        geom_hline(yintercept=median(seasAvgPrecip$avgPrecip), color="red")+
        scale_fill_manual(values = c("saddlebrown", "grey", "forestgreen"), name="tercile")+
        ylab("inches")
    
  