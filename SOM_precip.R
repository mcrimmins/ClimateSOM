# SOM for gridded precipitation 
# adapted from SOM_NARR.R
# MAC 6/26/20

library(raster)
library(lubridate)
library(reshape2)
library(kohonen)
library(tidyr)
library(ggplot2)
library(PBSmapping)
library(dplyr)
library(cowplot)
library(scales)
library(grid)
library(RColorBrewer)

ptm <- proc.time()

# set rasteroptions
rasterOptions(progress = 'text')

# functions
leap_every_year <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}
perc.rank <- function(x) trunc(rank(x))/length(x)

source("/home/crimmins/RProjects/StationPlots/rle2_function.R")

# map layers
states <- getData('GADM', country='United States', level=1)
  az<-subset(states, NAME_1=="Arizona")
  nm<-subset(states, NAME_1=="New Mexico")
  #aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico")
  aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico"| NAME_1=="California" | NAME_1=="Nevada" | NAME_1=="Utah" | NAME_1=="Colorado"| NAME_1=="Texas")
us <- getData('GADM', country='United States', level=0)
mx <- getData('GADM', country='Mexico', level=0)
#cn <- getData('GADM', country='Canada', level=0)
ecoreg<-rgdal::readOGR(dsn="~/RProjects/SOMs/monsoonPrecip/shapes", layer="us_eco_l3")
  ecoreg <- spTransform(ecoreg,crs(states))
# ELEVATION grid
elev<-raster("~/RProjects/SOMs/monsoonPrecip/shapes/PRISM_us_dem_4km_asc.asc")

# HUC4
huc4<-rgdal::readOGR(dsn="~/RProjects/SOMs/monsoonPrecip/shapes", layer="huc4clip")
# station points
# stations<-cbind.data.frame(c("Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque"),
#                            c(32.1145,33.4373,35.1404,36.0840,31.8053,35.0433),
#                            c(-110.9392,-112.0078,-111.6690,-115.1537,-106.3824,-106.6129))
stations<-cbind.data.frame(c("Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque"),
                           c(32.13130,33.42770,35.14410,36.07190,31.81111,35.04190),
                           c(-110.9552,-112.0038,-111.6663,-115.1634,-106.3758,-106.6155))

colnames(stations)<-c("stations","latitude","longitude")
coordinates(stations)= ~ longitude+latitude
# station codes
stationCodes<-cbind.data.frame(c("TUS","PHX","FLG","LAS","ELP","ABQ"),
                               c(32.13130,33.42770,35.14410,36.07190,31.81111,35.04190),
                               c(-110.9552,-112.0038,-111.6663,-115.1634,-106.3758,-106.6155))


# load PRISM for SW Region - from dailyDownloadPRISM.R
# prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2019_PRISM_daily_prcp.grd") 
 prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2020_PRISM_daily_prcp.grd") 
# PRISM percentiles from monsoonPRISM_misc.R
# prcp<-stack("/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2019_PRISM_daily_prcp.grd") 
# load PRISM mask from monsoonPRISM_misc.R
mask<-raster( "~/RProjects/SOMs/monsoonPrecip/SWUS_PRISM_MASK.grd")


# load CPC daily precip from ~/SWMonsoonTracker/NCEPGrids/createNCEPgrids.R
#prcp<- rotate(stack("/scratch/crimmins/cpc_global/processed/CPC_Daily_precip_global_1979_2019_NAClip.grd")) 
#    names(prcp)<-seq.Date(as.Date("1979-01-01"),as.Date("2019-12-31"),1)
    
# # load pre-processed raster stacks from /ClimPlot/NARR/processNARR.R
# gh500<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2019.grd")
# pwat<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2019.grd") 
# prcp<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2019.grd") 
# avgCAPE<-stack("/scratch/crimmins/NARR/processed/CAPE_daily_NARR_WUS_1979_2019.grd")
# maxCAPE<-stack("/scratch/crimmins/NARR/processed/maxCAPE_daily_NARR_WUS_1979_2019.grd")
# sfcDP<-stack("/scratch/crimmins/NARR/processed/DPT2m_daily_NARR_WUS_1979_2019.grd")
# uFlux<-stack("/scratch/crimmins/NARR/processed/WVUFLX_daily_NARR_WUS_1979_2019.grd")
# vFlux<-stack("/scratch/crimmins/NARR/processed/WVVFLX_daily_NARR_WUS_1979_2019.grd")

# crop to SW/NMX region
#e <- extent(-115,-102,25.5, 37.5)
# AZ and western NM
    e <- extent(-115.5,-106,31.3, 37.5) #(-115.5,31.3,-106,37.5)
# AZ only
#e <- extent(subset(states, NAME_1=="Arizona"))
#e<-extent(aznm)
#prcp <- (crop(prcp, e)) # can also add in rotate to convert lon to neg
# AZ and all of NM
    #e <- extent(-114.85,-103,31.3, 37) #

# dates - find and remove leap days
startYr<-1981
#dates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
dates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2020-12-31"),1))
  colnames(dates)<-"date"
  dates$month<-as.numeric(format(dates$date, "%m"))
  dates$day<-as.numeric(format(dates$date, "%d"))
  dates$year<-as.numeric(format(dates$date, "%Y"))
  dates$doy<-as.numeric(format(dates$date, "%j"))
  dates$doy_ly<-leap_every_year(dates$date) # this day of year without leap day shifts
  
# # JAS -- subset layers to months of interest
#   #mos<-c(6,7,8,9)
#   mos<-c(7,8,9)
#   subDates<-dates[which(dates$month %in% mos),]
#   #subLayers<-prcp # FOR PERC GRIDS
#   subLayers<-prcp[[which(dates$month %in% mos)]]
#     # crop to region
#     subLayers<-crop(subLayers,e)
#     elev<-crop(elev, e)
#     # apply mask
#     subLayers <- mask(subLayers, mask)
  
# June 15th-Sept 30 subset layers to months of interest
  subDates<-dates[which(dates$doy_ly %in% seq(167,274,1)),]
  #subLayers<-prcp # FOR PERC GRIDS
  subLayers<-prcp[[which(dates$doy_ly %in% seq(167,274,1))]]
    # crop to region
    subLayers<-crop(subLayers,e)
    elev<-crop(elev, e)
    # apply mask
    subLayers <- mask(subLayers, mask)    
    
    
  #subLayers<-prcp # for already subsampled percentiles
  
  # summarize to month and season for anomalies
    subDates$sumDate<-as.Date(paste0(subDates$year,"-",subDates$month,"-01"),format="%Y-%m-%d")
    sumMonth<- stackApply(subLayers, subDates$sumDate, fun = sum)
    sumSeas<-stackApply(subLayers, subDates$year, fun = sum)
    moAvgPrecip<-cellStats(sumMonth, 'mean')
      moAvgPrecip<-cbind.data.frame(unique(subDates$sumDate),moAvgPrecip)
    seasAvgPrecip<-cellStats(sumSeas, 'mean')
      seasAvgPrecip<-cbind.data.frame(unique(subDates$year),seasAvgPrecip)
    seasAvgPrecip$percRank<-perc.rank(seasAvgPrecip$seasAvgPrecip) 
    colnames(seasAvgPrecip)<-c("year","avgPrecip","percRank")
    # names
    seasAvgPrecip$anomName<-"normal"
    seasAvgPrecip$anomName[seasAvgPrecip$percRank<=0.333334] <- "dry"
    seasAvgPrecip$anomName[seasAvgPrecip$percRank>0.67] <- "wet"
    
    # ACIS stations --- update using 
    #load("stationPrecip.RData")
    load("stationPrecip_thru2020.RData")
      stationDaily<-data
    stationSumSeas<-stationDaily %>% group_by(year) %>% summarise_at(c(2:7), sum, na.rm = TRUE) 
      colnames(stationSumSeas)<-c("year","Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque")
      stationSumSeas[,c(2:7)]<-stationSumSeas[,c(2:7)]*25.4
    stationDaily<-stationDaily[,c(-1,-8,-9)]
      colnames(stationDaily)<-c("Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque")
      stationDaily<-stationDaily*25.4
    
    # PRISM get station extracts
    # stationSumSeas<-as.data.frame(t(raster::extract(sumSeas, stations)))
    #  colnames(stationSumSeas)<-stations$stations
    # get daily station data  
    stationDailyPRISM<-as.data.frame(t(raster::extract(subLayers, stations)))
      colnames(stationDailyPRISM)<-stations$stations
    # get watershed stats
    #huc4SumSeas<-as.data.frame(t(raster::extract(sumSeas, huc4, fun=mean)))
    #  colnames(huc4SumSeas)<-huc4$NAME
    #huc4daily<-as.data.frame(t(raster::extract(subLayers, huc4, fun=mean)))
    #  colnames(huc4daily)<-huc4$NAME  
    
    # dry year differences
      # dry3<-sumSeas[[head(sort(seasAvgPrecip$percRank, index.return=TRUE)$ix,3)]]
      #   fun <- function(x) { combn(x, 2, function(x) x[1] - x[2]) }
      # diffs <- calc(dry3, fun)
      # temp<-as.vector(matrix(do.call(paste0, as.data.frame(t(combn(names(dry3),2))))))
      #   temp<-gsub("[a-zA-Z ]", "", temp)
      #   temp<-substring(temp,2)
      # names(diffs)<-temp
      #   plot(diffs, zlim=c(-300,300))
      
  # JAS convert layers to dataframe
  # layers.df<-(as.data.frame(subLayers, long=TRUE, xy=TRUE))
  # colnames(layers.df)<-c("lon","lat","date","value")  
  # # long to wide
  # df.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
  
  # # June 15-Sept 30 -- convert layers to dataframe
  # layers.df<-(as.data.frame(subLayers[[which(subDates$year<=2000)]], long=TRUE, xy=TRUE))
  # colnames(layers.df)<-c("lon","lat","date","value")  
  # # long to wide
  # df.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
  #   # second chunk
  #   gc()
  #   layers.df<-(as.data.frame(subLayers[[which(subDates$year>2000)]], long=TRUE, xy=TRUE))
  #   colnames(layers.df)<-c("lon","lat","date","value")  
  #   gc()
  #     # long to wide
  #   df.wide2<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
  #   rm(layers.df)
  #   # combine
  #   df.wide<-rbind.data.frame(df.wide,df.wide2)
  #   rm(df.wide2)
  # 
  # # save raw precip
  # save(df.wide, file="~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JJAS_thru2020.RData")
  
  load("~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JJAS_thru2020.RData")
  # save perc precip
  # save(df.wide, file="~/RProjects/SOMs/monsoonPrecip/perc_AZwNM_PRISM_JAS.RData")
  #  load("~/RProjects/SOMs/monsoonPrecip/perc_AZwNM_PRISM_JAS.RData")
 
  # replace NAs with 0's if needed
  df.wide[is.na(df.wide)] <- 0
  
  # mask out on precip threshold
  #df.wide[df.wide<=5] <- NA
  
  # run SOM with all days or subset based on specific node
  idx<-seq(1,nrow(df.wide),1)
  #idx<-which(som.gh500$unit.classif==which(codeList=="1_1"))
  #whichNodes<-c("1_1","1_2","1_3","2_1","2_2","2_3","3_1")
 # idx<-which(som.gh500$unit.classif %in% match(codeList,whichNodes))

  
# #  SOM screening #####
# source("topo.error.R")
# ncols=c(2,3,4,3,5,4,6,7,5,4,8,9,6,5,10,11,6,8,5,13,7,14)
# nrows=c(2,2,2,3,2,3,2,2,3,4,2,2,3,4,2,2,4,3,5,2,4,2)
# qe<-c()
# te<-c()
# ptm <- proc.time()
# for(i in 1:length(nrows)){
#   som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
#                         #grid = somgrid(ncols[i], nrows[i], "rectangular"),
#                         grid = somgrid(ncols[i], nrows[i], topo="rectangular", neighbourhood.fct = c("gaussian")),
#                         #alpha = c(0.05, 0.001),
#                         radius = c(4,1),
#                         mode= "pbatch",
#                         cores = 7,
#                         rlen = 500,
#                         dist.fcts = "sumofsquares")
#   print(paste0(nrows[i],"x",ncols[i],"-round1"))
#   som.gh500.2 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
#                           #grid = somgrid(ncols, nrows, "rectangular"),
#                           grid = somgrid(ncols[i], nrows[i], topo="rectangular", neighbourhood.fct = c("gaussian")),
#                           #alpha = c(0.05, 0.001), # for online
#                           radius = c(3,0.33),
#                           mode= "pbatch", # "pbatch" or "online"
#                           #mode= "online", # "pbatch" or "online"
#                           #maxNA.fraction = 0.5,
#                           init= som.gh500$codes,
#                           cores = 7,
#                           rlen = 500,
#                           dist.fcts = "sumofsquares")
#   ## quantization error:
#   qe[i]<-mean(som.gh500.2$distances)
#   ## topographical error measures:
#   te[i]<-topo.error(som.gh500.2, "nodedist")
#   print(paste0(nrows[i],"x",ncols[i],"-round2"))
# }
# proc.time() - ptm
# 
# diagnostics<-cbind.data.frame(nrows,ncols,qe,te)
# plot(diagnostics$qe, diagnostics$te, xlim=c())
# text(diagnostics$qe, diagnostics$te, labels=paste0(diagnostics$nrows,"-",diagnostics$ncols))
# save(diagnostics, file = "~/RProjects/SOMs/monsoonPrecip/diagnostics_CP2_JJAS_1981_2020.RData")
# ######

#####  
# kohonen SOM
  # nrows=4
  # ncols=4
  # ptm <- proc.time()
  #   #som.gh500 <- som(as.matrix(df.wide[,2:ncol(df.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
  #   #set.seed(999) #keep set seed 16, 11 upper left/wet, 9,6 upper right dry/UL wet, 8 UL Wet/LR dry, 7 LR Wet/LL Dry, 5 LR wet/LL dry, 4/100 UL wet/UR dry, 101 wet UR/dry LL, 102/104 dry UL/wet UR, 103 UL Wet/LR Dry
  #   set.seed(128) # 4x4 set seed to 128 for inactive upper left
  #   som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
  #                         #grid = somgrid(ncols, nrows, "rectangular"),
  #                         grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
  #                         #alpha = c(0.05, 0.001), # for online
  #                         radius = c(4,1), # c(3,0.33)
  #                         mode= "pbatch", # "pbatch" or "online"
  #                         #mode= "online", # "pbatch" or "online"
  #                         #maxNA.fraction = 0.5,
  #                         cores = 7,
  #                         rlen = 700,
  #                         dist.fcts = "sumofsquares")
  # proc.time() - ptm
  # ## quantization error:
  # mean(som.gh500$distances)
  # ## topographical error measures:
  # source("topo.error.R")
  # topo.error(som.gh500, "nodedist")
  # plot(som.gh500, type="counts", shape = "straight", labels=counts)
  # #####
  
# #   ##### 
#   # CP2 - 2 phase SOM training
#   nrows=4
#   ncols=4
#   ptm <- proc.time()
#   #som.gh500 <- som(as.matrix(df.wide[,2:ncol(df.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
#   #set.seed(999) #keep set seed 16, 11 upper left/wet, 9,6 upper right dry/UL wet, 8 UL Wet/LR dry, 7 LR Wet/LL Dry, 5 LR wet/LL dry, 4/100 UL wet/UR dry, 101 wet UR/dry LL, 102/104 dry UL/wet UR, 103 UL Wet/LR Dry
#   set.seed(128) # 123 for 4x4 JAS, 128 for 4x4 JJAS
#   som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
#                         #grid = somgrid(ncols, nrows, "rectangular"),
#                         grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
#                         #alpha = c(0.05, 0.001), # for online
#                         radius = c(4,1), #(4,1) for 4x5
#                         mode= "pbatch", # "pbatch" or "online"
#                         #mode= "online", # "pbatch" or "online"
#                         #maxNA.fraction = 0.999,
#                         cores = 7,
#                         rlen = 5000, #5000
#                         dist.fcts = "sumofsquares")
#   print("completed phase 1, performing phase 2")
#   som.gh500.2 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
#                         #grid = somgrid(ncols, nrows, "rectangular"),
#                         grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
#                         #alpha = c(0.05, 0.001), # for online
#                         radius = c(2,0.33), # c(3,0.33) for 3x5
#                         mode= "pbatch", # "pbatch" or "online"
#                         #mode= "online", # "pbatch" or "online"
#                         #maxNA.fraction = 0.999,
#                         init= som.gh500$codes,
#                         cores = 7,
#                         rlen = 15000, #7000
#                         dist.fcts = "sumofsquares")
#       ## quantization error:
#       mean(som.gh500.2$distances)
#       ## topographical error measures:
#       source("topo.error.R")
#       topo.error(som.gh500.2, "nodedist")
#   proc.time() - ptm
#   som.gh500<-som.gh500.2
#   #####
# 
#  # test<-calc(subLayers, sd)
# #  test2<-calc(subLayers, mean)
# 
#   # save SOM output
#    save(som.gh500, file = "~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JJAS_SOM4x4_15K_CP2_1981_2020.RData")
#   # plot(som.gh500, type="changes") # changes, codes, counts, property, quality, mapping
  
  
  load("~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JJAS_SOM4x4_15K_CP2_1981_2020.RData") #1- rad(4,1):rlen:5000, #2-rad(3,0.33):rlen:7000
  nrows=4
  ncols=4
  
  #topo.error(som.gh500, "bmu")
  codebook<-as.data.frame(som.gh500$codes)
  code_grid<-as.data.frame(round(som.gh500$grid$pts,0))
  code_grid$mapUnit<-seq(1,nrow(code_grid))
  code_grid<-code_grid %>%
    unite(y,x, col="codes", sep="_") # add mapunit back in if needed
    
    
  
  # deal with code book
  codebook<-cbind(code_grid,codebook)
  codebook.long<-melt(codebook, id.vars = 1)
  # separate out mapunits  ----
  mapunits<-codebook.long[1:(nrows*ncols),]
  codebook.long<-codebook.long[-(1:(nrows*ncols)),]
  # ####
  codebook.long<-separate(codebook.long, variable, convert = TRUE, into = c("lat", "lon"), sep="_.")
  codebook.long$lat<-as.numeric(gsub("X", "", codebook.long$lat))
  codebook.long$lon<--codebook.long$lon
  codebook.long<-separate(codebook.long, codes, convert = FALSE, remove = FALSE, into = c("yRows","xCols"), sep="_") # switch position of yRows/xCols
  
  # assign days to nodes
  #nodes<-kohonen::map(som.gh500)
  #somTime<-as.data.frame(cbind(df.wide$date, nodes$unit.classif, nodes$distances))
  somTime<-as.data.frame(cbind(df.wide$date[idx], som.gh500$unit.classif, som.gh500$distances))
  colnames(somTime)<-c("date","mapUnit","errorDist")
  somTime$date<-as.character(somTime$date)
  somTime$date<-as.Date(somTime$date, format="X%Y.%m.%d")
  somTime$month<-as.numeric(format(somTime$date, format="%m"))
  somTime$year <-as.numeric(format(somTime$date, format="%Y"))
  somTime$day  <-as.numeric(format(somTime$date, format="%d"))
  #somTime<-separate(somTime, as.character(date), convert = TRUE, remove = FALSE, into = c("year","month","day"), sep=".")
  #somTime$date<-as.Date(paste(somTime$year,"-",somTime$day,"-",somTime$month,sep=""), format="%Y-%d-%m")
  somTime$doy<-as.numeric(format(somTime$date, "%j"))
  somTime$mapUnit<-as.integer(as.character(somTime$mapUnit))
  somTime$errorDist<-as.numeric(as.character(somTime$errorDist))
  # join codes to node table
  somTime<-left_join(somTime, code_grid)
  # get codeList
  codeList<-unique(codebook.long$codes)
  
  # expand grid to double
  #som.gh500.2 <- expandMap(som.gh500)
  #plot(som.gh500.2, type="counts", shape = "straight", labels=counts)
  summary(som.gh500)

# ACTIVITY CLASSIFICATIONS    
# add in activity categories - 4x symmetric
# activity<- rbind.data.frame(cbind(c("1_1"), c("Inactive"),c("Inactive")),
#             cbind(c("1_2","2_1","2_2"),rep("Active-low",3),rep("Active",3)),
#             cbind(c("1_3","3_1","3_2","2_3","3_3"), rep("Active-high",5),rep("Active",5)),
#             cbind(c("1_4","4_1","4_2","2_4","3_4","4_3","4_4"), rep("Widespread"),rep("Active")))
# colnames(activity)<-c("codes","activityCat","activityCat2")
# somTime<-left_join(somTime,activity)
# alt activity based on kmeans majority
# activityAlt<- rbind.data.frame(cbind(c("1_1"), c("Inactive")),
#                             cbind(c("1_2","2_1","2_2"),rep("Active-low",3)),
#                             cbind(c("1_3","3_1","3_2","2_3","3_3","1_4","2_4","4_2"), rep("Active-high")),
#                             cbind(c("4_1","3_4","4_3","4_4"), rep("Widespread")))
# colnames(activityAlt)<-c("codes","activityCatAlt")
# somTime<-left_join(somTime,activityAlt)
# LARGER ACTIVE HIGH CATEGORY
activity<- rbind.data.frame(cbind(c("1_1"), c("Inactive"),c("Inactive")),
                            cbind(c("1_2","2_1","2_2"),rep("Active-low",3),rep("Active",3)),
                            cbind(c("1_3","3_1","3_2","2_3","3_3","1_4","4_1","4_2","2_4"), rep("Active-high",9),rep("Active",9)),
                            cbind(c("3_4","4_3","4_4"), rep("Widespread"),rep("Active")))
colnames(activity)<-c("codes","activityCat","activityCat2")
somTime<-left_join(somTime,activity)
######



    
  # RASTERVIS mapping of SOM results
  library(rasterVis)
  library(colorspace)
  # create stack of codebook results
  cbRaster<-stack()
  for (i in 1:length(codeList)) {
    temp<-codebook.long[which(codebook.long$codes==codeList[i]),c(5,4,6)]
    temp<-rasterFromXYZ(temp)
    cbRaster<-stack(cbRaster, temp)
  }
  names(cbRaster)<-codeList

  col.titles<-paste0("",codeList," (",round((table(somTime$codes)/nrow(somTime))*100,1),"%, ",
                     round(table(somTime$codes)/length(unique(somTime$year)),1),"d/yr)")
  # location of precip centroid
  #cbRaster[cbRaster ==0] <- NA
  # xyCentroid<-list()
  # for (i in 1:length(codeList)) {
  #   temp1<-cbRaster[[i]]
  #   xyCentroid[[i]]<-colMeans(xyFromCell(temp1, which(temp1[]>1)))
  # }
  # xyCentroid<-as.data.frame(do.call("rbind", xyCentroid))
  # coordinates(xyCentroid)<-~x + y
  # # location of precip max
  # xyCentroid<-list()
  # for (i in 1:length(codeList)) {
  #   temp1<-cbRaster[[i]]
  #   idxMax<-which.max(temp1)
  #   xyCentroid[[i]]<-xyFromCell(temp1,idxMax)
  # }
  # xyCentroid<-as.data.frame(do.call("rbind", xyCentroid))
  # coordinates(xyCentroid)<-~x + y
  # location of closest regional mean value
  cbRaster[cbRaster ==0] <- NA
  xyCentroid<-list()
  for (i in 1:length(codeList)) {
    temp1<-cbRaster[[i]]
    xyCentroid[[i]]<-colMeans(xyFromCell(temp1, which(temp1[]>cellStats(temp1,mean))))
  }
  xyCentroid<-as.data.frame(do.call("rbind", xyCentroid))
  coordinates(xyCentroid)<-~x + y
  
  
  #####
  ## location of max value
      maxIdx<-which.max(cbRaster)
      maxIdx<-as.factor(maxIdx)
        rat <- levels(maxIdx)[[1]]
        rat[["codes"]]<- codeList[rat$ID]
        levels(maxIdx)<- rat
        mapTheme <- rasterTheme(region=brewer.pal(7,"Set1"))
      levelplot(maxIdx, par.settings=mapTheme, main="Node contributing max codebook value")+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  #####
  
  # ecoregion with max value
  ecoregNodeMax <- raster::extract(cbRaster, ecoreg, fun=max, na.rm=TRUE, df=TRUE)
  ecoregNodeMax<-cbind.data.frame(ecoregNodeMax,ecoreg$US_L3NAME)
  ecoregNodeMax<-na.omit(ecoregNodeMax)
  ecoregNodeMax[which.max(ecoregNodeMax$X1_1),c("ecoreg$US_L3NAME")]
    ecoName<-cbind.data.frame(ecoregNodeMax[apply(ecoregNodeMax[,c(2:17)],MARGIN=2,FUN=function(x) {which.max(x)}),18],
                                codeList)    
    colnames(ecoName)<-c("ecoregion","codes")  
      
  # activity text for maps    
  # temp<-(unique(somTime[,c("mapUnit","activityCat","codes")]))
  #     temp<-temp[order(temp[,1]),]
  #     text2add<-as.character(temp[,2])
  #     # color of text
  #     textCol<-c("green","yellow","orange","red",
  #                "yellow","yellow","orange","red",
  #                "orange","orange","orange","red",
  #                "red","red","red","red")
  #     # add ecoName
  #     text2add<-paste0(text2add<-as.character(temp[,2]),"-",ecoName$ecoregion)
  
  # map descriptors text
      temp<-(unique(somTime[,c("mapUnit","activityCat","codes")]))
      temp<-temp[order(temp[,1]),]
      textName<-c(" ","AZ","SE AZ","S AZ/NM",
                  "SE AZ","NW NM","AZ","NW AZ",
                  "SW NM"," SE AZ/SW NM","N AZ/NM","AZ",
                  "NW NM", "SW NM","W NM","AZ-NM")
      text2add<-textName
      
      
  #cbRaster[cbRaster < 0.254] <- NA  
  #cbRaster[cbRaster ==0] <- NA  
  at<-c(seq(0,30,0.3))
  #at<-c(seq(0,30,2.5))
  mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  #at<-c(seq(0,10,0.25))
  #mapTheme <- rasterTheme(region=rev(terrain_hcl(12)))
  #at<-c(seq(0,1,0.05))
  #mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
  pPrecip<-levelplot(cbRaster, contour=FALSE, margin=FALSE, layout=c(ncols,nrows), at=at,
                     names.attr=col.titles,
                     par.settings=mapTheme,
                     #par.settings = list(region=c("lightblue", "blue","green","green4","yellow","red", "red4"),
                      #                   axis.line = list(col = textCol[panel.number()])),
                     scales=list(draw=FALSE),
            main="Precip Patterns Jun 15-Sep 30 4x4 SOM - PRISM-daily 1981-2020")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))+
    layer(sp.points(stations, col="black", pch = 18, cex=0.5))+
    contourplot(elev, at=c(2000), labels=FALSE, lwd = 0.3, par.settings = GrTheme)+
    #layer(sp.polygons(ecoreg, col="black", pch = 18, cex=0.5))+
    #layer(sp.text(c(36.0840,-115.1537),"LAS"))+
    #layer(sp.points(xyCentroid[panel.number()],
    #              pch=20, cex=1, col="black"))+
    #layer(panel.rect(-115.5,31.3,-113.0,32.0,fill=textCol[panel.number()]))+
    #layer(panel.rect(-115.5,31.3,-113.0,32.0, border=textCol[panel.number()]))+
    layer(panel.text(-114.25, 31.55, text2add[panel.number()],col="black",cex=0.6))
    #layer(panel.text(stationCodes[panel.number(1),3], stationCodes[panel.number(1),2],
    #                 stationCodes[panel.number(1),1],col="black",cex=0.4))
  
    #layer(sp.polygons(huc4, col = 'gray20', lwd=0.75))
    #layer(sp.polygons(us, col = 'gray40', lwd=1))+
    #layer(sp.polygons(mx, col = 'gray40', lwd=1))
    
    # add lines
    #grid.ls(viewport=TRUE, grobs=FALSE) #plot_01.toplevel.vp::plot_01.panel.1.1.vp
    #grid.rect(vp = "plot_01.toplevel.vp::plot_01.panel.1.1.vp",
    #          gp = gpar(col = "red"))
     png("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/precip_JJAS_SOM_4x4.png", width = 10, height = 6, units = "in", res = 300L)
     pPrecip
     #####
     # # vert line - original
     #  grid.lines(x = unit(c(0.2465, 0.2465), "npc"),
     #             y = unit(c(0.75, 1), "npc"),
     #             default.units = "npc",
     #             arrow = NULL, name = NULL,
     #             gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
     #  grid.lines(x = unit(c(0.5, 0.5), "npc"),
     #             y = unit(c(0.497, 1), "npc"),
     #             default.units = "npc",
     #             arrow = NULL, name = NULL,
     #             gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")    
     #  grid.lines(x = unit(c(0.753, 0.753), "npc"),
     #             y = unit(c(0.246, 1), "npc"),
     #             default.units = "npc",
     #             arrow = NULL, name = NULL,
     #             gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp") 
     #  # horiz line 
     #  grid.lines(x = unit(c(0, 0.2465), "npc"),
     #             y = unit(c(0.75,0.75), "npc"),
     #             default.units = "npc",
     #             arrow = NULL, name = NULL,
     #             gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")  
     #  grid.lines(x = unit(c(0, 0.5), "npc"),
     #             y = unit(c(0.497,0.497), "npc"),
     #             default.units = "npc",
     #             arrow = NULL, name = NULL,
     #             gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
     #  grid.lines(x = unit(c(0, 0.753), "npc"),
     #             y = unit(c(0.246,0.246), "npc"),
     #             default.units = "npc",
     #             arrow = NULL, name = NULL,
     #             gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
     #  grid.text( label = "Inactive", x = 0.11, y = 0.80, rot = 90, default.units = "npc", gp=gpar(col="red"))
     #  grid.text( label = "Active-Low", x = 0.11, y = 0.58, rot = 90, default.units = "npc", gp=gpar(col="red"))
     #  grid.text( label = "Active-High", x = 0.11, y = 0.36, rot = 90, default.units = "npc", gp=gpar(col="red"))
     #  grid.text( label = "Widespread", x = 0.11, y = 0.14, rot = 90, default.units = "npc", gp=gpar(col="red"))
     #  #print(pPrecip, newpage = FALSE)
     #####  
      # vert line - 3 widespread
      grid.lines(x = unit(c(0.2465, 0.2465), "npc"),
                 y = unit(c(0.75, 1), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
      grid.lines(x = unit(c(0.5, 0.5), "npc"),
                 y = unit(c(0.497, 1), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")    
      grid.lines(x = unit(c(0.753, 0.753), "npc"),
                 y = unit(c(0.246, 0.5), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp") 
      
      grid.lines(x = unit(c(0.5, 0.5), "npc"),
                 y = unit(c(0,0.246), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
      
      # horiz line 
      grid.lines(x = unit(c(0, 0.2465), "npc"),
                 y = unit(c(0.75,0.75), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")  
      grid.lines(x = unit(c(0, 0.5), "npc"),
                 y = unit(c(0.497,0.497), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
      grid.lines(x = unit(c(0.5, 0.753), "npc"),
                 y = unit(c(0.246,0.246), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
      
      grid.lines(x = unit(c(0.753,1), "npc"),
                 y = unit(c(0.5,0.5), "npc"),
                 default.units = "npc",
                 arrow = NULL, name = NULL,
                 gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
      
      grid.text( label = "Inactive", x = 0.11, y = 0.80, rot = 90, default.units = "npc", gp=gpar(col="red"))
      grid.text( label = "Active-Low", x = 0.11, y = 0.58, rot = 90, default.units = "npc", gp=gpar(col="red"))
      grid.text( label = "Active-High", x = 0.11, y = 0.25, rot = 90, default.units = "npc", gp=gpar(col="red"))
      grid.text( label = "Widespread", x = 0.66, y = 0.025, rot = 0, default.units = "npc", gp=gpar(col="red"))
      
      dev.off() 
  
  # codes with elevation
  plot(elev,cbRaster, ylim=c(0,30))
  
  
  # spatial pattern correlation
  # library(pcaPP)
  # somTime$kendall<-NA
  # temp<-subLayers
  #   temp[is.na(temp)] <- 0
  #
  somTime$spearman<-NA
  somTime$pearson<-NA
  somTime$rmse<-NA
  somTime$mae<-NA
  for(i in 1:nrow(somTime)){
    somTime$pearson[i]<-cor(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]),
                              use = "na.or.complete", method="pearson")
    somTime$spearman[i]<-cor(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]),
        use = "na.or.complete", method="spearman")
    # somTime$kendall[i]<-cor(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]),
    #                          use = "na.or.complete", method="kendall")
    #somTime$kendall[i]<-cor.fk(values(cbRaster[[somTime$mapUnit[i]]]), values(temp[[i]]))
    
    somTime$rmse[i]<- sqrt(sum((values(subLayers[[i]])-values(cbRaster[[somTime$mapUnit[i]]]))^2, na.rm = TRUE)/length(values(subLayers[[i]])))
    somTime$mae[i]<- sum(abs(values(subLayers[[i]])-values(cbRaster[[somTime$mapUnit[i]]])), na.rm = TRUE)/length(values(subLayers[[i]]))
    print(i)
  }
  #rm(temp)
  
  mean(somTime$spearman, na.rm=TRUE)
  median(somTime$spearman, na.rm=TRUE)
  
  mean(somTime$rmse, na.rm=TRUE)
  median(somTime$rmse, na.rm=TRUE)
  
  mean(somTime$pearson, na.rm=TRUE)
  median(somTime$pearson, na.rm=TRUE)
  
  #mean(somTime$kendall, na.rm=TRUE)
  #median(somTime$kendall, na.rm=TRUE)
  #mean(somTime$pearson, na.rm=TRUE)
  
  ##### 
  # Kendall tau corrs of each day against each map unit
  # library(pcaPP)
  # temp<-subLayers
  # temp[is.na(temp)] <- 0
  # 
  # tauCodeBook <- data.frame(matrix(vector(), 0, length(codeList),
  #                        dimnames=list(c(), codeList)),
  #                 stringsAsFactors=F)
  # 
  # for(i in 1:nrow(somTime)){
  #     for(j in 1:length(codeList)){
  #       #tauCodeBook[i,j]<-cor.fk(values(cbRaster[[j]]), values(temp[[i]]))
  #       tauCodeBook[i,j]<-cor(values(cbRaster[[j]]), values(temp[[i]]),
  #                                use = "na.or.complete", method="spearman")
  #       
  #     }
  #   print(i)
  # }
  # somTime$maxTau<-max.col(tauCodeBook)
  # somTime$maxTauVal<-apply(tauCodeBook, 1, FUN=max)
  # somTime$maxTauCode<-codeList[somTime$maxTau]
  # somTime$codeDiff<-somTime$codes==somTime$maxTauCode
  #    # if tau - NA, replace with 1_1
  #       somTime$maxTau[is.na(somTime$maxTau)]<-1
  # somTime$unitDiff<-somTime$mapUnit-somTime$maxTau      
  
  #####
  

  # plot highest correlating days
  highestCorr<-somTime %>% group_by(codes) %>% slice_min(rmse, n=1) # lowest class error
  highestCorr<-somTime %>% group_by(codes) %>% slice_min(errorDist, n=1) # lowest class error
  highestCorr<-somTime %>% group_by(codes) %>% slice_max(spearman, n=1)
  highestCorr<-somTime %>% group_by(codes) %>% slice_max(pearson, n=1)
  highestCorr<-somTime %>% group_by(codes) %>% slice_min(spearman, n=1)
  highestCorr<-somTime %>% group_by(codes) %>% slice_min(pearson, n=1)
  corrRaster<-stack()
  for (i in 1:length(codeList)) {
  temp<-subLayers[[which(somTime$date==highestCorr$date[which(highestCorr$mapUnit==i)])]]
    corrRaster<-stack(corrRaster, temp)
  }
  corrRaster[corrRaster == 0] <- NA  
  at<-c(seq(0,50,1),100)
  mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  pHighCorr<-levelplot(corrRaster, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,
                     main="Highest Spearman June 15th-Sept 30th 4x4 SOM - PRISM-daily 1981-2020")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))

  # max Tau
  # highestCorr<-somTime %>% group_by(maxTauCode) %>% slice_max(maxTauVal, n=1)
  # corrRaster<-stack()
  # for (i in 1:length(codeList)) {
  #   temp<-subLayers[[which(somTime$date==highestCorr$date[which(highestCorr$maxTau==i)])]]
  #   corrRaster<-stack(corrRaster, temp)
  # }
  # corrRaster[corrRaster == 0] <- NA  
  # at<-c(seq(0,50,1),100)
  # mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  # pHighCorr<-levelplot(corrRaster, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,
  #                      main="Highest tau JAS 3x4 SOM - PRISM-daily 1981-2019")+
  #   layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  
  
  # inter-node correlations
   # library(corrplot)
   # corrmtx<-cor(t(codebook[,3:ncol(codebook)]), method = "spearman")
   # corrplot(corrmtx, type = "upper", 
   #          tl.col = "black", tl.srt = 45)
  
  # sum of square errors
  #i=4
  #test<-sum((values(subLayers[[i]])-values(cbRaster[[somTime$mapUnit[i]]]))^2, na.rm = TRUE) 
  
  # i<-4
  # plot(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]))
  # plot(stack(cbRaster[[somTime$mapUnit[i]]], subLayers[[i]]))
  #   plot(values(cbRaster[[somTime$mapUnit[i]]]), type="l")
  #   lines(values(subLayers[[i]]), col="red")
  
  # PRECIP METRICS FOR DAYS
  # dist of all values
  # quantile(as.matrix(df.wide[2:ncol(df.wide)]), probs = c(.25, .5, .75))
    
  # metrics of spatial autocorr?
  # spatial patterns of precip, metrics for each day
  somTime$percExtent<-(rowSums(df.wide[idx,2:ncol(df.wide)]>0)/(ncol(df.wide)-1))*100 # percent extent >0
  somTime$percExtent1<-(rowSums(df.wide[idx,2:ncol(df.wide)]>1)/(ncol(df.wide)-1))*100 # percent extent >1
  somTime$percExtent5<-(rowSums(df.wide[idx,2:ncol(df.wide)]>=5)/(ncol(df.wide)-1))*100 # percent extent >5 mm
  somTime$percExtent10<-(rowSums(df.wide[idx,2:ncol(df.wide)]>=10)/(ncol(df.wide)-1))*100 # percent extent >0
   #somTime$percExtent25<-(rowSums(df.wide[idx,2:ncol(df.wide)]>=25)/(ncol(df.wide)-1))*100 # percent extent >0
  somTime$maxPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, max) # max value of day
  somTime$meanPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, mean)# mean regional precip
  somTime$medPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, median) # max value of day  
  somTime$percZero<-(rowSums(df.wide[idx,2:ncol(df.wide)]==0)/(ncol(df.wide)-1))*100 # percent 0 precip  
  somTime$sumPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, sum) # max value of day 
  
  # location of precip max for each day
  xyCentroid<-list()
  for (i in 1:nlayers(subLayers)) {
    xyCentroid[[i]]<- xyFromCell(subLayers[[i]],which.max(values(subLayers[[i]])))
    print(i)
  }
  xyCentroid<-as.data.frame(do.call("rbind", xyCentroid))
  maxDaily<-cbind.data.frame(somTime[,c("codes","maxPrecip")],xyCentroid)
  colnames(maxDaily)<-c("codes","maxPrecip","lon","lat")
  ggplot(maxDaily, aes(lon,lat, color=maxPrecip))+
    geom_point()+
    geom_density_2d()+
    scale_color_gradient(low="yellow", high="blue",limits=c(0, 100), oob=squish)+
    facet_wrap(~codes)+
    ggtitle("Location of max daily precip by node")
  
  
  
  #####
  # add in climate indices
  load("~/RProjects/SOMs/monsoonPrecip/climInd_BSISO.RData")
  climInd<-climInd[,c("date","phase","amplitude","ONI")]
  climInd$date<-climInd$date+1 # align with PRISM doy
    colnames(climInd)[1]<-"date.c"
    # shift to match ONI
    climInd$ONI<-c(climInd$ONI[-1],NA)
  somTime<-merge(somTime,climInd, by.x="date",by.y="date.c")
   # ENSO phase
  somTime$ENSO<-"Neutral"
    somTime$ENSO<-ifelse(somTime$ONI <= -0.5, "La Nina", somTime$ENSO)   
    somTime$ENSO<-ifelse(somTime$ONI >= 0.5, "El Nino", somTime$ENSO)   
  #####  
  # make dummy date for plotting boxplots
  somTime$dummyDate<-as.Date(paste0("2000-",somTime$month,"-",somTime$day),format="%Y-%m-%d")
    
#####   
# PICK UP WITH PROCESSED DATA - 4/27/22
# save.image("~/RProjects/SOMs/saved_workspace_062922.RData")
  load("~/RProjects/SOMs/saved_workspace_062922.RData")  
#####    
      
  # transition probabilities
  library(gmodels)
    transProb<-cbind.data.frame(somTime[,c(1,2,4:8)], c(somTime[2:nrow(somTime),8],NA))
      colnames(transProb)[8]<-"day2code"
    transProb$day2code[which(transProb$month==9 & transProb$day==30)]<-NA
    crossStats<-CrossTable(transProb$codes,as.character(transProb$day2code), expected = TRUE)
    # viz table https://murraylax.org/rtutorials/viscrosstabs.html
    crossStats<-as.data.frame(table(transProb$codes,as.character(transProb$day2code)))
    ggplot(crossStats, aes(x=Var1, y=Freq, fill=Var2)) + 
      #geom_col() +
      geom_bar(position="fill", stat="identity") +
      guides(fill=guide_legend(title="Day 2"))+
      xlab("Day 1")+
      ylab("count")+
      ggtitle("Transition counts")
      #ylim(0,100)
  # angle between transition
    temp <- transProb %>% separate(codes, c("Y1","X1"), sep = "_", remove=FALSE)
    temp <- temp %>% separate(day2code, c("Y2","X2"), sep = "_", remove=FALSE)
    temp$Y1<-as.numeric(temp$Y1)*-1; temp$Y2<-as.numeric(temp$Y2)*-1;
    temp$X1<-as.numeric(temp$X1); temp$X2<-as.numeric(temp$X2);
    # arctan of coords
    temp$delta_x = temp$X2 - temp$X1
    temp$delta_y = temp$Y2 - temp$Y1
    temp$angle = (atan2(temp$delta_y, temp$delta_x)*180)/(pi)
    temp$angle = ifelse(temp$angle<0, 360+temp$angle, temp$angle)
    temp<-subset(temp, abs(delta_x)<=1 & abs(delta_y)<=1)
    temp<-subset(temp, !(abs(delta_x)==0 & abs(delta_y)==0))
    temp$angle<-temp$angle-5
    # ..count.. or ..ncount..
    ggplot(temp, aes(x = angle)) +
      geom_histogram(aes(y=..ncount../sum(..ncount..)),
                    binwidth = 1, boundary = 5, fill = 'grey', color = 'black', 
                     closed = "left") +
      scale_x_continuous(breaks=seq(0, 270, by=90), limits = c(-5, 355), labels = NULL) +
      coord_polar(start = -(90 * pi/180), direction=-1)+
      scale_y_continuous(limits=c(0, 0.02), oob=squish, labels=NULL)+
      facet_wrap(~codes)+
      xlab(NULL)+
      ylab(NULL)+
      theme(axis.text.x = element_text(size = 8))+
      ggtitle("Transition trajectory frequencies (relative proportions)")
    # code sequence counts
    transProb<-cbind.data.frame(somTime[,c(1,2,4:8)], c(somTime[2:nrow(somTime),8],NA),c(somTime[3:nrow(somTime),8],NA,NA))
      colnames(transProb)[8:9]<-c("day2code","day3code")
      transProb$multiDay<-paste(transProb$codes,transProb$day2code,transProb$day3code)
      transProb$twoDay<-paste(transProb$codes,transProb$day2code)
    multiDay<-as.data.frame(table(transProb$multiDay))
    twoDay<-as.data.frame(table(transProb$twoDay))
    # activity sequence counts
    transProb<-cbind.data.frame(somTime[,c("date","mapUnit","month","year","day","doy","activityCat")], c(somTime[2:nrow(somTime),"activityCat"],NA),c(somTime[3:nrow(somTime),"activityCat"],NA,NA))
      colnames(transProb)[8:9]<-c("day2code","day3code")
      transProb$multiDay<-paste(as.numeric(transProb$activityCat),transProb$day2code,transProb$day3code)
      transProb$twoDay<-paste(as.numeric(transProb$activityCat),transProb$day2code)
    multiDay<-as.data.frame(table(transProb$multiDay))
      twoDay<-as.data.frame(table(transProb$twoDay))
      ggplot(twoDay, aes(Var1,Freq))+
        geom_bar(stat="identity")+
        ggtitle("Activity transition frequencies")
      
  #####    
  # counts of nodes by day of year
    # add in max non 1_1 day
    maxNode<-somTime %>% group_by(month,day) %>% count(codes)
        maxNode<-spread(maxNode, key=codes, value=n)
        maxNode[, "max"] <- apply(maxNode[, 4:ncol(maxNode)], 1, which.max)
        maxNode$maxCount<-maxNode[,maxNode$max]
        maxNode$max<-colnames(maxNode[,4:(ncol(maxNode)-1)])[maxNode$max]
        #maxNode$max<-ifelse(maxNode$`1_1`<20, maxNode$max, NA)  
  
    countDOY<-somTime %>% group_by(month,day) %>% count(codes)
      countDOY <- countDOY %>% group_by(month,day) %>% summarize(maxCount=max(n),
                                                                 maxNode= codes[which.max(n)])
    countDOY$date<-as.Date(paste0(countDOY$month,"-",countDOY$day,"-2016"), format="%m-%d-%Y")      
    countDOY$wday<-wday(countDOY$date, label = T, week_start = 7)
    countDOY$week<-epiweek(countDOY$date)
    countDOY$maxNodeAlt<-maxNode$max
    
    countDOY %>%
      ggplot(aes(wday,-week, fill = maxCount)) +
      geom_tile(colour = "white")  + 
      geom_text(aes(label = maxNodeAlt), size = 3) +
      theme(aspect.ratio = 1/2,
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold", size = 15),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      scale_fill_gradient(low="lightblue", high="red") +
      facet_wrap(~month, nrow = 4, ncol = 1, scales = "free") +
      labs(title = "Most frequent nodes by day of year")
    
    # plot transitions of most common active days
    countDOY<-countDOY %>%
      separate(maxNodeAlt, c("row", "col"), "_")
    countDOY$row<-as.numeric(countDOY$row)
    countDOY$col<-as.numeric(countDOY$col)
    
    countDOY$row <- jitter(countDOY$row)
    countDOY$col <- jitter(countDOY$col)
    
  subset(countDOY,is.na(countDOY$row)==FALSE)%>%
    ggplot()+
      geom_point(aes(col, row, color=date))+
     #geom_path(aes(col, row, color=date))+
     scale_y_reverse(breaks=c(1,2,3,4), labels=c(1,2,3,4))+
     scale_x_continuous(breaks=c(1,2,3,4), labels=c(1,2,3,4))+
     geom_hline(yintercept = c(1.5,2.5,3.5))+
     geom_vline(xintercept = c(1.5,2.5,3.5))+
     expand_limits(x = 4.5, y = 4.5)+
      scale_color_date(low = "blue",high = "yellow")+
    #scale_color_gradient2(low="blue",mid = "green",high = "red")+
    ggtitle("most common active day progression")
   
  # make curvy
  # https://stackoverflow.com/questions/34473292/adding-slight-curve-or-bend-in-ggplot-geom-path-to-make-path-easier-to-read
    test<-data.frame(xspline(countDOY[22:40,c("row","col")], shape=-1, draw=F))

    
  # max node/day top 3 with counts
    temp<-somTime %>% group_by(month,day) %>% count(codes)  
    temp<-temp %>% group_by(month,day) %>% slice_max(n, n = 3)        
    temp$date<-as.Date(paste0("2000-",temp$month,"-",temp$day),"%Y-%m-%d")
    temp$percN<-(temp$n/length(unique(somTime$year)))*100
    
    ggplot(temp, aes(date,percN, fill=codes))+
      geom_bar(stat="identity")+
      ylab("% of days")+
      ggtitle("Most frequent node/day")+
      theme_bw()
    
    temp<-somTime %>% group_by(month,day,activityCat) %>% count(codes)  
    temp<-temp %>% group_by(month,day) %>%  slice_max(n, n = 1)        
    temp$date<-as.Date(paste0("2000-",temp$month,"-",temp$day),"%Y-%m-%d")
    temp$percN<-(temp$n/length(unique(somTime$year)))*100
    
    ggplot(temp, aes(date,percN, fill=codes))+
      geom_bar(stat="identity")
    
    #####  
    
  # distribution of daily precip values by node
  #extentPerc<-somTime[,c("codes","percExtent","percExtent5","percExtent10")]  
  extentPerc<-somTime[,c("codes","percExtent1","percExtent10")]    
    extentPerc<-melt(extentPerc, id.vars="codes")
    ggplot(extentPerc, aes(x=codes, y=value, fill=variable))+
      geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
      ggtitle("Distribution of Daily Extent Precip (%) by nodes")+
      theme(legend.position="bottom")+
      facet_wrap(~codes, scales="free_x")
    # tests of differences
    oneway.test(percExtent ~ codes, data = somTime, var.equal = TRUE)
    pairwise.t.test(somTime$percExtent10, somTime$codes)
    # distribution by activity
    #extentPerc<-somTime[,c("activityCat","percExtent","percExtent5","percExtent10")] 
    extentPerc<-somTime[,c("activityCat","percExtent1","percExtent10")] 
    extentPerc<-melt(extentPerc, id.vars="activityCat")
    ggplot(extentPerc, aes(x=activityCat, y=value, fill=variable))+
      geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
      ggtitle("Distribution of Daily Extent Precip (%) by nodes")+
      theme(legend.position="bottom")
      #facet_wrap(~codes, scales="free_x")
    # dist by codes/activity
    #extentPerc<-somTime[,c("codes","activityCat","percExtent","percExtent5","percExtent10")]
    extentPerc<-somTime[,c("codes","activityCat","percExtent10")] 
    extentPerc<-melt(extentPerc, id.vars=c("codes","activityCat"))
    ggplot(extentPerc, aes(x=reorder(codes, value), y=value, fill=variable))+
      geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2, outlier.shape = NA)+
      ggtitle("Distribution of Daily Extent Precip (%) by nodes")+
      theme(legend.position="bottom")+
      facet_wrap(~activityCat, scales="free_x",nrow = 1)
    #pairwise.t.test(somTime$percExtent1, somTime$activityCat)
    
    # within activity class
    temp<-subset(somTime, activityCat=="Active-high")
    pairwise.t.test(temp$percExtent10, temp$codes)
    kruskal.test(percExtent1 ~ codes, data = temp)                
    pairwise.wilcox.test(temp$percExtent10, temp$codes,
                         p.adjust.method = "BH")
    
    
    # dist of max precip by codes/activity
    extentPerc<-somTime[,c("codes","activityCat","maxPrecip")]  
    extentPerc<-melt(extentPerc, id.vars=c("codes","activityCat"))
    ggplot(extentPerc, aes(x=codes, y=value, fill=variable))+
      geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
      ggtitle("Distribution of Daily Max by nodes")+
      theme(legend.position="bottom")+
      facet_wrap(~activityCat, scales="free_x",nrow = 1)
    pairwise.t.test(somTime$maxPrecip, somTime$activityCat)
    pairwise.t.test(somTime$maxPrecip, somTime$codes)
    
    # dist of median precip by codes/activity
    extentPerc<-somTime[,c("codes","activityCat","medPrecip")]  
    extentPerc<-melt(extentPerc, id.vars=c("codes","activityCat"))
    ggplot(extentPerc, aes(x=codes, y=value, fill=variable))+
      geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
      ggtitle("Distribution of Daily Median by nodes")+
      theme(legend.position="bottom")+
      facet_wrap(~activityCat, scales="free_x",nrow = 1)
    pairwise.t.test(somTime$medPrecip, somTime$activityCat)
    pairwise.t.test(somTime$medPrecip, somTime$codes)
    # dist of median precip by codes/activity
    extentPerc<-somTime[,c("codes","activityCat","meanPrecip")]  
    extentPerc<-melt(extentPerc, id.vars=c("codes","activityCat"))
    ggplot(extentPerc, aes(x=codes, y=value, fill=variable))+
      geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
      ggtitle("Distribution of Daily Mean by nodes")+
      theme(legend.position="bottom")+
      facet_wrap(~activityCat, scales="free_x",nrow = 1)
    pairwise.t.test(somTime$meanPrecip, somTime$activityCat)
    pairwise.t.test(somTime$meanPrecip, somTime$codes)
    
    
    # k means clustering on extent
    #kmeanExtent<-kmeans(somTime[,c("percExtent","percExtent5","percExtent10")],4)
    kmeanExtent<-kmeans(somTime[,c("percExtent1","percExtent10")],4)
    tableAct<-(table(cbind.data.frame(somTime$activityCat,kmeanExtent$cluster)))
    tableCodes<-(table(cbind.data.frame(somTime$codes,kmeanExtent$cluster)))
      # kmean stats
        temp<-melt(kmeanExtent$centers)
        ggplot(temp, aes(Var1,value, fill=Var2))+
          geom_col()
        somTime$kmeans<-kmeanExtent$cluster
        ggplot(somTime, aes(kmeans))+
          geom_histogram()+
          facet_wrap(~codes)+
          ylim(0,500)
            
    
  ggplot(somTime, aes(as.factor(somTime$codes), percExtent1))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Daily Extent Precip (%) by nodes")
  ggplot(somTime, aes(as.factor(somTime$codes), percExtent5))+
      geom_boxplot(varwidth = TRUE)+
      ggtitle("Distribution of Daily Extent Precip>5mm (%) by nodes")    
  ggplot(somTime, aes(as.factor(somTime$codes), percExtent10))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Daily Extent Precip>10mm (%) by nodes") 
      
  ggplot(somTime, aes(as.factor(somTime$codes), maxPrecip))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Max Daily Precip (mm) by nodes")
  ggplot(somTime, aes(as.factor(somTime$codes), meanPrecip))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Mean Daily Precip (mm) by nodes")  
  ggplot(somTime, aes(as.factor(somTime$codes), medPrecip))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Median Daily Precip (mm) by nodes")  
  ggplot(somTime, aes(as.factor(somTime$codes), percZero))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Percent Zero Extent Precip (mm) by nodes")  
  ggplot(somTime, aes(as.factor(somTime$codes), rmse))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of RMSE (mm) by nodes") 
  ggplot(somTime, aes(as.factor(somTime$codes), mae))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of MAE (mm) by nodes")  
  ggplot(somTime, aes(as.factor(somTime$codes), errorDist))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of classification error (mm^2) by nodes") 

# summary stats by codes
  somTime %>% 
    group_by(codes) %>% 
    summarize(median = median(percExtent),
              mean = mean(percExtent),
              min = min(percExtent),
              max = max(percExtent))
  
  # summary stats by activity classes
  somTime %>% 
    group_by(activityCat) %>% 
    summarize(median = median(percExtent),
              mean = mean(percExtent),
              min = min(percExtent),
              max = max(percExtent))
  # percent of days in activit classes
  table(somTime$activityCat)/nrow(somTime)
  
    
  # plot 10 random maps from selected node to assess quality
  at<-c(seq(0.01,50,1),200)
  mapTheme <- rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
  temp<-subLayers[[idx]]
  temp[temp == 0] <- NA  
  pRand<-levelplot(temp[[sample(which(somTime$codes=="3_3"),10)]], contour=FALSE, at=at,
                     margin=FALSE, par.settings=mapTheme, 
                     main="Precip Patterns from 3_3 - PRISM-daily 1981-2020")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  
  # plot single day
  at<-c(seq(0,30,1),100)
  levelplot(temp[[846]], contour=FALSE,
            margin=FALSE, par.settings=mapTheme, at=at,
            main="Precip")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  
  # maps on composites of days in nodes ####
  # percentile data for composites
  # perc<-stack("/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2019_PRISM_daily_prcp.grd") 
  #    perc<-crop(perc,e)
      # apply mask
  #    temp<-mask(perc, mask)
  # composites on raw precip
  temp<-subLayers[[idx]]
  comp<-stackApply(temp, somTime$mapUnit, fun=max)
  #comp<-stackApply(temp, somTime$maxTau, fun=median)
  #comp2<-stackApply(temp, somTime$mapUnit, fun=mean)
  #comp<-(comp/comp2)*100
    tempNames<-as.data.frame(names(comp))
    colnames(tempNames)<-"index"
    tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
    comp<-subset(comp, order(as.numeric(tempNames$mapUnit)))
  names(comp)<-codeList   
  comp[comp ==0] <- NA 
  # plot 
  at<-c(seq(0,100,5),200)
  #at<-c(seq(0,200,5),2000)
  at<-c(seq(0,30,2.5))
  mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,
                     main="Composite Max Precip June 15th-Sep 30th 4x4 SOM - PRISM-daily 1981-2020", xlab=NULL, ylab=NULL)+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))+
    contourplot(elev, at=c(1000,2000,3000), labels=FALSE, lwd = 0.3, par.settings = GrTheme)
  # plot median 
  at<-c(seq(0,1,0.05))
  mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
  pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at, 
                   main="Composite Median Percentile Precip JAS 3x4 SOM - PRISM-daily 1981-2020", xlab=NULL, ylab=NULL)+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  #####
      ## location of max values on max composite
      maxIdx<-which.max(comp)
      maxIdx<-as.factor(maxIdx)
      rat <- levels(maxIdx)[[1]]
      rat[["codes"]]<- codeList[rat$ID]
      levels(maxIdx)<- rat
      mapTheme <- rasterTheme(region=brewer.pal(7,"Set1"))
      levelplot(maxIdx, par.settings=mapTheme, main="Node contributing max codebook value")+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  ####    
  
#### Prob of exceedance based on WMO thresholds  PoE    
  # composites on raw precip
    temp<-subLayers[[idx]]
    # set threshold
    comp<-stackApply(temp, somTime$mapUnit, fun=function(x, ...){(sum(x >= 20)/length(x))*100})
    # set names and reorder
    tempNames<-as.data.frame(names(comp))
    colnames(tempNames)<-"index"
    tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
    comp<-subset(comp, order(as.numeric(tempNames$mapUnit)))
    names(comp)<-codeList   
    comp[comp ==0] <- NA 
    # plot
    at<-c(seq(0,100,10))
    mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","orange3","red", "red4"))
    pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,
                     main="PoE > 20 Precip June 15th-Sep 30th 4x4 SOM - PRISM-daily 1981-2020", xlab=NULL, ylab=NULL)+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))+
      contourplot(elev, at=c(1000,2000,3000), labels=FALSE, lwd = 0.3, par.settings = GrTheme)
    #   # codes with elevation
    #plot(elev,comp, ylim=c(0,100))  
    #boxplot(comp)
    #bwplot(comp)
    
    #####
    # EXTENDED PoE Plot
    # activity text for maps    
    # temp<-(unique(somTime[,c("mapUnit","activityCat")]))
    # temp<-temp[order(temp[,1]),]
    # text2add<-as.character(temp[,2])
    # # color of text
    # textCol<-c("green","yellow","orange","red",
    #            "yellow","yellow","orange","red",
    #            "orange","orange","orange","red",
    #            "red","red","red","red")
    
    at<-c(seq(0,100,10))
    #names(comp)<-codeList
    col.titles<-paste0("",codeList," (",round((table(somTime$codes)/nrow(somTime))*100,1),"%, ",
                       round(table(somTime$codes)/length(unique(somTime$year)),1),"dy/yr)")
    at<-c(seq(0,100,10))
    mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","orange3","red", "red4"))
    pComp<-levelplot(comp, contour=FALSE, margin=FALSE, layout=c(ncols,nrows), at=at,
                     names.attr=col.titles,
                     par.settings=mapTheme,
                     scales=list(draw=FALSE),
                     main="PoE > 10mm Precip June 15th-Sep 30th 4x4 SOM - PRISM-daily 1981-2020", xlab=NULL, ylab=NULL)+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))+
    layer(sp.points(stations, col="black", pch = 18, cex=0.5))+
    contourplot(elev, at=c(1000,2000,3000), labels=FALSE, lwd = 0.3, par.settings = GrTheme)+
      #layer(sp.text(c(36.0840,-115.1537),"LAS"))+
      #layer(sp.points(xyCentroid[panel.number()],
      #              pch=20, cex=1, col="black"))+
      #layer(panel.rect(-115.5,31.3,-113.0,32.0,fill=textCol[panel.number()]))+
      #layer(panel.rect(-115.5,31.3,-113.0,32.0, border=textCol[panel.number()]))+
    layer(panel.text(-114.25, 31.55, text2add[panel.number()],col="black",cex=0.4))
    #layer(panel.text(stationCodes[panel.number(1),3], stationCodes[panel.number(1),2],
    #                 stationCodes[panel.number(1),1],col="black",cex=0.4))
    
    #layer(sp.polygons(huc4, col = 'gray20', lwd=0.75))
    #layer(sp.polygons(us, col = 'gray40', lwd=1))+
    #layer(sp.polygons(mx, col = 'gray40', lwd=1))
    
    # add lines
    #grid.ls(viewport=TRUE, grobs=FALSE) #plot_01.toplevel.vp::plot_01.panel.1.1.vp
    #grid.rect(vp = "plot_01.toplevel.vp::plot_01.panel.1.1.vp",
    #          gp = gpar(col = "red"))
    png("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/precipPoE10mm_SOM_4x4.png",
        width = 10, height = 6.8, units = "in", res = 300L)
    pComp
    # vert line 
    grid.lines(x = unit(c(0.2465, 0.2465), "npc"),
               y = unit(c(0.75, 1), "npc"),
               default.units = "npc",
               arrow = NULL, name = NULL,
               gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
    grid.lines(x = unit(c(0.5, 0.5), "npc"),
               y = unit(c(0.497, 1), "npc"),
               default.units = "npc",
               arrow = NULL, name = NULL,
               gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")    
    grid.lines(x = unit(c(0.753, 0.753), "npc"),
               y = unit(c(0.246, 1), "npc"),
               default.units = "npc",
               arrow = NULL, name = NULL,
               gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp") 
    # horiz line 
    grid.lines(x = unit(c(0, 0.2465), "npc"),
               y = unit(c(0.75,0.75), "npc"),
               default.units = "npc",
               arrow = NULL, name = NULL,
               gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")  
    grid.lines(x = unit(c(0, 0.5), "npc"),
               y = unit(c(0.497,0.497), "npc"),
               default.units = "npc",
               arrow = NULL, name = NULL,
               gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
    grid.lines(x = unit(c(0, 0.753), "npc"),
               y = unit(c(0.246,0.246), "npc"),
               default.units = "npc",
               arrow = NULL, name = NULL,
               gp=gpar(lwd=2, col="red"), draw = TRUE, vp = "plot_01.toplevel.vp::plot_01.figure.vp")
    grid.text( label = "Inactive", x = 0.11, y = 0.80, rot = 90, default.units = "npc", gp=gpar(col="red"))
    grid.text( label = "Active-Low", x = 0.11, y = 0.58, rot = 90, default.units = "npc", gp=gpar(col="red"))
    grid.text( label = "Active-High", x = 0.11, y = 0.36, rot = 90, default.units = "npc", gp=gpar(col="red"))
    grid.text( label = "Widespread", x = 0.11, y = 0.14, rot = 90, default.units = "npc", gp=gpar(col="red"))
    #print(pComp, newpage = FALSE)  
    dev.off() 
    #####
      
    #####
    # composites based on activity classes
    temp<-subLayers[[idx]]
    # set threshold
    comp<-stackApply(temp, somTime$activityCat2, fun=function(x, ...){(sum(x >= 1)/length(x))*100})
    # set names and reorder
    #comp<-subset(comp, order(c(4,1,2,3)))
    #comp<-subset(comp, order(c(2,1)))
    comp[comp ==0] <- NA 
    # plot
    at<-c(seq(0,100,10))
    mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","orange3","red", "red4"))
    pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(2,1), at=at,
                     main="PoE > 1mm Precip Jun 15th-Sep 30th 4x4 SOM - PRISM-daily 1981-2020", xlab=NULL, ylab=NULL)+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))+
      contourplot(elev, at=c(1000,2000,3000), labels=FALSE, lwd = 0.3, par.settings = GrTheme)
    
    #
    
    #####
    # tabulation of activity classes
    activityTab<-as.data.frame(table(somTime$activityCat2, somTime$year))
      colnames(activityTab)<-c("Class","Year","Count")
      activityTab<-spread(activityTab,Class,Count)
      activityTab$InactiveAnom<-scale(activityTab$Inactive, scale=FALSE)
      activityTab$ActiveAnom<-scale(activityTab$Active, scale=FALSE)
    #table(somTime$activityCat)
    #table(somTime$codes)
#####
# SEASONAL TIMING and LENGTH     
      # cumulative timing of activity classes
      temp<-somTime[,1:18]
      temp$ActiveYes<-ifelse(temp$activityCat2=="Active",1,0)
      # active based on median extent
      temp$ActiveYesMed<-ifelse(temp$percExtent1>=median(temp$percExtent1),1,0)
      ggplot(temp, aes(activityCat,percExtent1))+
        geom_boxplot()
      plot(subLayers[[which(somTime$date=="2010-09-15")]])
      
      # cumulative ActiveYes
      cumActive <- temp %>% 
        dplyr::group_by(year, doy) %>% 
        dplyr::summarise(value = sum(ActiveYes)) %>%
        dplyr::mutate(csum = cumsum(value)) 
      temp$ActiveYesCsum<-cumActive$csum
      # cumulative ActiveYes
      cumActive <- temp %>% 
        dplyr::group_by(year, doy) %>% 
        dplyr::summarise(value = sum(ActiveYesMed)) %>%
        dplyr::mutate(csum = cumsum(value))
      temp$ActiveYesMedCsum<-cumActive$csum
      
      avgActiveYes <- temp %>%
        dplyr::group_by(doy) %>% # 
        dplyr::summarise(value = mean(ActiveYesCsum))
      
      ggplot(temp, aes(doy,ActiveYesCsum,color=as.factor(year)))+
        geom_line()+
      geom_point(data=avgActiveYes,aes(doy,value, color="Avg"))
      
      # onset/demise dates based on runs of 3 days 1_1 active
      # try with rle2
      temp<-somTime[,1:18]
        temp$ActiveYes<-ifelse(temp$activityCat2=="Active",1,0)
        # active based on median extent
        temp$ActiveYesMed<-ifelse(temp$percExtent1>=median(temp$percExtent1),1,0)
        # do by year
        yrs<-seq(1981,2020,1)
        df_total<-data.frame()
        for(i in 1:length(yrs)){
          temp2<-subset(temp, year==yrs[i])
          rle.seq1 <- rle(temp2$ActiveYes)
            index <- which(rle.seq1$values==1 & rle.seq1$lengths >=3) # based on 3-day run
            temp2$day3run<-NA
            temp2$day3run[cumsum(rle.seq1$lengths)[index]]<-1
          newindex <- ifelse(index > 1, index-1,0)
            starts = cumsum(rle.seq1$lengths)[newindex]+1
            if(0 %in% newindex)starts=c(1,starts)
            starts        
            temp2$day3run<-NA
            temp2$day3run[starts]<-1
          tempOut <- temp2$day3run
          df_total <- rbind(df_total,data.frame(tempOut))
        }
        temp$day3run<-df_total$tempOut

      # get onset date by year
        tempOnset<- temp %>% group_by(year) %>%
                    summarize(first = date[min(which(!is.na(day3run)))])
        tempOnset$fakeDate<-as.Date(paste0("2000-",format(tempOnset$first,"%m"),"-",format(tempOnset$first,"%d")))
      # median median threshold dates
      df_total<-data.frame()
      for(i in 1:length(yrs)){
          temp2<-subset(temp, year==yrs[i])
          rle.seq1 <- rle(temp2$ActiveYesMed)
            index <- which(rle.seq1$values==1 & rle.seq1$lengths >=3) # based on 3-day run
            temp2$day3runMed<-NA
            temp2$day3runMed[cumsum(rle.seq1$lengths)[index]]<-1
          newindex <- ifelse(index > 1, index-1,0)
            starts = cumsum(rle.seq1$lengths)[newindex]+1
            if(0 %in% newindex)starts=c(1,starts)
            starts        
            temp2$day3runMed<-NA
            temp2$day3runMed[starts]<-1
          tempOut <- temp2$day3runMed
          df_total <- rbind(df_total,data.frame(tempOut))
      }
      temp$day3runMed<-df_total$tempOut    
            
        # get onset date by year
        tempOnset2<- temp %>% group_by(year) %>%
          summarize(first = date[min(which(!is.na(day3runMed)))])
        tempOnset2$fakeDate<-as.Date(paste0("2000-",format(tempOnset2$first,"%m"),"-",format(tempOnset2$first,"%d")))
        
        # combine data frames
        tempOnset$type<-"SOM"
        tempOnset2$type<-"MedianThreshold"
        tempOnset<-rbind.data.frame(tempOnset[,c("year","fakeDate","type")],tempOnset2[,c("year","fakeDate","type")])
        
        ggplot(tempOnset, aes(fakeDate, fill=type))+
          geom_boxplot()+
          ggtitle("Distribution of Onset using Active 3-day")
        
        onsets<-spread(tempOnset,type,fakeDate)
        ggplot(onsets, aes(MedianThreshold,SOM, label=year))+
          geom_point()+
          geom_text()+
          geom_abline(intercept = 0, slope = 1)+
          ggtitle("3-day Onset Dates")
         
        mean(onsets$SOM)
        mean(onsets$MedianThreshold)
        median(onsets$SOM)
        median(onsets$MedianThreshold)
        sd(onsets$SOM)
        sd(onsets$MedianThreshold) 
   
      # ending dates
        runLength<-3
        temp<-somTime[,1:18]
          # recode to 0/1 - inactive
          #temp$ActiveYes<-ifelse(temp$activityCat2=="Active",0,1)
          #temp$ActiveYesMed<-ifelse(temp$percExtent1>=median(temp$percExtent1),0,1)
          # recode to 0/1 -active
          temp$ActiveYes<-ifelse(temp$activityCat2=="Active",1,0)
          temp$ActiveYesMed<-ifelse(temp$percExtent1>=median(temp$percExtent1),1,0)
        # do by year
        yrs<-seq(1981,2020,1)
        df_total<-data.frame()
        for(i in 1:length(yrs)){
          temp2<-subset(temp, year==yrs[i])
          rle.seq1 <- rle(temp2$ActiveYes)
          index <- which(rle.seq1$values==1 & rle.seq1$lengths >=runLength) # based on 3-day run
          temp2$day3run<-NA
          temp2$day3run[cumsum(rle.seq1$lengths)[index]]<-1
            # get beginning of seq
            # newindex <- ifelse(index > 1, index-1,0)
            # starts = cumsum(rle.seq1$lengths)[newindex]+1
            # if(0 %in% newindex)starts=c(1,starts)
            # #starts        
            # temp2$day3run<-NA
            # temp2$day3run[starts]<-1
          tempOut <- temp2$day3run
          df_total <- rbind(df_total,data.frame(tempOut))
        }
        temp$day3run<-df_total$tempOut
        
        # get onset date by year
        tempEnd<- temp %>% group_by(year) %>%
          summarize(first = date[max(which(!is.na(day3run)))])
        tempEnd$fakeDate<-as.Date(paste0("2000-",format(tempEnd$first,"%m"),"-",format(tempEnd$first,"%d")))
        
        # median median threshold dates
        df_total<-data.frame()
        for(i in 1:length(yrs)){
          temp2<-subset(temp, year==yrs[i])
          rle.seq1 <- rle(temp2$ActiveYesMed)
          index <- which(rle.seq1$values==1 & rle.seq1$lengths >=runLength) # based on 3-day run
          temp2$day3runMed<-NA
          temp2$day3runMed[cumsum(rle.seq1$lengths)[index]]<-1
            # get at beginning of seq
            # newindex <- ifelse(index > 1, index-1,0)
            # starts = cumsum(rle.seq1$lengths)[newindex]+1
            # if(0 %in% newindex)starts=c(1,starts)
            # starts        
            # temp2$day3runMed<-NA
            # temp2$day3runMed[starts]<-1
          tempOut <- temp2$day3runMed
          df_total <- rbind(df_total,data.frame(tempOut))
        }
        temp$day3runMed<-df_total$tempOut    
        
        
        # get onset date by year
        tempEnd2<- temp %>% group_by(year) %>%
          summarize(first = date[max(which(!is.na(day3runMed)))])
        tempEnd2$fakeDate<-as.Date(paste0("2000-",format(tempEnd2$first,"%m"),"-",format(tempEnd2$first,"%d")))
        
        # combine data frames
        tempEnd$type<-"SOM"
        tempEnd2$type<-"MedianThreshold"
        tempEnd<-rbind.data.frame(tempEnd[,c("year","fakeDate","type")],tempEnd2[,c("year","fakeDate","type")])
        
        ggplot(tempEnd, aes(fakeDate, fill=type))+
          geom_boxplot()+
          ggtitle("Distribution of End using Active 3-day")
        
       ends<-spread(tempEnd,type,fakeDate)
        ggplot(ends, aes(MedianThreshold,SOM, label=year))+
          geom_point()+
          geom_text()+
          geom_abline(intercept = 0, slope = 1)+
          ggtitle("3-day Ends Dates")
        
        mean(ends$SOM)
        mean(ends$MedianThreshold)
        median(ends$SOM)
        median(ends$MedianThreshold)
        sd(ends$SOM)
        sd(ends$MedianThreshold) 
        
        # bursts/breaks
        runLength<-3
        temp<-somTime[,1:18]
        # recode to 0/1 - inactive
        #temp$ActiveYes<-ifelse(temp$activityCat2=="Active",0,1)
        #temp$ActiveYesMed<-ifelse(temp$percExtent1>=median(temp$percExtent1),0,1)
        # recode to 0/1 -active
        temp$ActiveBursts<-ifelse(temp$activityCat2=="Active",1,0)
          #temp$ActiveBurstsMed<-ifelse(temp$percExtent1>=median(temp$percExtent1),1,0)
        temp$ActiveBreaks<-ifelse(temp$activityCat2=="Active",0,1)
          #temp$ActiveBreaksMed<-ifelse(temp$percExtent1>=median(temp$percExtent1),0,1)
        # do by year
        yrs<-seq(1981,2020,1)
        df_total<-data.frame()
        for(i in 1:length(yrs)){
          temp2<-subset(temp, year==yrs[i])
          # subset to active period onset to end
          beg<-as.Date(paste0(yrs[i],"-",format(onsets$SOM[which(onsets$year==yrs[i])],"%m-%d")))
          end<-as.Date(paste0(yrs[i],"-",format(ends$SOM[which(ends$year==yrs[i])],"%m-%d")))
          temp2<- temp2[temp2$date >= beg & temp2$date <= end, ]
          # bursts
          rle.seq1 <- rle(temp2$ActiveBursts)
            rle.seq1<-cbind.data.frame(rle.seq1$lengths,rle.seq1$values)
            rle.seq1<-subset(rle.seq1,rle.seq1$`rle.seq1$values`==1)
          # bursts
          rle.seq2 <- rle(temp2$ActiveBreaks)
            rle.seq2<-cbind.data.frame(rle.seq2$lengths,rle.seq2$values)
            rle.seq2<-subset(rle.seq2,rle.seq2$`rle.seq2$values`==1)
            
          df_total <- rbind(df_total,data.frame(yrs[i],max(rle.seq1$`rle.seq1$lengths`),
                                                mean(rle.seq1$`rle.seq1$lengths`),
                                                max(rle.seq2$`rle.seq2$lengths`),
                                                mean(rle.seq2$`rle.seq2$lengths`)))
        }
        colnames(df_total)<-c("year","maxBurst","meanBurst","maxBreak","meanBreak")
        ##### end bursts/breaks
        
        # seasonal precip within seas (onset/end)
        #subDates$sumDate<-as.Date(paste0(subDates$year,"-",subDates$month,"-01"),format="%Y-%m-%d")
        subDates2<-subDates
          #subDates2$seasFlag<-NA
          df_total2<-data.frame()
          yrs<-seq(1981,2020,1)
          for(i in 1:length(yrs)){
            beg<-as.Date(paste0(yrs[i],"-",format(onsets$SOM[which(onsets$year==yrs[i])],"%m-%d")))
            end<-as.Date(paste0(yrs[i],"-",format(ends$SOM[which(ends$year==yrs[i])],"%m-%d")))
            temp2<- cbind.data.frame(seq.Date(beg,end, by = "days"),yrs[i])
            df_total2 <- rbind(df_total2,data.frame(temp2))
          }
        colnames(df_total2)<-c("date","seasFlag")
        subDates2<-merge(subDates2,df_total2, by="date", all=TRUE)  
        # seasonal totals within season  
        sumSeasAdj<-stackApply(subLayers, subDates2$seasFlag, fun = sum)
          seasAvgPrecipAdj<-cellStats(sumSeasAdj[[2:nlayers(sumSeasAdj)]], 'mean')
          seasAvgPrecipAdj<-cbind.data.frame(unique(subDates$year),seasAvgPrecipAdj)
          seasAvgPrecipAdj$percRank<-perc.rank(seasAvgPrecipAdj$seasAvgPrecip) 
          colnames(seasAvgPrecipAdj)<-c("year","avgPrecip","percRank")
          # names
          seasAvgPrecipAdj$anomName<-"normal"
          seasAvgPrecipAdj$anomName[seasAvgPrecipAdj$percRank<=0.333334] <- "dry"
          seasAvgPrecipAdj$anomName[seasAvgPrecipAdj$percRank>0.67] <- "wet"
        
        # season length and stats
        temp2<-as.data.frame(table(temp$year,temp$activityCat2))
          temp2<-subset(temp2, Var2=="Active")
        seasStats<-cbind.data.frame(onsets[,c("year","SOM")],ends[,c("SOM")],seasAvgPrecip[,c("avgPrecip","percRank","anomName")],
                                    temp2[,c("Freq")],df_total[,c("maxBurst","meanBurst","maxBreak","meanBreak")],seasAvgPrecipAdj[,c("avgPrecip")])
          colnames(seasStats)[c(1:3,7,12)]<-c("year","onset","end","ActiveDays","seasAvgPrecipAdj")
        seasStats$length<-seasStats$end-seasStats$onset 
        seasStats$ActIntens<-seasStats$ActiveDays/as.numeric(seasStats$length) 
        seasStats$length<-as.numeric(seasStats$length)
        seasStats$diffPrec<-seasStats$avgPrecip-seasStats$seasAvgPrecipAdj
        
        # add in yearly activity cats
        temp<-as.data.frame(table(somTime$year,somTime$activityCat))
        temp<-spread(temp,Var2,Freq)
        seasStats<-cbind.data.frame(seasStats,temp[,c(2:5)])
          seasStats$allActive<-rowSums(seasStats[,c(17:19)])
        
        ggplot(seasStats, aes(avgPrecip,length,label=year))+
          geom_point()+
          geom_text()+
          ggtitle("Season length vs total avg precip")
        
        ggplot(seasStats, aes(ActiveDays,length,label=year))+
          geom_point()+
          geom_text()+
          ggtitle("Season length vs active days")
        
        ggplot(seasStats, aes(year,ActIntens, label=year))+
          geom_point()+
          geom_line()+
          geom_text()+
          ggtitle("Activity Intensity")
      
        ggplot(seasStats, aes(avgPrecip,seasAvgPrecipAdj,label=year))+
          geom_point()+
          geom_text()+
          geom_abline(intercept = 0, slope = 1)+
          ggtitle("Seas Total Precip vs Within Seas Total")
        
        ggplot(seasStats, aes(seasAvgPrecipAdj,length,label=year))+
          geom_point()+
          geom_text()+
          #geom_abline(intercept = 0, slope = 1)+
          ggtitle("Within Seas Total Precip vs length")
        
        ggplot(seasStats, aes(avgPrecip-seasAvgPrecipAdj,length,label=year))+
          geom_point()+
          geom_text()+
          ggtitle("Season length vs diff in seas-adj total precip")
        
        vtable::st(seasStats)
        tapply(somTime$percExtent,somTime$activityCat2, summary)
        
        ggplot(seasStats, aes(year,avgPrecip-seasAvgPrecipAdj,label=year))+
          geom_point()+
          geom_text()+
          ggtitle("diff in seas-adj total precip")
       
        library("PerformanceAnalytics")
        #chart.Correlation(somTime[,c("percExtent","percExtent5","percExtent10","maxPrecip","meanPrecip","medPrecip")], histogram=TRUE, pch=19)
        chart.Correlation(seasStats[,c("avgPrecip","ActiveDays","maxBurst","meanBurst","maxBreak","meanBreak",
                                       "seasAvgPrecipAdj","length","ActIntens","Inactive","allActive","Active-low","Active-high","Widespread")], histogram=TRUE, pch=19)
        # scatters of correlated vars
        ggplot(seasStats, aes(avgPrecip,Widespread,label=year))+
          geom_point()+
          geom_text()+
          geom_smooth(method='lm')+
          ggtitle("Seas Precip vs Widespread days")
        ggplot(seasStats, aes(avgPrecip,ActiveDays,label=year))+
          geom_point()+
          geom_text()+
          geom_smooth(method='lm')+
          ggtitle("Seas Precip vs Active days")
        ggplot(seasStats, aes(avgPrecip,`Active-high`,label=year))+
          geom_point()+
          geom_text()+
          geom_smooth(method='lm')+
          ggtitle("Seas Precip vs Active-high days")
        ggplot(seasStats, aes(avgPrecip,`Active-low`,label=year))+
          geom_point()+
          geom_text()+
          geom_smooth(method='lm')+
          ggtitle("Seas Precip vs Active-low days")
        
        
        m<-lm(avgPrecip~`Active-high`+Widespread, data=seasStats)
          summary(m)
        m<-lm(avgPrecip~allActive, data=seasStats)
          summary(m)
        m<-lm(avgPrecip~`Active-high`, data=seasStats)
          summary(m)
        m<-lm(avgPrecip~`Active-low`, data=seasStats)
          summary(m)  
          
        
        # add to activity frequency plot
        # activity chart 
        countDOY<-somTime %>% group_by(dummyDate) %>% count(activityCat)
        # moving average smooth
        countDOY<- countDOY %>%
          group_by(activityCat) %>%
          mutate(avg5 = zoo::rollapplyr(
            data = n,
            width = 10,
            FUN = mean,
            by.column = FALSE,
            fill = NA, align="center"))
        countDOY$avg5 <- ifelse(is.na(countDOY$avg5), countDOY$n, countDOY$avg5)
        
        # polynomial smooth
        countDOY<-somTime %>% group_by(dummyDate) %>% count(activityCat)
        #countDOY<-complete(countDOY, activityCat, dummyDate)
        
        # moving average smooth
        countDOY<- countDOY %>%
          group_by(activityCat) %>%
          mutate(smooth=as.numeric(smooth(n, kind = "3RS3R")))
          #mutate(smooth=ksmooth(dummyDate,n, "normal",bandwidth = 1))
          #mutate(smoothCount=KernSmooth::locpoly(x=dummyDate, y=n, bandwidth=KernSmooth::dpill(dummyDate,n,gridsize=n()), gridsize = n())$y)
        #countDOY$avg5 <- ifelse(is.na(countDOY$avg5), countDOY$n, countDOY$avg5)
        
       
        
   
        
        
        ggplot(countDOY, aes(fill=as.factor(activityCat), y=n, x=dummyDate)) + 
          #geom_bar(position="stack", stat="identity")
          geom_bar(position="fill", stat="identity")+
          scale_fill_brewer(type = "qual",
                            palette = "Paired",
                            direction = 1,
                            aesthetics = "fill")+
          ggtitle("Activity frequency by day through season")
        
        
        p1<-ggplot(countDOY, aes(fill=as.factor(activityCat), y=n, x=dummyDate)) + 
          #geom_bar(position="stack", stat="identity")
          geom_bar(position="fill", stat="identity")+
          scale_fill_brewer(type = "qual",
                            palette = "Paired",
                            direction = 1,
                            aesthetics = "fill")+
          geom_hline(yintercept = 0.5)+
          ggtitle("Activity frequency by day through season (10-day centerd mean)")
        
        temp<-gather(seasStats[,2:3],var,dates)
        p2<-ggplot(temp, aes(dates, fill=var))+
          geom_boxplot()+
          ggtitle("Onset and End (3-day active period)")
        plot_grid(p1,p2, nrow=2, align="v")
        
        # multiplot with activity and seasonal precip
        # activity by year
        countDOY<-somTime %>% group_by(year) %>% count(activityCat)
        p1<-ggplot(countDOY, aes(fill=activityCat, y=n, x=year)) + 
          #geom_bar(position="stack", stat="identity")
          geom_bar(position="fill", stat="identity")+
          scale_fill_brewer(type = "qual",
                            palette = "Paired",
                            direction = 1,
                            aesthetics = "fill")+
          ggtitle("Activity frequency by year through season")
        p2<-ggplot(seasStats, aes(year,avgPrecip, fill=as.factor(seasStats$anomName)) )+
          geom_bar(stat = 'identity')+
          ggtitle("Regional Average Total Precip (Jun 15-Sep30, mm)")+
          geom_hline(yintercept=mean(seasAvgPrecip$avgPrecip), color="black")+
          geom_hline(yintercept=median(seasAvgPrecip$avgPrecip), color="red")+
          scale_fill_manual(values = c("saddlebrown", "grey", "forestgreen"), name="tercile")
        # season length
        ggplot(seasStats, aes(year, length))+
          geom_bar(stat="identity")
        # onset/end  
        temp<-gather(seasStats[,1:3],var,dates, -year)
        p3<-ggplot(temp, aes(year,dates, color=var))+
          geom_line()+
          geom_point()+
          ggtitle("Season Onset and End by year")
        plot_grid(p1,p2, nrow=2, align="v")
        
       # widespread days vs seasAvgPrecip
        temp<-somTime[,1:18]
        temp<-merge(temp,seasStats[,c(1,4,6)],by="year")
        temp<-subset(temp, activityCat=="Widespread")
        table(temp$codes,temp$anomName)
        # just widespread days by year
        temp<-subset(temp, activityCat=="Widespread")
        countDOY<-temp %>% group_by(year) %>% count(codes)
        p1<-ggplot(countDOY, aes(fill=codes, y=n, x=year)) + 
          geom_bar(position="stack", stat="identity")+
          #geom_bar(position="fill", stat="identity")+
          scale_fill_brewer(type = "qual",
                            palette = "Paired",
                            direction = 1,
                            aesthetics = "fill")+
          ggtitle("Widespread-nodes frequency by year through season")
         # widespread nodes vs seasonal precip
          countDOY<-spread(countDOY,codes,n)
          countDOY<-cbind.data.frame(countDOY,seasStats[,c("avgPrecip")])
            colnames(countDOY)[length(countDOY)]<-"avgPrecip"
            countDOY[is.na(countDOY)] <- 0
            library("PerformanceAnalytics")
            #chart.Correlation(somTime[,c("percExtent","percExtent5","percExtent10","maxPrecip","meanPrecip","medPrecip")], histogram=TRUE, pch=19)
            chart.Correlation(countDOY[,c(2:length(countDOY))],method="spearman")  
        
# END OF SEASONAL TIMING/LENGTH        
#####
      

# contribution of days to seasonal totals - from SeasStats
  temp<- somTime %>% group_by(year) %>%
                  mutate(totalMeans=sum(meanPrecip))
  temp<-merge(temp,seasStats[,c(1,4,6)],by="year")
  #meanSumSeas<-merge(meanSumSeas, seasStats[,c("year","avgPrecip")],by="year")            
  temp$seasProp<-temp$meanPrecip/temp$totalMeans          
    # find day with max contrib/year
      maxContrib<-temp %>% group_by(year) %>%
                  filter(seasProp == max(seasProp, na.rm=TRUE))
      ggplot(maxContrib, aes(codes,seasProp))+
        geom_boxplot(varwidth = TRUE)+
        ggtitle("Contribution to seasonal total precip by node")
      ggplot(maxContrib, aes(year,seasProp, fill=codes))+
        geom_bar(stat = "identity")+
        geom_hline(yintercept = 1/108)+
        ggtitle("Contribution to seasonal total precip by node")
      # average daily contribution by node
      nodeContrib<-temp %>% group_by(codes) %>%
                    summarize(avg=mean(seasProp),
                              med=median(seasProp),
                              max=max(seasProp),
                              count=n(),
                              activity=first(activityCat))
      nodeContrib$avgAnn<-nodeContrib$avg*(nodeContrib$count/40)
      nodeContrib$avgCt<-nodeContrib$count/40
      
      nodeContrib <- nodeContrib %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      nodeContrib$row<-as.numeric(nodeContrib$row); nodeContrib$col<-as.numeric(nodeContrib$col); 
      p2<-ggplot(nodeContrib, aes(x=col,y=-row))+
        geom_tile(aes(fill = avgAnn*100), size=1.5)+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        #geom_text(aes(label = paste0(round(avgAnn*100,1)," (",round(avgCt,1),")")), size=6)+
        geom_text(aes(label = paste0(round(avgAnn*100,1))), size=6)+
        scale_fill_gradient(low = "#f7fcf5", high = "#006d2c", na.value = NA, name="%")+
        #scale_fill_gradient2(low = "blue",mid="grey", midpoint = 6.25,
        #                     high = "red", na.value = NA, name="% contribution")+
        ggtitle("  Avg total % contrib to seasonal precip")+
        #ggtitle("Median daily % contrib to seasonal total precip")+
        theme_bw()+
        scale_y_continuous(breaks=c(-1,-2,-3,-4),
                           labels=c("1", "2", "3","4"))+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank())+
        theme(legend.position="bottom")
      
      p1<-ggplot(nodeContrib, aes(x=col,y=-row))+
        geom_tile(aes(fill = avg*100), size=1.5)+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        #geom_text(aes(label = paste0(round(avg*100,1)," (",round(avgCt,1),")")), size=6)+
        geom_text(aes(label = paste0(round(avg*100,1))), size=6)+
        scale_fill_gradient(low = "#f7fbff", high = "#2171b5", na.value = NA, name="%")+
        #scale_fill_gradient2(low = "blue",mid="grey", midpoint = 6.25,
        #                     high = "red", na.value = NA, name="% contribution")+
        ggtitle("  Avg daily % contrib to seasonal total precip")+
        #ggtitle("Median daily % contrib to seasonal total precip")+
        theme_bw()+
        scale_y_continuous(breaks=c(-1,-2,-3,-4),
                         labels=c("1", "2", "3","4"))+
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank())+
        theme(legend.position="bottom")
      #plot_grid(p1,p2, ncol = 2, align = "h", labels=c("(a)","(b)"))
      
      # add in activity lines
      #gb <- ggplot_build(p1) 
      p1<- p1 + geom_segment(
                        aes(x=0.5, xend=1.5, y=-1.5, yend=-1.5), 
                        color="black")+
           geom_segment(
                        aes(x=1.5, xend=1.5, y=-0.5, yend=-1.5), 
                        color="black")+
           geom_segment( 
                     aes(x=2.5, xend=2.5, y=-2.5, yend=-0.5), 
                     color="black")+
           geom_segment( 
                     aes(x=0.5, xend=2.5, y=-2.5, yend=-2.5), 
                     color="black")+
           geom_segment( 
                     aes(x=2.5, xend=2.5, y=-3.5, yend=-4.5), 
                     color="black")+
           geom_segment( 
                     aes(x=3.5, xend=3.5, y=-2.5, yend=-3.5), 
                     color="black")+
           geom_segment( 
                     aes(x=2.5, xend=3.5, y=-3.5, yend=-3.5), 
                     color="black")+
           geom_segment( 
                     aes(x=3.5, xend=4.5, y=-2.5, yend=-2.5), 
                     color="black")
 
      p2<- p2 + geom_segment(
              aes(x=0.5, xend=1.5, y=-1.5, yend=-1.5), 
              color="black")+
              geom_segment(
                aes(x=1.5, xend=1.5, y=-0.5, yend=-1.5), 
                color="black")+
              geom_segment( 
                aes(x=2.5, xend=2.5, y=-2.5, yend=-0.5), 
                color="black")+
              geom_segment( 
                aes(x=0.5, xend=2.5, y=-2.5, yend=-2.5), 
                color="black")+
              geom_segment( 
                aes(x=2.5, xend=2.5, y=-3.5, yend=-4.5), 
                color="black")+
              geom_segment( 
                aes(x=3.5, xend=3.5, y=-2.5, yend=-3.5), 
                color="black")+
              geom_segment( 
                aes(x=2.5, xend=3.5, y=-3.5, yend=-3.5), 
                color="black")+
              geom_segment( 
                aes(x=3.5, xend=4.5, y=-2.5, yend=-2.5), 
                color="black")           
      ## FIGURE 4??
      plot_grid(p1,p2, ncol = 2, align = "h", labels=c("(a)","(b)"))  
      
      
      
      # seasonal contribution
      nodeContrib<-temp %>% group_by(codes,year) %>%
                    summarize(sumProp=sum(seasProp),
                              count=n())
      nodeContrib<-nodeContrib %>% group_by(codes) %>% 
                    summarize(avgSeasProp=mean(sumProp))
      nodeContrib <- nodeContrib %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      nodeContrib$row<-as.numeric(nodeContrib$row); nodeContrib$col<-as.numeric(nodeContrib$col); 
      ggplot(nodeContrib, aes(x=col,y=-row))+
        geom_tile(aes(fill = avgSeasProp*100))+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        geom_text(aes(label = round(avgSeasProp*100,1)), size=6)+
        #scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="% contribution")+
        scale_fill_gradient2(low = "lightblue",mid="grey", midpoint = 1,
                             high = "red", na.value = NA, name="% contribution")+
        ggtitle("Avg % contrib to seasonal total precip (based on annual avg)")+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
      # TERCILES -- average daily contribution by node
      nodeContrib<-temp %>% group_by(codes,anomName) %>%
        summarize(avg=mean(seasProp),
                  med=median(seasProp),
                  max=max(seasProp),
                  count=n(),
                  activity=first(activityCat))
      nodeContrib$anomYrs<-ifelse(nodeContrib$anomName=="wet",14,13)
      nodeContrib$avgAnn<-nodeContrib$avg*(nodeContrib$count/nodeContrib$anomYrs)
      nodeContrib <- nodeContrib %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      nodeContrib$row<-as.numeric(nodeContrib$row); nodeContrib$col<-as.numeric(nodeContrib$col); 
      ggplot(nodeContrib, aes(x=col,y=-row))+
        geom_tile(aes(fill = avg*100,color=activity), size=3)+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        geom_text(aes(label = round(avg*100,1)), size=6)+
        scale_colour_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a","#33a02c"))+
        #scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="% contribution")+
        scale_fill_gradient2(low = "blue",mid="grey", midpoint = 3,
                             high = "red", na.value = NA, name="% contribution")+
        ggtitle("Avg % contrib to seasonal total precip (based on avg rate*avg days)")+
        facet_wrap(.~anomName)+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
      # by activity class (KEY STATS FOR PAPER)
      nodeContrib<-temp %>% group_by(activityCat) %>%
        summarize(avg=mean(seasProp),
                  med=median(seasProp),
                  max=max(seasProp),
                  count=n())
      nodeContrib$avgAnn<-nodeContrib$avg*(nodeContrib$count/40)
      nodeContrib$avgCt<-nodeContrib$count/40
      ggplot(nodeContrib, aes(x=activityCat, y=1))+
        geom_tile(aes(fill = avgAnn*100))+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        geom_text(aes(label = paste0(round(avgAnn*100,1)," (",round(avgCt,1)," days)")), size=6)+
        #scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="% contribution")+
        scale_fill_gradient2(low = "lightblue",mid="grey", midpoint = 1,
                             high = "red", na.value = NA, name="% contribution")+
        ggtitle("Avg % contrib to seasonal total precip (based on avg rate*avg days)")+
        theme_bw()+
        theme(
              axis.text.y=element_blank())
      
      # USING MEAN DAILY PRECIP - contribution of days to seasonal totals - from SeasStats
      temp<- somTime %>% group_by(year) %>%
        mutate(totalMeans=sum(meanPrecip))
      temp<-merge(temp,seasStats[,c(1,4,6)],by="year")
      #meanSumSeas<-merge(meanSumSeas, seasStats[,c("year","avgPrecip")],by="year")            
      temp$seasProp<-temp$meanPrecip/temp$totalMeans          
      # find day with max contrib/year
      maxContrib<-temp %>% group_by(year) %>%
        filter(meanPrecip == max(meanPrecip, na.rm=TRUE))
      ggplot(maxContrib, aes(codes,meanPrecip))+
        geom_boxplot(varwidth = TRUE)+
        ggtitle("Contribution to seasonal total precip by node")
      ggplot(maxContrib, aes(year,meanPrecip, fill=codes))+
        geom_bar(stat = "identity")+
        #geom_hline(yintercept = 1/108)+
        ggtitle("Max daily mean precip by node")
      # average daily contribution by node
      nodeContrib<-temp %>% group_by(codes) %>%
        summarize(avg=mean(meanPrecip),
                  med=median(meanPrecip),
                  max=max(meanPrecip),
                  count=n())
      nodeContrib$avgAnn<-nodeContrib$avg*(nodeContrib$count/40)
      nodeContrib$avgCt<-nodeContrib$count/40
      nodeContrib <- nodeContrib %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      nodeContrib$row<-as.numeric(nodeContrib$row); nodeContrib$col<-as.numeric(nodeContrib$col); 
      ggplot(nodeContrib, aes(x=col,y=-row))+
        geom_tile(aes(fill = avgAnn))+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        geom_text(aes(label = paste0(round(avgAnn,1)," (",round(avgCt,1),")")), size=6)+
        #scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="% contribution")+
        scale_fill_gradient2(low = "blue",mid="grey", midpoint = 5,
                             high = "red", na.value = NA, name="contribution in mm")+
        ggtitle("Avg seasonal total precip in mm (based on avg rate*avg days)")+
        #ggtitle("Avg daily max precip")+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
      # TERCILES -- average daily contribution by node
      nodeContrib<-temp %>% group_by(codes,anomName) %>%
        summarize(avg=mean(meanPrecip),
                  med=median(meanPrecip),
                  max=max(meanPrecip),
                  count=n())
      nodeContrib$anomYrs<-ifelse(nodeContrib$anomName=="wet",14,13)
      nodeContrib$avgAnn<-nodeContrib$avg*(nodeContrib$count/nodeContrib$anomYrs)
      nodeContrib$avgCt<-nodeContrib$count/nodeContrib$anomYrs
      nodeContrib <- nodeContrib %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      nodeContrib$row<-as.numeric(nodeContrib$row); nodeContrib$col<-as.numeric(nodeContrib$col); 
        # add in new label
        anomLab<- nodeContrib %>% group_by(anomName) %>% summarize(sum=sum(avgAnn))
        nodeContrib<-merge(nodeContrib, anomLab, by="anomName")
        nodeContrib$anomLab<-paste0(nodeContrib$anomName," (seas avg: ",round(nodeContrib$sum,0)," mm)")
      ggplot(nodeContrib, aes(x=col,y=-row))+
        geom_tile(aes(fill = avg))+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        geom_text(aes(label = paste0(round(avgAnn,0)," (",round(avgCt,1),")")), size=4)+
        scale_fill_gradient(low = "#ffeda0", high = "#f03b20", na.value = NA, name="contribution (mm)")+
        #scale_fill_gradient2(low = "blue",mid="yellow", midpoint = 7.5,
        #                     high = "red", na.value = NA, name="contribution (mm)")+
        ggtitle("Avg contrib to seasonal total precip in mm (days)")+
        facet_wrap(.~anomLab)+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
      # by activity class
      nodeContrib<-temp %>% group_by(activityCat,anomName) %>%
        summarize(avg=mean(meanPrecip),
                  med=median(meanPrecip),
                  max=max(meanPrecip),
                  count=n())
      nodeContrib$anomYrs<-ifelse(nodeContrib$anomName=="wet",14,13)
      nodeContrib$avgAnn<-nodeContrib$avg*(nodeContrib$count/nodeContrib$anomYrs)
      nodeContrib$avgCt<-nodeContrib$count/nodeContrib$anomYrs
        # add in new label
        anomLab<- nodeContrib %>% group_by(anomName) %>% summarize(sum=sum(avgAnn))
        nodeContrib<-merge(nodeContrib, anomLab, by="anomName")
        nodeContrib$anomLab<-paste0(nodeContrib$anomName," (seas avg: ",round(nodeContrib$sum,0)," mm)")
      ggplot(nodeContrib, aes(x=activityCat, y=1))+
        geom_tile(aes(fill = avgAnn))+
        #geom_text(aes(label = codes), size=6)+
        #geom_text(aes(label = n), size=6)+
        geom_text(aes(label = paste0(round(avgAnn,1),"\n (",round(avgCt,1),")")), size=6)+
        scale_fill_gradient(low = "#ffeda0", high = "#f03b20", na.value = NA, name="contribution (mm)")+
        #scale_fill_gradient2(low = "lightblue",mid="grey", midpoint = 1,
        #                     high = "red", na.value = NA, name="contribution (mm)")+
        ggtitle("Avg % contrib to seasonal total precip (based on avg rate*avg days)")+
        facet_wrap(.~anomLab)+
        theme_bw()+
        theme(
          axis.text.y=element_blank())
      # USING MEAN PRECIP FOR SELECT COMPARISON YEARS
      temp<- somTime %>% group_by(year) %>%
        mutate(totalMeans=sum(meanPrecip))
      temp<-merge(temp,seasStats[,c(1,4,6)],by="year")
      #meanSumSeas<-merge(meanSumSeas, seasStats[,c("year","avgPrecip")],by="year")            
      temp$seasProp<-temp$meanPrecip/temp$totalMeans
      # SELECT COMPARISON YEARS
      temp<-subset(temp, year %in% c(1984,2020))
      # average daily contribution by node
      nodeContrib<-temp %>% group_by(codes,year) %>%
        summarize(avg=mean(meanPrecip),
                  med=median(meanPrecip),
                  max=max(meanPrecip),
                  count=n())
      nodeContrib$avgAnn<-nodeContrib$avg*(nodeContrib$count)
      nodeContrib$avgCt<-nodeContrib$count
      nodeContrib <- nodeContrib %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      nodeContrib$row<-as.numeric(nodeContrib$row); nodeContrib$col<-as.numeric(nodeContrib$col); 
        # add in new label
        yearLab<- nodeContrib %>% group_by(year) %>% summarize(sum=sum(avgAnn))
          nodeContrib<-merge(nodeContrib, yearLab, by="year")
          nodeContrib$yearLab<-paste0(nodeContrib$year," (seas avg: ",round(nodeContrib$sum,0)," mm)")
          ggplot(nodeContrib, aes(x=col,y=-row))+
            geom_tile(aes(fill = avgAnn))+
            #geom_text(aes(label = codes), size=6)+
            #geom_text(aes(label = n), size=6)+
            geom_text(aes(label = paste0(round(avgAnn,0)," (",round(avgCt,1),")")), size=4)+
            scale_fill_gradient(low = "#ffeda0", high = "#f03b20", na.value = NA, name="contribution (mm)")+
            #scale_fill_gradient2(low = "blue",mid="yellow", midpoint = 7.5,
            #                     high = "red", na.value = NA, name="contribution (mm)")+
            ggtitle("Avg contrib to seasonal total precip in mm (days)")+
            facet_wrap(.~yearLab)+
            theme_bw()+
            theme(axis.text.x=element_blank(),
                  axis.text.y=element_blank())
      
      # heat map of node occurrence for select year
          temp<-subset(somTime, year==2020)
          temp$codes<-as.factor(temp$codes)
          temp$codes<-factor(temp$codes, levels=c("1_1","1_2","2_2","2_1","1_4","1_3","2_4","2_3","3_3","3_2","3_1","4_2","4_1","3_4","4_3","4_4"))
          
          ggplot(temp, aes(date,codes, fill=activityCat))+
            geom_point(shape=22, size=2)+
            scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a","#33a02c"), drop=FALSE)+
            scale_y_discrete(drop=FALSE)+
            theme_bw()+
            ggtitle("Daily pattern classification - 2020")
          
        # heat map of all patterns by day
         
           nodeCount<-somTime %>% group_by(doy,codes) %>%
            summarize(count=n(),
                      activity=first(activityCat))
           nodeCount$codes<-as.factor(nodeCount$codes)
            nodeCount$codes<-factor(nodeCount$codes, levels=c("1_1","1_2","2_2","2_1","1_4","1_3","2_4","2_3","3_3","3_2","3_1","4_2","4_1","3_4","4_3","4_4"))
           nodeCount$activity<-as.factor(nodeCount$activity)
            nodeCount$activity<-factor(nodeCount$activity, levels=c("Widespread","Active-high","Active-low","Inactive"))
           # add in dates
           nodeCount$date<-as.Date(nodeCount$doy, origin = "2000-01-01")
           nodeCount$prop<-(nodeCount$count/40)*100
           
           # ggplot(nodeCount, aes(x=doy,y=codes))+
           #   geom_tile(aes(fill = count))+
           #   scale_fill_gradient2(low="#c7e9b4", mid="#41b6c4", high="#225ea8", midpoint = 5,
           #                        na.value = NA, name="count",limits=c(0, 10), oob=squish)+
           #   theme_bw()+
           #   ggtitle("Pattern frequency by day of year")
           
          p2<- ggplot(nodeCount, aes(date,codes, color=prop))+
             geom_point(shape=15, size=2)+
             scale_color_gradient2(low="#ffeda0", mid="#fd8d3c", high="#800026", midpoint = 12.5,
                                  na.value = NA,limits=c(2.5, 25), oob=squish,
                                  name="Frequency (% of days)",
                                  labels = c("5", "10", "15", "20",">25"),
                                  breaks = c(5,10,15,20,25))+
             ylab("node/pattern")+
             xlab("day of year")+
             theme_bw()+
             theme(panel.grid.major.y = element_blank())+
             ggtitle("Pattern frequency by day of year")+
             facet_grid(activity~.,scales = "free", space = "free")+
             theme(strip.background = element_blank(), #remove background for facet labels
                   panel.border = element_rect(colour = "black", fill = NA), #add black border
                   panel.spacing = unit(0, "lines"),
                   strip.text.y.right = element_text(angle = 0),
                   legend.position="bottom") 
           
           # add to activity frequency plot
           # activity chart 
           countDOY<-somTime %>% group_by(dummyDate) %>% count(activityCat)
           # moving average smooth
           countDOY<- countDOY %>%
             group_by(activityCat) %>%
             mutate(avg5 = zoo::rollapplyr(
               data = n,
               width = 10,
               FUN = mean,
               by.column = FALSE,
               fill = NA, align="center"))
           countDOY$avg5 <- ifelse(is.na(countDOY$avg5), countDOY$n, countDOY$avg5)
           #countDOY$avg5 <- countDOY$avg5*100
           
           # ggplot(countDOY, aes(fill=as.factor(activityCat), y=n, x=dummyDate)) + 
           #   #geom_bar(position="stack", stat="identity")
           #   geom_bar(position="fill", stat="identity")+
           #   scale_fill_brewer(type = "qual",
           #                     palette = "Paired",
           #                     direction = 1,
           #                     aesthetics = "fill")+
           #   ggtitle("Activity frequency by day through season")
           
           
           p1<-ggplot(countDOY, aes(fill=as.factor(activityCat), y=n, x=dummyDate)) + 
             #geom_bar(position="stack", stat="identity")
             geom_bar(position="fill", stat="identity", width=1.1)+
             scale_fill_brewer(type = "qual",
                               palette = "Paired",
                               direction = 1,
                               aesthetics = "fill", name="Activity class")+
             geom_hline(yintercept = 0.5)+
             ggtitle("Activity frequency by day through season")+
             ylab("Proportion of days")+
             xlab("")+
             theme_bw()
           # FIGURE 3 ???
           plot_grid(p1,p2, ncol = 1, align = "v", axis = 'lr')
           
          
      # closer look at 3_4 and 4_4
        temp<- somTime %>% group_by(year) %>%
            mutate(totalMeans=sum(meanPrecip))
        temp<-merge(temp,seasStats[,c(1,4,6)],by="year")
        #meanSumSeas<-merge(meanSumSeas, seasStats[,c("year","avgPrecip")],by="year")            
        temp$seasProp<-temp$meanPrecip/temp$totalMeans
        # SELECT COMPARISON YEARS
        temp<-subset(temp, codes %in% c("1_4","4_4"))
        nodeContrib<-temp %>% group_by(codes,year) %>%
                      summarize(avg=mean(meanPrecip),
                                med=median(meanPrecip),
                                max=max(meanPrecip),
                                sumPrecip=sum(meanPrecip),
                                sumProp=sum(seasProp),
                                count=n())  
        ggplot(nodeContrib, aes(year,count, fill=codes))+
          geom_bar(stat = "identity", color="grey")+  # position="dodge",
          #geom_hline(yintercept = 1/108)+
          ggtitle("Frequency of nodes 1_4 and 4_4")+
          facet_wrap(.~codes, nrow=2)+
          geom_quantile(quantiles = 0.5)
        ggplot(nodeContrib, aes(year,sumPrecip, fill=codes))+
          geom_bar(stat = "identity", color="grey")+  # position="dodge",
          #geom_hline(yintercept = 1/108)+
          ggtitle("Total precip seas from nodes 1_4 and 4_4")+
          facet_wrap(.~codes, nrow=2)+
          geom_quantile(quantiles = 0.5)+
          geom_smooth(method = "lm")

        # closer look at all widespread days
        temp<- somTime %>% group_by(year) %>%
          mutate(totalMeans=sum(meanPrecip))
        temp<-merge(temp,seasStats[,c(1,4,6)],by="year")
        #meanSumSeas<-merge(meanSumSeas, seasStats[,c("year","avgPrecip")],by="year")            
        temp$seasProp<-temp$meanPrecip/temp$totalMeans
        # SELECT COMPARISON YEARS
        temp<-subset(temp, activityCat=="Widespread")
        nodeContrib<-temp %>% group_by(codes,year) %>%
          summarize(avg=mean(meanPrecip),
                    med=median(meanPrecip),
                    max=max(meanPrecip),
                    sumPrecip=sum(meanPrecip),
                    sumProp=sum(seasProp),
                    count=n())  
        ggplot(nodeContrib, aes(year,sumPrecip, fill=codes))+
          geom_bar(stat = "identity", color="grey")+  # position="dodge",
          #geom_hline(yintercept = 1/108)+
          ggtitle("Frequency of nodes")+
          facet_wrap(.~codes, nrow=2)+
          geom_quantile(quantiles = 0.5)+
          geom_smooth(method = "lm")
        
        
  ##### 
  # get seasonal average map
  totLayer<-stackApply(subLayers, somTime$year, fun=sum)
      seasAvg<-calc(totLayer, mean)
      
  # average total and percent contribution to seasonal total
      comp<-stackApply(subLayers, somTime$mapUnit, fun=sum)
      tempNames<-as.data.frame(names(comp))
      colnames(tempNames)<-"index"
      tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
      comp<-subset(comp, order(c(as.numeric(tempNames$mapUnit))))
      names(comp)<-codeList
      avgTotNode<-comp/length(seasAvgPrecip$anomName)
      percNodeSeas<-(comp/calc(comp, sum, na.rm=TRUE))*100
      names(percNodeSeas)<-codeList    
      
  # Composites of seas totals in given years
  yr<-1984
  temp<-subLayers[[which(somTime$year==yr)]]
  comp<-stackApply(temp, somTime$mapUnit[which(somTime$year==yr)], fun=sum)
    tempNames<-as.data.frame(names(comp))
    colnames(tempNames)<-"index"
    tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
    # find missing
    emptyNodes<-setdiff(seq(from=1,to=nrows*ncols, by=1), (as.numeric(tempNames$mapUnit)))
    emptyLyr<-comp[[1]]; emptyLyr[emptyLyr >= 0] <- NA
    if(length(emptyNodes)!=0){
      comp <- stack(comp, stack(replicate(length(emptyNodes), emptyLyr)))      
    }
    comp<-subset(comp, order(c(as.numeric(tempNames$mapUnit),emptyNodes)))
    names(comp)<-codeList
    
    # add in counts of days in node
    tempCode<-as.data.frame(table(somTime$mapUnit[which(somTime$year==yr)]))
      tempCode<-merge(code_grid,tempCode,by.x="mapUnit",by.y="Var1", all=TRUE) 
      tempCode<-cbind.data.frame(tempCode,as.data.frame(round(table(somTime$codes)/length(unique(somTime$year)),1)))
        colnames(tempCode)[5]<-"avg"
      tempCode$anom<-round(tempCode$Freq-tempCode$avg,1)
      #col.titles<-paste0("",tempCode$codes,", ",tempCode$Freq," (",tempCode$anom,")")  
      col.titles<-paste0(tempCode$Freq," (",tempCode$anom,")")  
    
    # plot composite total 
    #at<-c(seq(0,100,5),250)
    at<-c(seq(0,250,5))
    
    mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
    pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at, names.attr=col.titles, 
                     main=paste0("Composite Precipitation ",yr), xlab=NULL, ylab=NULL)+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    # plot composite total anomaly
    anom<-comp-avgTotNode
    at<-c(-100,seq(-50,50,10),100)
    mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
    pComp<-levelplot(comp-avgTotNode, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                     main=paste0("Composite Precipitation Anomaly ",yr), xlab=NULL, ylab=NULL)+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    
    # percent of seasonal total
    percNode<-(comp/calc(comp, sum, na.rm=TRUE))*100
      names(percNode)<-codeList
    at<-seq(0,100,2.5)
    mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pPerc<-levelplot(percNode, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                       main=paste0("Percent of seasonal total by node: ",yr), xlab=NULL, ylab=NULL)+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    plot(percNode-percNodeSeas)
      
    # seasonal total
      at<-c(seq(0,500,25),1000)
      total<-calc(comp, sum, na.rm=TRUE)
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pTotal<-levelplot(total, contour=FALSE, margin=FALSE, par.settings=mapTheme, at=at,  
                       main=paste0("Seasonal total: ",yr), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    # percent of average
      percAvg<-(total/seasAvg)*100
      at<-c(seq(0,200,10),500)
      mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
      pAnom<-levelplot(percAvg, contour=FALSE, margin=FALSE, par.settings=mapTheme, at=at,  
                        main=paste0("Seasonal total: ",yr), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    # anomaly
      anomSeas<-(total-seasAvg)
      at<-c(seq(-200,200,10))
      mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
      pAnom<-levelplot(anomSeas, contour=FALSE, margin=FALSE, par.settings=mapTheme, at=at,  
                       main=paste0("Seasonal total Anomaly: ",yr), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))  
          # check against node version
          # at<-c(seq(-200,200,10))
          # mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
          # pAnom<-levelplot(calc(comp-avgTotNode,sum), contour=FALSE, margin=FALSE, par.settings=mapTheme, at=at,  
          #                  main=paste0("Seasonal total Anomaly: ",yr), xlab=NULL, ylab=NULL)+
          #   layer(sp.polygons(aznm, col = 'gray40', lwd=1))  
          # plot((total-seasAvg)-calc(comp-avgTotNode,sum))
      
    # normalized contribution
      # normalize by count?
      counts<-as.vector(table(somTime$codes[which(somTime$year==yr)]))
      norm<-comp/counts
      names(norm)<-codeList
      at<-seq(0,100,2.5)
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pNorm<-levelplot(norm, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                       main=paste0("Average total by node: ",yr), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      
      
    # compare two seasonal totals
      seasTotals<-stackApply(subLayers, somTime$year, fun=sum)
      at<-seq(0,100,2.5)
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pTot<-levelplot(stack(seasTotals[[c(16,32)]]), contour=FALSE, margin=FALSE, par.settings=mapTheme, 
                       main=paste0("Percent of seasonal total by node: ",yr), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      
      mapTheme<-rasterTheme(region = c("purple", "white", "orange"))
      at<-seq(-600,600,30)
      pDiff<-levelplot(seasTotals[[32]]-seasTotals[[16]], contour=FALSE, margin=FALSE, par.settings=mapTheme, at=at, 
                      main="Seasonal total: 2012-1996", xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      
      
  #####
  
      # average seasonal total by node
      
      
      
      # average total and percent contribution to seasonal total
      comp<-stackApply(subLayers, somTime$mapUnit, fun=sum)
      tempNames<-as.data.frame(names(comp))
      colnames(tempNames)<-"index"
      tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
      comp<-subset(comp, order(c(as.numeric(tempNames$mapUnit))))
      names(comp)<-codeList
        avgNode<-comp/length(seasAvgPrecip$anomName)
        percNode<-(comp/calc(comp, sum, na.rm=TRUE))*100
      names(percNode)<-codeList
      # normalize by count?
      counts<-as.vector(table(somTime$mapUnit))
      norm<-comp/counts
      # 
      at<-seq(0,50,1)
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      #at<-c(seq(0,24,1),50)
      #mapTheme<-rasterTheme(region = c("red","white","blue"))
      pPerc<-levelplot(percNode, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                       main=paste0("Percent of seasonal total by node"), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      at<-seq(0,80,1)
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pAvg<-levelplot(avgNode, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,   
                       main=paste0("Avg seasonal total by node (mm)"), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      # total JAS precip
      total<-calc(avgNode, sum)
      pTotal<-levelplot(total, contour=FALSE, margin=FALSE, par.settings=mapTheme,    
                      main=paste0("Avg seasonal total for region (mm)"), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      # norm JAS precip
      at<-c(seq(0,30,0.3))
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pNorm<-levelplot(norm, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,   
                      main=paste0("Avg precip by node (mm)"), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
            
      #test<-calc(percNode, sum, na.rm=TRUE)
  # daily average percent contribution
  
      
  # wet/norm/dry proportions    
      anomName<-"wet"
      #anomYrs<-seasAvgPrecip$`unique(subDates$year)`[which(seasAvgPrecip$anomName==anomName)]
      anomYrs<-seasAvgPrecip$year[which(seasAvgPrecip$anomName==anomName)]
      temp<-subLayers[[which(somTime$year %in% anomYrs)]] 
      comp<-stackApply(temp, somTime$mapUnit[which(somTime$year %in% anomYrs)], fun=sum)  
      tempNames<-as.data.frame(names(comp))
      colnames(tempNames)<-"index"
      tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
      # find missing
      emptyNodes<-setdiff(seq(from=1,to=nrows*ncols, by=1), (as.numeric(tempNames$mapUnit)))
      emptyLyr<-comp[[1]]; emptyLyr[emptyLyr >= 0] <- NA
      if(length(emptyNodes)!=0){
        comp <- stack(comp, stack(replicate(length(emptyNodes), emptyLyr)))      
      }
      comp<-subset(comp, order(c(as.numeric(tempNames$mapUnit),emptyNodes)))
      names(comp)<-codeList
      # normalize by count?
      counts<-as.vector(table(somTime$mapUnit[which(somTime$year %in% anomYrs)]))
      norm<-comp/counts
      # plot 
      at<-c(seq(0,125,5))
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pComp<-levelplot(comp/length(anomYrs), contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                       main=paste0("Average Composite Precipitation: ",anomName), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      # percent of seasonal total
      percNode<-(comp/calc(comp, sum, na.rm=TRUE))*100
      names(percNode)<-codeList
      at<-seq(0,50,1)
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pPerc<-levelplot(percNode, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                       main=paste0("Percent of seasonal total by node: ",anomName), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      # average precip per node
      at<-seq(0,35,0.5)
      mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pNorm<-levelplot(norm, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows),at=at,  
                       main=paste0("Average precip by node: ",anomName), xlab=NULL, ylab=NULL)+
        layer(sp.polygons(aznm, col = 'gray40', lwd=1))
      
      # diffs - wet/norm/dry proportions 
      anomName<-"dry"
        #anomYrs<-seasAvgPrecip$`unique(subDates$year)`[which(seasAvgPrecip$anomName==anomName)]
        anomYrs<-seasAvgPrecip$year[which(seasAvgPrecip$anomName==anomName)]
        temp<-subLayers[[which(somTime$year %in% anomYrs)]] 
        comp<-stackApply(temp, somTime$mapUnit[which(somTime$year %in% anomYrs)], fun=sum)  
        tempNames<-as.data.frame(names(comp))
        colnames(tempNames)<-"index"
        tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
        # find missing
        emptyNodes<-setdiff(seq(from=1,to=nrows*ncols, by=1), (as.numeric(tempNames$mapUnit)))
        emptyLyr<-comp[[1]]; emptyLyr[emptyLyr >= 0] <- NA
        if(length(emptyNodes)!=0){
          comp <- stack(comp, stack(replicate(length(emptyNodes), emptyLyr)))      
        }
        comp<-subset(comp, order(c(as.numeric(tempNames$mapUnit),emptyNodes)))
        names(comp)<-codeList 
        # percent of seasonal total - CHANGE NAME HERE
        #percNode_wet<-((comp/length(anomYrs))/calc(comp/length(anomYrs), sum, na.rm=TRUE))*100
        #names(percNode_wet)<-codeList
        # or just seasonal total
        percNode_dry<-(comp/length(anomYrs))
        names(percNode_dry)<-codeList
        
        ##
        at<-c(-60,seq(-30,30,2.5),60)
        pDiff<-levelplot(percNode_wet-percNode_dry, contour=FALSE, margin=FALSE, par.settings=RdBuTheme,layout=c(ncols,nrows),at=at,  
                         main=paste0("Average Seasonal total by node: Wet-Dry"), xlab=NULL, ylab=NULL)+
          layer(sp.polygons(aznm, col = 'gray40', lwd=1))
        
      
  ####    
      
      
      
  # Clustering of nodes
  ## Show the U matrix
  Umat <- plot(som.gh500, type="dist.neighbours", main = "SOM neighbour distances")
  ## use hierarchical clustering to cluster the codebook vectors
  som.hc <- cutree(hclust(object.distances(som.gh500, "codes")),8)
  add.cluster.boundaries(som.gh500, som.hc)
  
  # map clustered nodes ----
  
  #####
  
  # DIAGNOSTICS
  plot(som.gh500, type="changes") # changes, codes, counts, property, quality, mapping
  
  # summary plots -- appears to plot opposite up/down from SOM plot
  counts <- plot(som.gh500, type="counts", shape = "straight", labels=counts)
  codes <- plot(som.gh500, type="codes", shape = "straight")
  similarities <- plot(som.gh500, type="quality", palette.name = terrain.colors)
  plot(som.gh500, type="dist.neighbours", main = "SOM neighbour distances")
  plot(som.gh500)
  
# distributions of codebook values
  temp<-as.data.frame(t(som.gh500$codes[[1]]))
  colnames(temp)<-codeList
  temp<-gather(temp, codes, precip)
  temp<-merge(temp, activity, by="codes")
  ggplot(temp, aes(x=reorder(codes, precip), y=precip, fill=activityCat))+
    geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
    ggtitle("Distribution of precipitation values in each node")+
    theme(legend.position="bottom")+
    facet_wrap(~activityCat, scales="free_x",nrow = 1)+
    theme_bw()
  
  # sammon mapping
  library(MASS)
  gh500.codes <- som.gh500$codes
  dis <- dist(as.matrix(som.gh500$codes[[1]]))
  gh500.sam <- sammon(dis)
  plot(gh500.sam$points, type="n")
  text(gh500.sam$points,labels=codeList)
  ##  Polygon version of map ----
  library(sp)
    temp<-list()
    ctr<-1
    for(j in 1:(nrows-1)){
      for(k in 1:(ncols-1)){
        temp[ctr]<-Polygons(list(Polygon(cbind(c(gh500.sam$points[(k+(ncols*j)-ncols),1], gh500.sam$points[(k+(ncols*j)-ncols)+1,1],
                                           gh500.sam$points[k+(ncols*j)+1,1],gh500.sam$points[k+(ncols*j),1]),
                                         c(gh500.sam$points[(k+(ncols*j)-ncols),2], gh500.sam$points[(k+(ncols*j)-ncols)+1,2],
                                           gh500.sam$points[k+(ncols*j)+1,2],gh500.sam$points[k+(ncols*j),2])))),paste0(ctr))
        ctr<-ctr+1
      }
    }  
    sr1<-SpatialPolygons(temp)  
    plot(gh500.sam$points, type="n")
    #text(gh500.sam$points,labels=as.character(1:nrow(code_grid)))
    text(gh500.sam$points,labels=codeList)
    #text(gh500.sam$points,labels=activity$activityCat)
    plot(sr1, add=TRUE)
  # ----
    
    
    
    
    
  # GGPLOT versions
  # counts --- plotting incorrectly
  counts<-somTime %>% group_by(codes) %>% count(codes)
    counts <- counts %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    counts$row<-as.numeric(counts$row); counts$col<-as.numeric(counts$col); 
    counts$avg<-(counts$n/length(seq(1981,2019,1)))
pCt<-ggplot(counts, aes(x=col,y=-row))+
      geom_tile(aes(fill = (n/nrow(somTime)*100)))+
      #geom_text(aes(label = codes), size=6)+
      #geom_text(aes(label = n), size=6)+
      geom_text(aes(label = round((n/nrow(somTime)*100),0)), size=6)+
      scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="% days")+
      #scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 10,
      #                     high = "red", na.value = NA, name="% days")+
      ggtitle("Percent of days/node")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
pAvgDays<-ggplot(counts, aes(x=col,y=-row))+
  geom_tile(aes(fill = avg))+
  #geom_text(aes(label = codes), size=6)+
  #geom_text(aes(label = n), size=6)+
  geom_text(aes(label = round(avg,1)), size=6)+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="# days")+
  #scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 10,
  #                     high = "red", na.value = NA, name="% days")+
  ggtitle("Avg # of days/node")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())


  # SOM distances
    ndist <- unit.distances(som.gh500$grid)
    cddist <- as.matrix(object.distances(som.gh500, type = "codes"))
    cddist[abs(ndist - 1) > .001] <- NA
    neigh.dists <- colMeans(cddist, na.rm = TRUE)
    som_grid <- som.gh500[[4]]$pts %>%
      as_tibble %>% 
      mutate(id=row_number())
    som_grid$codes<-paste0(som_grid$y,"_",som_grid$x)
    som_grid <- som_grid %>% mutate(dist=neigh.dists)
pDist<-ggplot(som_grid, aes(x=x,y=-y))+
      geom_tile(aes(fill = dist))+
      geom_text(aes(label = codes), size=6)+
      scale_fill_gradient(low = "lightblue", high = "red", na.value = NA, name="distance")+
      ggtitle("Neighborhood Distance")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
# SOM Spearman corr by day
spear<-somTime %>% group_by(codes) %>% summarise(spearR = mean(spearman, na.rm=TRUE))
  spear <- spear %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  spear$row<-as.numeric(spear$row); spear$col<-as.numeric(spear$col);
  pRho<-ggplot(spear, aes(x=col,y=-row))+
    geom_tile(aes(fill = spearR))+
    geom_text(aes(label = round(spear$spearR,2)), size=6)+
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="mean Spearman")+
    ggtitle("Mean Spearman Corr by node")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())  
# boxplot of spearman values
  ggplot(somTime, aes(as.factor(somTime$codes), spearman))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Spearman Rho by nodes")
  ggplot(somTime, aes(as.factor(somTime$codes), pearson))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Pearson r by nodes") 
  ggplot(somTime, aes(as.factor(somTime$codes), kendall))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Kendall tau by nodes") 
  # facet wrapped
  ggplot(somTime, aes(as.factor(somTime$codes), spearman))+
    facet_wrap(~as.factor(somTime$codes),scales="free_x")+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Spearman Rho by nodes")
  
  # SOM node quality
    sim<-cbind.data.frame(som.gh500$unit.classif,som.gh500$distances)
      colnames(sim)<-c("node","distance")
    sim <- sim %>% group_by(node) %>% summarise(dist = mean(distance))
    sim$code<-codeList
    sim <- sim %>% separate(code, c("row","col"), sep = "_", remove=FALSE)
      sim$row<-as.numeric(sim$row); sim$col<-as.numeric(sim$col); 
  pQ<-ggplot(sim, aes(x=col,y=-row))+
        geom_tile(aes(fill = dist))+
        geom_text(aes(label = code), size=6)+
        scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="mean error")+
        ggtitle("SOM Node Quality")+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
  plot_grid(pCt,pDist,pQ,pRho, ncol = 1, align = "v")
  
  # relationships between SOM metrics
  ggplot(somTime, aes(x=percExtent10,y=maxPrecip))+
    geom_point()+
    facet_wrap(~as.factor(codes))+
    ggtitle("Daily Precip Extent vs Max Precip")
  ggplot(somTime, aes(x=meanPrecip,y=spearman))+
    geom_point()+
    facet_wrap(~as.factor(codes))+
    ggtitle("Daily Mean Precip vs Pearson")
    #xlim(0,5)
  ggplot(somTime, aes(x=percExtent,y=meanPrecip, color=maxPrecip))+
    geom_point()+
    facet_wrap(~as.factor(codes))+
    ggtitle("Daily Precip Extent vs Mean Precip")+
    scale_color_gradient(low = "blue", high = "red")
  ggplot(somTime, aes(as.factor(somTime$codes), meanPrecip))+
    facet_wrap(~as.factor(somTime$codes),scales="free_x")+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of mean daily precip by nodes")
  
  
  #
  
  # SOM code vectors
  
  # SOM code vectors distributions
  ggplot(codebook.long, aes(as.factor(codebook.long$codes), value))+
    geom_boxplot(varwidth = TRUE)+
    xlab("Node")+
    ylab("mm")+
    ggtitle("Distribution of Codebook Vectors (precip in mm)")
  
  temp<-merge(codebook.long, activity, by="codes")
  temp<-temp %>% 
    group_by(codes, activityCat) %>%
    summarize(n_thresh = sum(value>10),
              n_cells = n(),
              perc_cov=(n_thresh/n_cells)*100)
  ggplot(temp, aes(x=reorder(codes, perc_cov), y=perc_cov, fill=activityCat))+
    #geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
    geom_bar(stat="identity")+
    ggtitle("percent coverage of codebook values >10mm in each node")+
    theme(legend.position="bottom")+
    facet_wrap(~activityCat, scales="free_x",nrow = 1)+
    theme_bw()
  
  # SOM node counts by year - facet wrap indiv heat maps
  countsYr<-somTime %>% group_by(codes,year) %>% count(codes)
  countsYr <- countsYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  countsYr$row<-as.numeric(countsYr$row); countsYr$col<-as.numeric(countsYr$col); 
    ggplot(countsYr, aes(x=col, y=-row))+
      geom_tile(aes(fill=n))+
      geom_text(aes(label = n), size=4)+
      facet_wrap(~year)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 25,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Count of days in each node by year")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  # time series
  temp<-subset(countsYr, codes=="4_4")
  ggplot(temp, aes(x=year, y=n, color=as.factor(codes)))+
    geom_line()
  
  # SOM node anomaly
  countAnom<-merge(countsYr, counts, by="codes")
  countAnom$anomCT<-countAnom$n.x/countAnom$avg
  ggplot(countAnom, aes(x=col.x, y=-row.x))+
    geom_tile(aes(fill=anomCT))+
    #geom_text(aes(label = anomCT), size=4)+
    facet_wrap(~year)+
    scale_fill_gradient2(low = "purple",mid="white", midpoint = 1,
                         high = "orange", na.value = NA, name="ratio to avg",limits=c(0, 2), oob=squish)+
    ggtitle("Anom of days in each node by year")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  
  
  # SOM node counts by month
  countsMo<-somTime %>% group_by(codes,month) %>% count(codes)
  countsMo <- countsMo %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  countsMo$row<-as.numeric(countsMo$row); countsMo$col<-as.numeric(countsMo$col); 
    ggplot(countsMo, aes(x=col, y=-row))+
      geom_tile(aes(fill=n))+
      geom_text(aes(label = n), size=4)+
      facet_wrap(~month)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 350,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Count of days in each node by month")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
  # Expected Counts - SOM node counts by month
  temp2<-somTime %>% group_by(month) %>% count(month)
    temp2$prop<-temp2$n/sum(temp2$n)
  temp3<-somTime %>% group_by(codes) %>% count(codes)
  countsMo<-somTime %>% group_by(codes,month) %>% count(codes)
  countsMo <- countsMo %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  countsMo$row<-as.numeric(countsMo$row); countsMo$col<-as.numeric(countsMo$col); 
    countsMo<-merge(countsMo,temp2, by="month")
    countsMo$anom<-(countsMo$n.x/countsMo$n.y)*100
    countsMo<-merge(countsMo, temp3, by="codes")
    countsMo$expCt<-countsMo$n*countsMo$prop
    countsMo$anomCt<-countsMo$n.x-countsMo$expCt
    # proportion test prop.test()
    propPval<- lapply(seq_along(countsMo$n.x),
                      function(i) prop.test(countsMo$n.x[i],countsMo$n[i],p=countsMo$prop[i],alternative = "two.sided")$p.value)
    binomPval<- lapply(seq_along(countsMo$n.x),
                       function(i) binom.test(countsMo$n.x[i],countsMo$n[i],p=countsMo$prop[i],alternative = "two.sided")$p.value)
    countsMo$propPval<-unlist(propPval)
    countsMo$binomPval<-unlist(binomPval)
    countsMo$sig<-countsMo$propPval<=0.05
    countsMo$sigBinom<-countsMo$binomPval<=0.05
    countsMo$anomCtLabel<-ifelse(countsMo$sig==TRUE,
                                  paste0(round(countsMo$anomCt,1),"*"),
                                  paste0(round(countsMo$anomCt,1)))
  ggplot(countsMo, aes(x=col, y=-row))+
    geom_tile(aes(fill=anomCt))+
    geom_text(aes(label = anomCtLabel), size=4)+
    facet_wrap(~month)+
    scale_fill_gradient2(low = "orange",mid="grey", midpoint = 0,
                         high = "purple", na.value = NA, name="count", limits=c(-50, 50), oob=squish)+
    ggtitle("Anom of expected days in each node by Month")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  
  
  # SOM node counts by MJO phase - CHECK DAY OF YEAR ALIGNMENT; 
  temp<-subset(somTime, amplitude>=1)
  temp2<-temp %>% group_by(phase) %>% count(phase)
    temp2$prop<-temp2$n/sum(temp2$n)
  temp3<-temp %>% group_by(codes) %>% count(codes)
  countsMJO<-temp %>% group_by(codes,phase) %>% count(codes)
  countsMJO <- countsMJO %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  countsMJO$row<-as.numeric(countsMJO$row); countsMJO$col<-as.numeric(countsMJO$col); 
    countsMJO<-merge(countsMJO,temp2, by="phase")
    countsMJO$anom<-(countsMJO$n.x/countsMJO$n.y)*100
    countsMJO<-merge(countsMJO, temp3, by="codes")
    countsMJO$expCt<-countsMJO$n*countsMJO$prop
    countsMJO$anomCt<-countsMJO$n.x-countsMJO$expCt
  # proportion test prop.test()
    propPval<- lapply(seq_along(countsMJO$n.x),
                      function(i) prop.test(countsMJO$n.x[i],countsMJO$n[i],p=countsMJO$prop[i],alternative = "two.sided")$p.value)
    binomPval<- lapply(seq_along(countsMJO$n.x),
                       function(i) binom.test(countsMJO$n.x[i],countsMJO$n[i],p=countsMJO$prop[i],alternative = "two.sided")$p.value)
    countsMJO$propPval<-unlist(propPval)
    countsMJO$binomPval<-unlist(binomPval)
    countsMJO$sig<-countsMJO$propPval<=0.05
    countsMJO$sigBinom<-countsMJO$binomPval<=0.05
    countsMJO$anomCtLabel<-ifelse(countsMJO$sig==TRUE,
                                  paste0(round(countsMJO$anomCt,1),"*"),
                                  paste0(round(countsMJO$anomCt,1)))
  # plot 
  ggplot(countsMJO, aes(x=col, y=-row))+
    geom_tile(aes(fill=anomCt))+
    geom_text(aes(label = anomCtLabel), size=4)+
    facet_wrap(~phase, ncol = 4)+
    scale_fill_gradient2(low = "lightblue",mid="white", midpoint = 0,
                         high = "red", na.value = NA, name="count")+
    ggtitle("Anom of expected days in each node by BSISO (amp>1, *pval<0.05) Phase")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())

  # SOM node counts by ENSO phase
  temp2<-somTime %>% group_by(ENSO) %>% count(ENSO)
    temp2$prop<-temp2$n/sum(temp2$n)
    temp3<-somTime %>% group_by(codes) %>% count(codes)
  countsENSO<-somTime %>% group_by(codes,ENSO) %>% count(codes)
  countsENSO <- countsENSO %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  countsENSO$row<-as.numeric(countsENSO$row); countsENSO$col<-as.numeric(countsENSO$col); 
    countsENSO<-merge(countsENSO,temp2, by="ENSO")
    countsENSO$anom<-(countsENSO$n.x/countsENSO$n.y)*100
    countsENSO<-merge(countsENSO, temp3, by="codes")
    countsENSO$expCt<-countsENSO$n*countsENSO$prop
    countsENSO$anomCt<-countsENSO$n.x-countsENSO$expCt
  # proportion test prop.test()
    propPval<- lapply(seq_along(countsENSO$n.x),
                      function(i) prop.test(countsENSO$n.x[i],countsENSO$n[i],p=countsENSO$prop[i],alternative = "two.sided")$p.value)
    binomPval<- lapply(seq_along(countsENSO$n.x),
                       function(i) binom.test(countsENSO$n.x[i],countsENSO$n[i],p=countsENSO$prop[i],alternative = "two.sided")$p.value)
    countsENSO$propPval<-unlist(propPval)
    countsENSO$binomPval<-unlist(binomPval)
    countsENSO$sig<-countsENSO$propPval<=0.05
    countsENSO$sigBinom<-countsENSO$binomPval<=0.05
    countsENSO$anomCtLabel<-ifelse(countsENSO$sig==TRUE,
                                  paste0(round(countsENSO$anomCt,1),"*"),
                                  paste0(round(countsENSO$anomCt,1)))  
  # ggplot  
  ggplot(countsENSO, aes(x=col, y=-row))+
    geom_tile(aes(fill=anomCt))+
    geom_text(aes(label = anomCtLabel), size=4)+
    facet_wrap(~ENSO)+
    scale_fill_gradient2(low = "lightblue",mid="white", midpoint = 0,
                         high = "red", na.value = NA, name="count")+
    ggtitle("Anom of expected days in each node by ENSO Phase (*pval<0.05)")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  
  
  # classification error by node-year
  errorYr<-somTime %>% 
           group_by(codes,year) %>%
           summarise(error=mean(errorDist))
  errorYr <- errorYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  errorYr$row<-as.numeric(errorYr$row); errorYr$col<-as.numeric(errorYr$col); 
  ggplot(errorYr, aes(x=col, y=-row))+
    geom_tile(aes(fill=round(error,0)))+
    #geom_text(aes(label = error), size=2)+
    facet_wrap(~year)+
    scale_fill_gradient2(low = "blue",mid="yellow", 
                         high = "red", midpoint = 300000, 
                         na.value = NA, name="class. error")+
    ggtitle("Classification error by node-year")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
      # error by year
      errorYr<-somTime %>% 
        group_by(year) %>%
        summarise(error=median(errorDist))
      ggplot(errorYr, aes(x=year, y=error))+
        geom_bar(stat = "identity", fill="red")+
        ggtitle("Median error of all nodes by year")

  # classification error by node-year
    spearYr<-somTime %>% 
        group_by(codes,year) %>%
        summarise(error=mean(spearman))
      spearYr <- spearYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      spearYr$row<-as.numeric(spearYr$row); spearYr$col<-as.numeric(spearYr$col); 
      ggplot(spearYr, aes(x=col, y=-row))+
        geom_tile(aes(fill=error))+
        #geom_text(aes(label = error), size=2)+
        facet_wrap(~year)+
        scale_fill_gradient2(low = "blue",mid="yellow", 
                             high = "red", midpoint = 0.5, 
                             na.value = NA, name="mean rho")+
        ggtitle("Mean Spearman Corr by node-year")+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
      # error by year
      spearYr<-somTime %>% 
        group_by(year) %>%
        summarise(error=mean(spearman, na.rm=TRUE))
      ggplot(spearYr, aes(x=year, y=error))+
        geom_bar(stat = "identity", fill="red")+
        ggtitle("Mean Spearman Rho of all nodes by year")     
      
        
  # seasonal total precip by node
  seasPrecip<-somTime %>% group_by(year) %>% summarise(totPrecip=sum(meanPrecip))    
  sumYr<-somTime %>% 
    group_by(codes,year) %>%
    summarise(sumPrecip=sum(meanPrecip), countNodes=n())
  sumYr <- sumYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  sumYr$row<-as.numeric(sumYr$row); sumYr$col<-as.numeric(sumYr$col); 
    sumYr<-merge(sumYr,seasPrecip, by="year")
    sumYr$perc<-sumYr$sumPrecip/sumYr$totPrecip
  ggplot(sumYr, aes(x=col, y=-row))+
    geom_tile(aes(fill=round(sumPrecip,0)))+
    #geom_text(aes(label = error), size=2)+
    facet_wrap(~year)+
    scale_fill_gradient2(low = "blue",mid="yellow", 
                         high = "red", midpoint = 20, 
                         na.value = NA, name="Sum Precip")+
    ggtitle("Sum of Precip by node-year")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  # % of seasonal precip by node
  ggplot(sumYr, aes(x=col, y=-row))+
    geom_tile(aes(fill=round(perc*100,1)))+
    #geom_text(aes(label = error), size=2)+
    facet_wrap(~year)+
    scale_fill_gradient2(low = "blue",mid="yellow", 
                         high = "red", midpoint = 15, 
                         na.value = NA, name="%", limits=c(0,30))+
    ggtitle("% of seasonal precip by node-year")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  ####
  
  
  # more daily/seasonal precip metrics
  temp<-somTime[,1:18]
  temp<-merge(temp,seasStats[,c(1,4,6)],by="year") # uses terciles
  temp<-temp %>%
    group_by(codes,year) %>%
    summarise(countNodes=n(), anomName=first(anomName))
  temp <- temp %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  temp$row<-as.numeric(temp$row); temp$col<-as.numeric(temp$col); 
  # alt version from somTime
  # sumYr2<-somTime %>% 
  #   group_by(year) %>%
  #   summarise(sumPrecip=sum(meanPrecip))
  # sumYr2$ltAvg<-mean(sumYr2$sumPrecip)
  # sumYr2$seasAnom<-sumYr2$sumPrecip-sumYr2$ltAvg
  # sumYr2<- merge(sumYr, sumYr2, by='year')
  # sumYr2$anomName<-"normal"
  #   sumYr2$anomName[sumYr2$seasAnom<(-10)] <- "dry"
  #   sumYr2$anomName[sumYr2$seasAnom>(10)] <- "wet"
  # sumYr2 <- sumYr2 %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  #   sumYr2$row<-as.numeric(sumYr2$row); sumYr2$col<-as.numeric(sumYr2$col); 
  # plot  
  ggplot(temp, aes(x=col, y=-row))+
      geom_tile(aes(fill=countNodes))+
      #geom_text(aes(label = countNodes), size=2)+
      facet_grid(as.factor(anomName)~year)+
      #facet_wrap(as.factor(anomName)~year)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", 
                           high = "red", midpoint = 25, 
                           na.value = NA, name="Count")+
      ggtitle("Node counts/year - Wet/Normal/Dry")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
  # counts/anom - swap out sumYr2
    anomCount<-temp %>% group_by(codes, anomName) %>% summarise(counts=sum(countNodes))
    anomCount <- anomCount %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    anomCount$row<-as.numeric(anomCount$row); anomCount$col<-as.numeric(anomCount$col); 
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=counts))+
      geom_text(aes(label = counts), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 400,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Count of days in each node by Seas Precip Tercile")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
  # anom of expected counts by anomaly category
    # SOM node counts by anom category
    # temp2<-sumYr2 %>% group_by(anomName) %>% summarise(n=sum(countNodes))
    #   temp2$prop<-temp2$n/sum(temp2$n)
    # temp3<-sumYr2 %>% group_by(codes) %>% summarise(n=sum(countNodes))
    # adapted to temp df
    temp2<-temp %>% group_by(anomName) %>% summarise(n=sum(countNodes))
      temp2$prop<-temp2$n/sum(temp2$n)
    temp3<-temp %>% group_by(codes) %>% summarise(n=sum(countNodes))
   
    anomCount<-temp %>% group_by(codes, anomName) %>% summarise(counts=sum(countNodes))
    anomCount <- anomCount %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    anomCount$row<-as.numeric(anomCount$row); anomCount$col<-as.numeric(anomCount$col); 
    
    anomCount<-merge(anomCount,temp2, by="anomName")
    anomCount$anom<-(anomCount$counts/anomCount$n)*100
    anomCount<-merge(anomCount, temp3, by="codes")
    anomCount$expCt<-anomCount$n.y*anomCount$prop
    anomCount$anomCt<-anomCount$counts-anomCount$expCt
    # proportion test prop.test()
    propPval<- lapply(seq_along(anomCount$counts),
                      function(i) prop.test(anomCount$counts[i],anomCount$n.y[i],p=anomCount$prop[i],alternative = "two.sided")$p.value)
    binomPval<- lapply(seq_along(anomCount$counts),
                      function(i) binom.test(anomCount$counts[i],anomCount$n.y[i],p=anomCount$prop[i],alternative = "two.sided")$p.value)
    anomCount$propPval<-unlist(propPval)
    anomCount$binomPval<-unlist(binomPval)
    anomCount$sig<-anomCount$propPval<=0.05
    anomCount$sigBinom<-anomCount$binomPval<=0.05
    anomCount$anomCtLabel<-ifelse(anomCount$sig==TRUE,
                                   paste0(round(anomCount$anomCt,1),"*"),
                                   paste0(round(anomCount$anomCt,1))) 
    # add in odds ratio
    anomCount$odds<-anomCount$counts/anomCount$expCt
    anomCount$anomOddsLabel<-ifelse(anomCount$sig==TRUE,
                                  paste0(round(anomCount$odds,2),"*"),
                                  paste0(round(anomCount$odds,2))) 
    
    # remove 'normal'
    anomCount<-subset(anomCount, anomName %in% c("wet","dry"))
    
    # ggplot  
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=anomCt))+
      geom_text(aes(label = anomCtLabel), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "lightblue",mid="white", midpoint = 0,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Anom of expected days in each node by Precip Anom Tercile (*pval<0.05)")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
    # ggplot - odds ratio
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=odds))+
      geom_text(aes(label = anomOddsLabel), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "lightblue",mid="white", midpoint = 1,
                           high = "red", na.value = NA, name="Odds ratio")+
      ggtitle("Odds ratio of expected days in each node by Precip Anom Tercile (*pval<0.05)")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
    
    # do proportion tests/odds ratio for 
    temp<-somTime[,1:18]
    temp<-merge(temp,seasStats[,c(1,4,6)],by="year") # uses terciles
    temp<-temp %>%
      group_by(activityCat,year) %>%
      summarise(activity=n(), anomName=first(anomName))
      # group stats
      temp2<-temp %>% group_by(anomName) %>% summarise(n=sum(activity))
      temp2$prop<-temp2$n/sum(temp2$n)
      temp3<-temp %>% group_by(activityCat) %>% summarise(n=sum(activity))
      # get anoms
      anomCount<-temp %>% group_by(activityCat, anomName) %>% summarise(counts=sum(activity))
      # merge
      anomCount<-merge(anomCount,temp2, by="anomName")
        anomCount$anom<-(anomCount$counts/anomCount$n)*100
      anomCount<-merge(anomCount, temp3, by="activityCat")
        anomCount$expCt<-anomCount$n.y*anomCount$prop
        anomCount$anomCt<-anomCount$counts-anomCount$expCt
        # proportion test prop.test()
        propPval<- lapply(seq_along(anomCount$counts),
                          function(i) prop.test(anomCount$counts[i],anomCount$n.y[i],p=anomCount$prop[i],alternative = "two.sided")$p.value)
        binomPval<- lapply(seq_along(anomCount$counts),
                           function(i) binom.test(anomCount$counts[i],anomCount$n.y[i],p=anomCount$prop[i],alternative = "two.sided")$p.value)
        anomCount$propPval<-unlist(propPval)
        anomCount$binomPval<-unlist(binomPval)
        anomCount$sig<-anomCount$propPval<=0.05
        anomCount$sigBinom<-anomCount$binomPval<=0.05
        anomCount$anomCtLabel<-ifelse(anomCount$sig==TRUE,
                                      paste0(round(anomCount$anomCt,1),"*"),
                                      paste0(round(anomCount$anomCt,1))) 
        # add in odds ratio
        anomCount$odds<-anomCount$counts/anomCount$expCt
        anomCount$anomOddsLabel<-ifelse(anomCount$sig==TRUE,
                                        paste0(round(anomCount$odds,2),"*"),
                                        paste0(round(anomCount$odds,2)))   
        # ggplot of activity cat vs anoms
        ggplot(anomCount, aes(x=activityCat, y=1))+
          geom_tile(aes(fill=odds))+
          geom_text(aes(label = anomOddsLabel), size=4)+
          facet_wrap(~anomName)+
          scale_fill_gradient2(low = "lightblue",mid="white", midpoint = 1,
                               high = "red", na.value = NA, name="count")+
          ggtitle("Odds ratio of expected days in each node by Precip Anom Tercile (*pval<0.05)")+
          theme_bw()+
          theme(axis.text.y=element_blank())
        
        
    
    # counts by seasonal anomaly
    temp<-somTime[,c("year","codes")]
    temp<-merge(temp, seasAvgPrecip, by="year")
    anomCount<-temp %>% group_by(codes, anomName) %>% count(codes)
    anomCount <- anomCount %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    anomCount$row<-as.numeric(anomCount$row); anomCount$col<-as.numeric(anomCount$col);
    # avgs by tercile
    anomCount$avgdays<-ifelse(anomCount$anomName=="wet", round(anomCount$n/14,1),round(anomCount$n/13,1))
    # ggplot  
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=avgdays))+
      geom_text(aes(label = avgdays), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 25,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Average days for each seas anom tercile")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
    # get difference from average for each tercile
    codeCount<-temp %>% group_by(codes) %>% count(codes)
      codeCount$avg<-codeCount$n/40
    anomCount<-merge(anomCount,codeCount,by="codes")
    anomCount$ratio<-round(anomCount$avgdays/anomCount$avg,2)
    # ggplot  
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=ratio))+
      geom_text(aes(label = ratio), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "#67a9cf",mid="#f7f7f7", midpoint = 1,
                           high = "#ef8a62", na.value = NA, name="count")+
      ggtitle("Ratio of anom avg/avg for each seas anom tercile")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
    
    
    # diff between wet and dry year average frequencies
   anomDiffs<-dcast(anomCount,codes~anomName)
   anomDiffs$diffs<-round(anomDiffs$wet-anomDiffs$dry,1)
   anomDiffs <- anomDiffs %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
   anomDiffs$row<-as.numeric(anomDiffs$row); anomDiffs$col<-as.numeric(anomDiffs$col); 
    
   ggplot(anomDiffs, aes(x=col, y=-row))+
     geom_tile(aes(fill=diffs))+
     geom_text(aes(label = diffs), size=4)+
     scale_fill_gradient2(low = "brown",mid="grey", midpoint = 0,
                          high = "green", na.value = NA, name="count")+
     ggtitle("Difference in avg days Wet-Dry Anoms")+
     theme_bw()+
     theme(axis.text.x=element_blank(),
           axis.text.y=element_blank())
    # counts of codes by year
   temp<-somTime[,c("year","codes")]
   yrCount<-temp %>% group_by(codes,year) %>% count(codes)
   yrCount<-dcast(yrCount,codes~year)
   temp2<-yrCount[,c("codes","1996","2012")]
   temp2<-yrCount$`1996`-yrCount$`2012` 
   temp2$codes<-temp$codes
   
   
  # trends in nodes plot; plot with precip 
  
  
  # plot map units
  ggplot(somTime, aes(doy, year)) + 
    geom_tile(aes(fill = mapUnit), colour = "grey") + 
    geom_text(aes(label = codes), size=1.50)+
    scale_fill_gradient2(low = "lightblue", mid = "green",
                         high = "orange", midpoint = (ncols*nrows)/2, space = "Lab",
                         na.value = "grey50", guide = "colourbar")+
    ggtitle("Daily Precip Type Classification")
  # plot error
  ggplot(somTime, aes(doy, year)) + 
    geom_tile(aes(fill = errorDist), colour = "grey") + 
    scale_fill_gradient2(low = "white", mid = "yellow",
                         high = "red", space = "Lab",
                         na.value = "grey50", guide = "colourbar")+
    ggtitle("Daily Precip Type Classification Error")
  # plot daily spearman
  ggplot(somTime, aes(doy, year)) + 
    geom_tile(aes(fill = spearman), colour = "grey") + 
    scale_fill_gradient2(low = "lightblue", mid = "yellow",
                         high = "red", space = "Lab",
                         na.value = "grey50", guide = "colourbar")+
    ggtitle("Daily Map Correlation (Rho)")
  
  #####
  # CREATE COMPOSITES on nodes ../ClimPlot/NARR/ downloadNARR.R and processNARR.R
  NARR<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2020.grd")
  #NARR<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2019.grd") 
  #NARR<-stack("/scratch/crimmins/NARR/processed/anom/GH500_daily_NARR_WUS_1979_2019_pentMean_anomaly.grd")
  #NARR<-stack("/scratch/crimmins/NARR/processed/anom/PWAT_daily_NARR_WUS_1979_2019_pentMean_anomaly.grd")
  
  # dates - find and remove leap days
  startYr<-1979
  compDates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2020-12-31"),1))
    colnames(compDates)<-"date"
  compLayers<-NARR[[match(somTime$date,compDates$date)]]
  compMean<-stackApply(compLayers, somTime$mapUnit, fun=mean)
    compMean<-subset(compMean, order(as.numeric(sapply(strsplit(names(compMean), "_"), tail, 1))))
    names(compMean)<-codeList   
  compSD<-stackApply(compLayers, somTime$mapUnit, fun=sd)
    compSD<-subset(compSD, order(as.numeric(sapply(strsplit(names(compSD), "_"), tail, 1))))
    names(compSD)<-codeList     
 
     # t-test anom
   # library(MKinfer)
    pvalAnom<-stack()
    for (i in 1:(nrows*ncols)){
      tempStack<-compLayers[[which(somTime$mapUnit==i)]]
      # ttest pval
      fun=function(x) { if (is.na(x[1])){ NA } else { t.test(x)$p.value}}
      pval <- calc(tempStack, fun)
      # boot ttest
      #fun=function(x) { if (is.na(x[1])){ NA } else { boot.t.test(x, R=1000)$p.value}}
      #pval <- calc(tempStack, fun)
      names(pval)<-paste0("mapUnit_",i)
      pvalAnom<-stack(pvalAnom,pval)
    }
    names(pvalAnom)<-codeList
    
    pvalAnom <- reclassify(pvalAnom, c(-Inf,0.05,1, 0.05,Inf,NA))

  
  ####  
    
  # plot 
  #   at<-seq(0,60,5)
  #   mapTheme<-rasterTheme(region=brewer.pal(8,"BrBG")) 
  # pComp<-levelplot(compMean, contour=FALSE, margin=FALSE,layout=c(ncols,nrows), 
  #                main="PWAT Composite Patterns JAS 2x3 SOM - NARR 1981-2019", par.settings=mapTheme)+
  #     #contourplot(compSD,linetype = "dashed")+
  #     layer(sp.polygons(az, col = 'gray40', lwd=1))
  # 
  #   library(gridExtra)
  #   grid.arrange(pPrecip, pComp, ncol=1)
    
  # ggplot to plot rasters of diff resolutions
  # First, to a SpatialPointsDataFrame
  #compMean_pts <- rasterToPoints(compMean, spatial = TRUE)
  compMean_pts <- rasterToPoints(compMean, spatial = TRUE)
  cbRaster_pts <- rasterToPoints(cbRaster, spatial = TRUE)
  pvalAnom_pts <- rasterToPoints(pvalAnom, spatial = TRUE)
  # Then to a 'conventional' dataframe
  compMean_df  <- data.frame(compMean_pts)
    compMean_df<- melt(compMean_df, id.vars=c("x", "y"))
    compMean_df<-subset(compMean_df, variable!="optional")
  cbRaster_df  <- data.frame(cbRaster_pts)
    cbRaster_df<- melt(cbRaster_df, id.vars=c("x", "y"))
    cbRaster_df<-subset(cbRaster_df, variable!="optional")
  pvalAnom_df  <- data.frame(pvalAnom_pts)
    pvalAnom_df<- melt(pvalAnom_df, id.vars=c("x", "y"))
    pvalAnom_df<-subset(pvalAnom_df, variable!="optional")  
    rm(compMean_pts)
    rm(cbRaster_pts)
    rm(pvalAnom_pts)
    
  # fortify polygons
    #usPoly <- fortify(us, region="NAME_0")
    library(maps)
    usPoly <- map_data("state")
    mxPoly <- map_data(database = "world", regions = "Mexico")
  
  # GH500 Anoms   
  ggplot() +
        geom_raster(data = cbRaster_df , aes(x = x, y = y, fill = value))+
        geom_point(data = pvalAnom_df, aes(x = x, y = y, color = value))+
        geom_polygon(data=usPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
        geom_polygon(data=mxPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
        stat_contour(data = compMean_df , aes(x = x, y = y, z = value, color = ..level..))+
        scale_fill_gradientn(colors = c("lightblue", "blue", "green","yellow","red"))+
        #scale_color_continuous(type = "viridis") +
        scale_color_gradient2(midpoint = 0, low = "blue", mid = "grey", high = "red", name = "GH" )+
        #scale_color_gradient2(midpoint = 25, low = "brown", mid = "lightgrey", high = "green", name = "PWAT(mm)" )+
        coord_cartesian(xlim = c(-125, -100), ylim = c(25, 50))+
        facet_wrap(~variable)
  
  # PWAT Anoms   
  p<-ggplot() +
    geom_raster(data = cbRaster_df , aes(x = x, y = y, fill = value))+
    geom_point(data = pvalAnom_df, aes(x = x, y = y, color = value))+
    geom_polygon(data=usPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
    geom_polygon(data=mxPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
    stat_contour(data = compMean_df , aes(x = x, y = y, z = value, color = ..level..))+
    scale_fill_gradientn(colors = c("lightblue", "blue", "green","yellow","red"))+
    #scale_color_continuous(type = "viridis") +
    scale_color_gradient2(midpoint = 0, low = "brown", mid = "grey", high = "green", name = "PWAT Anom" )+
    #scale_color_gradient2(midpoint = 0, low = "blue", mid = "grey", high = "red", name = "GH" )+
    #scale_color_gradient2(midpoint = 25, low = "brown", mid = "lightgrey", high = "green", name = "PWAT(mm)" )+
    coord_cartesian(xlim = c(-125, -100), ylim = c(25, 50))+
    facet_wrap(~variable)
  
  # GH500   
  ggplot() +
    geom_raster(data = cbRaster_df , aes(x = x, y = y, fill = value))+
    #geom_point(data = pvalAnom_df, aes(x = x, y = y, color = value))+
    geom_polygon(data=usPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
    geom_polygon(data=mxPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
    stat_contour(data = compMean_df , aes(x = x, y = y, z = value, color = ..level..))+
    scale_fill_gradientn(colors = c("lightblue", "blue", "green","yellow","red"))+
    #scale_color_continuous(type = "viridis") +
    scale_color_gradient2(midpoint = 5800, low = "blue", mid = "yellow", high = "red", name = "GH" )+
    #scale_color_gradient2(midpoint = 5800, low = "brown", mid = "lightgrey", high = "green", name = "PWAT(mm)" )+
    coord_cartesian(xlim = c(-125, -100), ylim = c(25, 50))+
    facet_wrap(~variable)
  

#####
# RASTERVIS COMPOSITE MAPS - 00Z grids from ~/ClimPlot/NARR/downloadProcess00ZNARR.R
  # GH500<-stack("/scratch/crimmins/NARR/processed/00z/GH500_00z_NARR_WUS_1979_2019.grd")
  # PWAT<-stack("/scratch/crimmins/NARR/processed/00z/PWAT_00z_NARR_WUS_1979_2019.grd")
  # PRCP<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2019.grd")
  
  # daily means/totals
  GH500<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2020.grd")
  PWAT<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2020.grd")
  PRCP<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2020.grd")
  
  # create mean composites
  # dates - find and remove leap days
  startYr<-1979
    compDates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2020-12-31"),1))
    colnames(compDates)<-"date"
    compLayers<-GH500[[match(somTime$date,compDates$date)]]
      compGH500<-stackApply(compLayers, somTime$mapUnit, fun=mean)
      compGH500<-subset(compGH500, order(as.numeric(sapply(strsplit(names(compGH500), "_"), tail, 1))))
      names(compGH500)<-codeList  
    compLayers<-PWAT[[match(somTime$date,compDates$date)]]
      compPWAT<-stackApply(compLayers, somTime$mapUnit, fun=mean)
      compPWAT<-subset(compPWAT, order(as.numeric(sapply(strsplit(names(compPWAT), "_"), tail, 1))))
      names(compPWAT)<-codeList
    compLayers<-PRCP[[match(somTime$date,compDates$date)]]
      compPRCP<-stackApply(compLayers, somTime$mapUnit, fun=median)
      compPRCP<-subset(compPRCP, order(as.numeric(sapply(strsplit(names(compPRCP), "_"), tail, 1))))
      names(compPRCP)<-codeList
      
  # plot
      library(rnaturalearth)
        countries<-ne_countries(type = 'countries', scale = 'small')
      #compPRCP[compPRCP ==0] <- NA
       compPRCP[compPRCP <0.254] <- NA
      at<-c(seq(0,15,0.25))
      #   mapTheme<-rasterTheme(region=brewer.pal(8,"BrBG")) 
       mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
       pComp<-levelplot(compPRCP, contour=FALSE, margin=FALSE,layout=c(ncols,nrows), at=at,
                        xlim=c(-130,-90), ylim=c(20,42.5),
                      main="NARR GH500/Precip Composite Patterns Jun 15th-Sep 30th 4x4 SOM (1981-2020)", par.settings=mapTheme)+
              #layer(sp.polygons(aznm, col = 'gray50', lwd=0.5))+
              #layer(sp.polygons(states, col = 'gray50', lwd=0.5, fill='grey90'))+
              layer(sp.polygons(countries, col = 'gray50', lwd=0.5, fill='grey90'))+
             
         
         levelplot(compPRCP, contour=FALSE, margin=FALSE,layout=c(ncols,nrows), at=at,
                   xlim=c(-130,-90), ylim=c(20,42.5),
                   main="NARR GH500/Precip Composite Patterns Jun 15th-Sep 30th 4x4 SOM (1981-2020)", par.settings=mapTheme)+
         
         layer(sp.polygons(states, col = 'gray50', lwd=0.5))+
         
              contourplot(compPWAT,linetype="solid", at=c(31.75), labels=FALSE, col='darkorange',lwd=0.75)+
              #contourplot(compPWAT,linetype="solid", at=c(30), labels=FALSE, col='darkred',lwd=0.75)+
              #contourplot(compGH500,linetype="dotdash", at=seq(5700,5850,25), labels=FALSE, col='gray50',lwd=0.5)+
              #contourplot(compGH500,linetype="dotted", at=c(seq(5880,6000,5)), labels=FALSE, col='gray28',lwd=0.5)
              contourplot(compGH500,linetype="dotdash", at=seq(5700,6000,10), labels=FALSE, col='black',lwd=0.5)

        png("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/SOM4x4_composites.png", width = 11, height = 8.5, units = "in", res = 300L)
        #grid.newpage()
        print(pComp, newpage = FALSE)
        dev.off()   
       
#####  
    
#####
# plot composites on intensity classification
        # RASTERVIS COMPOSITE MAPS
        #GH500<-stack("/scratch/crimmins/NARR/processed/00z/GH500_00z_NARR_WUS_1979_2019.grd")
        #PWAT<-stack("/scratch/crimmins/NARR/processed/00z/PWAT_00z_NARR_WUS_1979_2019.grd")
        #PRCP<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2019.grd")
        
        # daily means/totals
        GH500<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2020.grd")
        PWAT<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2020.grd")
        PRCP<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2020.grd")
        
        # create mean composites
        # dates - find and remove leap days
        startYr<-1979
        compDates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2020-12-31"),1))
        colnames(compDates)<-"date"
        compLayers<-GH500[[match(somTime$date,compDates$date)]]
        compGH500<-stackApply(compLayers, somTime$activityCat, fun=mean)
        #compGH500<-subset(compGH500, order(as.numeric(sapply(strsplit(names(compGH500), "_"), tail, 1))))
            #compGH500<-subset(compGH500, order(c(4,1,2,3)))
        #names(compGH500)<-codeList  
        compLayers<-PWAT[[match(somTime$date,compDates$date)]]
        compPWAT<-stackApply(compLayers, somTime$activityCat, fun=mean)
        #compPWAT<-subset(compPWAT, order(as.numeric(sapply(strsplit(names(compPWAT), "_"), tail, 1))))
            #compPWAT<-subset(compPWAT, order(c(4,1,2,3)))
        #names(compPWAT)<-codeList
        compLayers<-PRCP[[match(somTime$date,compDates$date)]]
        compPRCP<-stackApply(compLayers, somTime$activityCat, fun=median)
        #compPRCP<-subset(compPRCP, order(as.numeric(sapply(strsplit(names(compPRCP), "_"), tail, 1))))
            #compPRCP<-subset(compPRCP, order(c(4,1,2,3)))
        #names(compPRCP)<-codeList
           
        # plot
        library(rnaturalearth)
        countries<-ne_countries(type = 'countries', scale = 'small')
        #compPRCP[compPRCP ==0] <- NA
        compPRCP[compPRCP <0.254] <- NA
        at<-c(seq(0,6,0.25))
        #   mapTheme<-rasterTheme(region=brewer.pal(8,"BrBG")) 
        mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
        pComp<-levelplot(compPRCP, contour=FALSE, margin=FALSE,layout=c(4,1), at=at,
                         main="NARR GH500/Precip Composite Patterns June 15th-Sep 30th Activity Cat (1981-2020)", par.settings=mapTheme)+
          #layer(sp.polygons(aznm, col = 'black', lwd=0.5))+
          layer(sp.polygons(countries, col = 'grey90', lwd=0.5))+
          levelplot(compPRCP, contour=FALSE, margin=FALSE,layout=c(4,1), at=at,
                    main="NARR GH500/Precip Composite Patterns June 15th-Sep 30th Activity Cat (1981-2020)", par.settings=mapTheme)+
          
          layer(sp.polygons(states, col = 'gray50', lwd=0.5))+
          
          contourplot(compPWAT,linetype="solid", at=c(31.75), labels=FALSE, col='darkorange',lwd=0.75)+
          #contourplot(compPWAT,linetype="solid", at=c(30), labels=FALSE, col='darkred',lwd=0.75)+
          #contourplot(compGH500,linetype="dotdash", at=seq(5700,5850,25), labels=FALSE, col='gray50',lwd=0.5)+
          #contourplot(compGH500,linetype="dotted", at=c(seq(5880,6000,5)), labels=FALSE, col='gray28',lwd=0.5)
          contourplot(compGH500,linetype="dotdash", at=seq(5700,6000,12.5), labels=FALSE, col='black',lwd=0.5)
        
        
        
        png("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/SOM4x4_activity_composites.png", width = 11, height = 4, units = "in", res = 300L)
        #grid.newpage()
        print(pComp, newpage = FALSE)
        dev.off()   
#####        
        
        
        
#####
# flash flood data
        load("~/RProjects/SOMs/monsoonPrecip/USGS_FFlood.RData")
        hucs<-merge(hucs, somTime[,c(1,8)], by="date")
        temp<-hucs@data
        ggplot(temp, aes(as.factor(temp$codes), as.numeric(as.character(temp$Peak.Q..cm))))+
          geom_boxplot(varwidth = TRUE)+
          ggtitle("Distribution of Max USGS Peak Flows by nodes")
        # plot on map
        temp <- temp %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
        temp$Lat <- as.numeric(as.character(temp$Lat)); temp$Lon <- as.numeric(as.character(temp$Lon))
        temp$Peak.Q..cm <- as.numeric(as.character(temp$Peak.Q..cm))
        library(maps)
        usPoly <- map_data("state")
        ggplot()+
        geom_polygon(data=usPoly,aes(long, lat, group = group),
                     fill = NA, col = "black", size = 0.2)+
         geom_point(data=temp, aes(x=Lon,y=Lat, color=Peak.Q..cm))+
          facet_wrap(row~col)+
          coord_cartesian(xlim = c(-116,-105), ylim = c(30, 38))+
          scale_colour_gradient(low="lightblue",high="red", name="cfs")+
          ggtitle("USGS Flash Flood Peak Flows (1986-2015)")
 
 # seasonal total precipitation      
 #####       
 ggplot(seasAvgPrecip, aes(year,avgPrecip, fill=as.factor(seasAvgPrecip$anomName)) )+
          geom_bar(stat = 'identity')+
          ggtitle("Regional Average Total Precip (Jun 15-Sep30, mm)")+
          geom_hline(yintercept=mean(seasAvgPrecip$avgPrecip), color="black")+
          geom_hline(yintercept=median(seasAvgPrecip$avgPrecip), color="red")+
          scale_fill_manual(values = c("saddlebrown", "grey", "forestgreen"), name="tercile")
    # count of nodes vs seas Avg precip (monsoon vs non monsoon days)
      temp<-somTime %>%  group_by(codes,year) %>% count(codes)
      temp<-subset(temp, codes=="1_1")
        plot(seasAvgPrecip$percRank, temp$n, main="1_1 day count vs seas avg precip")
          text(seasAvgPrecip$percRank, temp$n, seasAvgPrecip$year, pos=2)
          abline(lm(temp$n~seasAvgPrecip$percRank))
    # stepwise regression of node counts vs seas avg precip
          temp<-somTime %>%  group_by(codes,year) %>% count(codes)
          temp<-dcast(temp, formula = year~codes, value.var = "n")
          temp$precip<-seasAvgPrecip$avgPrecip
          temp<-temp[,-1]
          library(MASS)
          # Fit the full model 
          full.model <- lm(precip ~., data = temp)
          # Stepwise regression model
          library(leaps)          
          models <- regsubsets(precip ~., data = temp, nvmax = 5,
                               method = "seqrep")
          reg.summary = summary(models)
          reg.summary$rsq
          
 # make dummy date for plotting boxplots
 somTime$dummyDate<-as.Date(paste0("2000-",somTime$month,"-",somTime$day),format="%Y-%m-%d")
 
 p1<-ggplot(somTime, aes(dummyDate, percExtent1, group=dummyDate))+
          geom_boxplot(varwidth = TRUE)+
          ggtitle("Precip Coverage>1 (%) by day")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   scale_x_date(labels = date_format("%m/%d"))+
   xlab("Day of Year")+
   ylab("% Extent")

 ggplot(somTime, aes(dummyDate, percExtent5, group=dummyDate))+
   geom_boxplot(varwidth = TRUE)+
   ggtitle("Precip Coverage>5mm (%) by day")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   scale_x_date(labels = date_format("%m/%d"))+
   xlab("Day of Year")+
   ylab("% Extent")
 
 p2<-ggplot(somTime, aes(dummyDate, percExtent10, group=dummyDate))+
   geom_boxplot(varwidth = TRUE)+
   ggtitle("Precip Coverage>10mm (%) by day")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   scale_x_date(labels = date_format("%m/%d"))+
   xlab("Day of Year")+
   ylab("% Extent")
 
 ggplot(somTime, aes(dummyDate, meanPrecip, group=dummyDate))+
   geom_boxplot(varwidth = TRUE)+
   ggtitle("Mean Precip(mm) by day")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   scale_x_date(labels = date_format("%m/%d"))+
   xlab("Day of Year")+
   ylab("mm")
    
 ggplot(somTime, aes(dummyDate, medPrecip, group=dummyDate))+
   geom_boxplot(varwidth = TRUE)+
   ggtitle("Median Precip(mm) by day")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   scale_x_date(labels = date_format("%m/%d"))+
   xlab("Day of Year")+
   ylab("mm")    
        
 p3<-ggplot(somTime, aes(dummyDate, maxPrecip, group=dummyDate))+
   geom_boxplot(varwidth = TRUE)+
   ggtitle("Max Precip(mm) by day")+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   scale_x_date(labels = date_format("%m/%d"))+
   xlab("Day of Year")+
   ylab("mm")    
 # cowplot multi
 plot_grid(p1,p2,p3, align="v", ncol=1)

 library("PerformanceAnalytics")
 #chart.Correlation(somTime[,c("percExtent","percExtent5","percExtent10","maxPrecip","meanPrecip","medPrecip")], histogram=TRUE, pch=19)
 chart.Correlation(somTime[,c("percExtent1","percExtent10","maxPrecip","meanPrecip","medPrecip")], histogram=TRUE, pch=19)
  # rank of values
   temp<-subLayers
   temp[temp ==0] <- NA 
   #histPrecip<-hist(getValues(temp), main="Histogram of all precip values (mm)")
   qPrecip<-quantile(getValues(temp), probs=c(0.1,0.5,0.75,0.90), na.rm=TRUE)
 # active vs non-active days
   # Ellis et al median of coverage as threshold of active vs non-active days
   length( which( somTime$percExtent1 <= median(somTime$percExtent1)))
   length( which( somTime$codes == "1_1"))
   inactiveDays<-somTime %>%
     group_by(year) %>%
     summarize(bloMedian=sum(percExtent1 <= median(somTime$percExtent1)),
               code1_1=sum(codes == "1_1"))
   ggplot(inactiveDays, aes(bloMedian,code1_1))+
     geom_point()+
     geom_text(aes(label=year),hjust=0, vjust=0)+
     geom_abline(intercept = 0, slope = 1)
   cor(inactiveDays$bloMedian,inactiveDays$code1_1)
   # below median classification error by activity class
   somTime$bloMedian<-ifelse(somTime$percExtent1 <= median(somTime$percExtent1),"Inactive","Active")
   table(somTime$bloMedian,somTime$activityCat)
   table(somTime$bloMedian,somTime$activityCat2)
   ggplot(somTime, aes(x=bloMedian,y=percExtent1, fill=activityCat))+
     geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
     ggtitle("Distribution of Daily Extent Precip (%) (compare with median threshold)")+
     theme(legend.position="bottom")+
     facet_wrap(~codes, scales="free_x")
   # somTime[which(somTime$activityCat=="Inactive" & somTime$bloMedian=="Active"),]
   
# grouped extent and mean/max day of year
   temp<-somTime[,c("dummyDate","percExtent1","percExtent10")]
   temp<-melt(temp, id.vars = c("dummyDate"))
   p1<-ggplot(temp, aes(dummyDate, value, group=interaction(dummyDate,variable), color=variable))+
     geom_boxplot(outlier.shape = NA)+
     ggtitle("Precip Coverage>1 (%) by day")+
     #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
     scale_x_date(labels = date_format("%m/%d"))+
     xlab("Day of Year")+
     ylab("% Extent")
   temp<-somTime[,c("dummyDate","meanPrecip","maxPrecip")]
   temp<-melt(temp, id.vars = c("dummyDate"))
   p2<-ggplot(temp, aes(dummyDate, value, group=interaction(dummyDate,variable), color=variable))+
     geom_boxplot(outlier.shape = NA)+
     ggtitle("Mean/Max Precip(mm) by day")+
     #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
     scale_x_date(labels = date_format("%m/%d"))+
     xlab("Day of Year")+
     ylab("mm")+  
      ylim(0,100)
   
   # another version of a doy chart  -- COMBINED Activity/SOM DOY
   countDOY<-somTime %>% group_by(dummyDate) %>% count(codes)
   # reorder factors if necessary
   countDOY$codes<-as.factor(countDOY$codes)
   countDOY$codes<-factor(countDOY$codes, levels=c("1_1","1_2","2_1","2_2","1_3","3_1","3_2","2_3","3_3","1_4","4_1","4_2","2_4","3_4","4_3","4_4"))
      colourCount = length(unique(countDOY$codes))
      getPalette = colorRampPalette(brewer.pal(9, "Paired"))
      #countDOY<-left_join(countDOY,activity)
   
  countCat<-somTime %>% group_by(dummyDate) %>% count(activityCat) 
    countCat$activityCat<-factor(countCat$activityCat, levels=c("Inactive","Active-low","Active-high","Widespread"))
         
  # set colors manual
  colQual<-c("#deebf7",
    "#9ecae1","#4292c6","#08519c",
    "#e5f5e0","#a1d99b","#41ab5d","#006d2c","#00441b",
    "#fff7bc","#fee391","#fec44f","#fe9929","#cc4c02","#993404","#662506")
  
  library(grid) #unit(c(top, right, bottom, left), units)
  # color of text by activity Cat?
  # https://stackoverflow.com/questions/43478980/ggplot2-different-text-colors-for-each-legend-label
  p<-ggplot(countDOY, aes(fill=codes, y=n, x=dummyDate)) + 
     #geom_bar(position="stack", stat="identity")
     geom_bar(position="fill", stat="identity")+
     #scale_fill_brewer(type = "qual",palette = "Paired",direction = 1,aesthetics = "fill")+
     #scale_fill_manual(values = getPalette(colourCount))+
     # scale_fill_manual(values = colQual)+
     scale_fill_manual(values = colQual[order(levels(countDOY$codes))])+
     ggtitle("Node frequency by day through season")+
     geom_bar(data=countCat,aes(color=activityCat, y=n, x=dummyDate, fill=NA),position="fill", stat="identity")+
      scale_color_manual(values = c(NA, "grey", "orange","red"), guide=FALSE)+
     #scale_fill_manual(values = colQual)
      scale_fill_manual("", 
                    breaks = c(levels(countDOY$codes)),
                    #values = c(colQual[order(levels(countDOY$codes))]),
                    values = c(colQual),
                    labels = c(levels(countDOY$codes)))+
    theme(legend.position="bottom")
    #theme(plot.margin=unit(c(0,0,0.2,0),"npc"))
 p<- p+guides(fill=guide_legend(nrow=1))+
   annotation_custom(rectGrob(x=unit(0.5, "npc"),y=unit(-0.05, "npc"),height = unit(0.20, "npc"), width = unit(0.1, "npc"),  
                               gp=gpar(col="red", fill=NA)))
   
 # Code to override clipping
 gt <- ggplot_gtable(ggplot_build(p))
 gt$layout$clip[gt$layout$name == "panel"] <- "off"
 #grid.draw(gt)
 
 png("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/testDOY.png", width = 11, height = 8.5, units = "in", res = 300L)
 #grid.newpage()
 grid.draw(gt)
 dev.off() 
 
    # activity chart 
    countDOY<-somTime %>% group_by(dummyDate) %>% count(activityCat)
     ggplot(countDOY, aes(fill=as.factor(activityCat), y=n, x=dummyDate)) + 
     #geom_bar(position="stack", stat="identity")
     geom_bar(position="fill", stat="identity")+
     scale_fill_brewer(type = "qual",
                       palette = "Paired",
                       direction = 1,
                       aesthetics = "fill")+
     ggtitle("Activity frequency by day through season")
     
     # single out a code
     temp<-somTime
     temp$oneCode<-ifelse(temp$codes=="4_4", "code4_3","other")
     countDOY<-temp %>% group_by(dummyDate) %>% count(oneCode)
     ggplot(countDOY, aes(fill=as.factor(oneCode), y=n, x=dummyDate)) + 
       #geom_bar(position="stack", stat="identity")
       geom_bar(position="fill", stat="identity")+
       scale_fill_brewer(type = "qual",
                         palette = "Paired",
                         direction = 1,
                         aesthetics = "fill")+
       ggtitle("Activity frequency by day through season")
   
   
   # another version of a doy chart
   countDOY<-somTime %>% group_by(year) %>% count(codes)
   ggplot(countDOY, aes(fill=as.factor(codes), y=n, x=year)) + 
     #geom_bar(position="stack", stat="identity")
     geom_bar(position="fill", stat="identity")+
     scale_fill_manual(values = getPalette(colourCount))+
     ggtitle("Node frequency by year through season")
   # activity by year
   countDOY<-somTime %>% group_by(year) %>% count(activityCat)
   ggplot(countDOY, aes(fill=activityCat, y=n, x=year)) + 
     #geom_bar(position="stack", stat="identity")
     geom_bar(position="fill", stat="identity")+
     scale_fill_manual(values = getPalette(colourCount))+
     ggtitle("Activity frequency by year through season")
   
   # mean first/last day of occurrence 
   somTime$dummyDate<-as.Date(paste0("2000-",somTime$month,"-",somTime$day),format="%Y-%m-%d")
   firstDoy<-somTime %>% group_by(codes,year) %>% slice_min(dummyDate, n=1) # 
   firstDoy<-firstDoy %>% group_by(codes) %>% summarise(firstdoyMean=mean(dummyDate),
                                                        firstdoySD=sd(dummyDate))
   lastDoy<-somTime %>% group_by(codes,year) %>% slice_max(dummyDate, n=1) # 
   lastDoy<-lastDoy %>% group_by(codes) %>% summarise(lastdoyMean=mean(dummyDate),
                                                        lastdoySD=sd(dummyDate))
   doyStats<-cbind.data.frame(firstDoy,lastDoy[,2:3])
   doyStats<-cbind.data.frame(
      melt(doyStats,id.vars="codes",measure.vars = c("firstdoyMean","lastdoyMean")),
      melt(doyStats,id.vars="codes",measure.vars = c("firstdoySD","lastdoySD")))
   doyStats<-doyStats[,c(1,2,3,6)]
    colnames(doyStats)[3:4]<-c("mean","sd")
   
   ggplot(doyStats,aes(codes,mean, color=variable))+
     geom_point()+
     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2)+
     facet_wrap(~codes,scales="free_x")+
     ggtitle("Mean/sd first and last doy for each node")
   
   ggplot(doyStats, aes(mean,as.factor(codes),color=variable))+
     geom_point()+
     geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd), width=.2)+
     ggtitle("Mean First/Last day of occurrence by node")
   
   # make dummy date for plotting boxplots
   somTime$dummyDate<-as.Date(paste0("2000-",somTime$month,"-",somTime$day),format="%Y-%m-%d")
   ggplot(somTime, aes(x=codes, y=dummyDate))+
     geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
     ggtitle("Distribution of Date of occurrence for each node")+
     theme(legend.position="bottom")+
     facet_wrap(~codes, scales="free_x")
   # nodes as factors
  
   
   
   
 #####       
   
#####
   # create intensity index
   # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
   library(factoextra)
   res.pca <- prcomp(somTime[,c(13:17)], scale = TRUE)
   fviz_eig(res.pca)
   fviz_pca_var(res.pca,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE     # Avoid text overlapping
   )
   somTime$pc1<-res.pca$x[,1]
   somTime$pc2<-res.pca$x[,2]
   
#####   
   
   
        
#####        
        
        
#####
# station metrics using ACIS
        # bind station datas
        somTime<-cbind.data.frame(somTime,stationDaily)
        
        stationTemp<-melt(somTime[,c("codes","Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso")])
        stationTemp$variable = factor(stationTemp$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
        
        ggplot(stationTemp, aes(as.factor(stationTemp$codes), value, color=variable))+
          geom_boxplot(varwidth = FALSE, outlier.shape = 20)+
          facet_wrap(~codes,scales = "free")+
          ggtitle("Distribution of Daily ACIS Station Precip by nodes")+
          ylim(0,90)
           
    # percent of seas total by node
        stationTemp<-melt(somTime[,c("codes","Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque")])
          #stationTotal<-stationTemp %>%  group_by(codes, variable) %>% summarise(sumPrecip=sum(value))
          stationPerc<- stationTemp %>%
                  group_by(variable) %>%
                  mutate(sumPrecip=sum(value, na.rm = TRUE)) %>%
                  group_by(codes, add=TRUE) %>%
                  summarise(perTotal=sum(value, na.rm = TRUE)/min(sumPrecip))
          stationPerc$variable = factor(stationPerc$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
          
          
          ggplot(stationPerc, aes(variable,perTotal, fill=variable))+
            geom_bar(stat='identity')+
            facet_wrap(~codes)+
            ylab("avg proportion of seas total")+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())+
            ggtitle("Average proportion of seasonal total precip by node (ACIS Stations)")
          
# station metrics using PRISM extracts
          # bind station datas
          somTime<-cbind.data.frame(somTime,stationDailyPRISM)
          
          stationTemp<-melt(somTime[,c("codes","Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso")])
          stationTemp$variable = factor(stationTemp$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
          
          ggplot(stationTemp, aes(as.factor(stationTemp$codes), value, color=variable))+
            geom_boxplot(varwidth = FALSE, outlier.shape = 20)+
            facet_wrap(~codes,scales = "free")+
            ggtitle("Distribution of Daily PRISM Station Precip by nodes")+
            ylim(0,90)
         
          stationTemp<-melt(somTime[,c("activityCat","Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso")])
          stationTemp$variable = factor(stationTemp$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
           
          ggplot(stationTemp, aes(as.factor(stationTemp$activityCat), value, color=variable))+
            geom_boxplot(varwidth = FALSE, outlier.shape = NA)+
            facet_wrap(~activityCat, nrow=1, scales = "free")+
            ggtitle("Distribution of Daily PRISM Station Precip by Activity Category")+
            ylim(0,20)
          
          # percent of seas total by node
          stationTemp<-melt(somTime[,c("codes","Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque")])
          #stationTotal<-stationTemp %>%  group_by(codes, variable) %>% summarise(sumPrecip=sum(value))
          stationPerc<- stationTemp %>%
            group_by(variable) %>%
            mutate(sumPrecip=sum(value, na.rm = TRUE)) %>%
            group_by(codes, add=TRUE) %>%
            summarise(perTotal=sum(value, na.rm = TRUE)/min(sumPrecip))
          stationPerc$variable = factor(stationPerc$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
          
          
          ggplot(stationPerc, aes(variable,perTotal, fill=variable))+
            geom_bar(stat='identity')+
            facet_wrap(~codes)+
            ylab("avg proportion of seas total")+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())+
            ggtitle("Average proportion of seasonal total precip by node (PRISM)")
          
          
          # pie chart
          ggplot(stationPerc, aes(x="",y=perTotal,fill=codes))+
            geom_bar(stat='identity', width = 1)+
            coord_polar("y", start=0)+
            facet_wrap(~variable, nrow=2)+
            theme_void()
          
          
          # percent of seas total by node
          stationTemp<-melt(somTime[,c("activityCat","Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque")])
          #stationTotal<-stationTemp %>%  group_by(codes, variable) %>% summarise(sumPrecip=sum(value))
          stationPerc<- stationTemp %>%
            group_by(variable) %>%
            mutate(sumPrecip=sum(value, na.rm = TRUE)) %>%
            group_by(activityCat, add=TRUE) %>%
            summarise(perTotal=sum(value, na.rm = TRUE)/min(sumPrecip))
          stationPerc$variable = factor(stationPerc$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
          
          ggplot(stationPerc, aes(variable,perTotal, fill=variable))+
            geom_bar(stat='identity')+
            facet_wrap(~activityCat, nrow=1, scales = "free")+
            ylab("avg proportion of seas total")+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())+
            ggtitle("Average proportion of seasonal total precip by node (PRISM)")+
            ylim(0,1)           
          
          # pie chart
          ggplot(stationPerc, aes(x="",y=perTotal,fill=activityCat))+
            geom_bar(stat='identity', width = 1)+
            coord_polar("y", start=0)+
            facet_wrap(~variable, nrow=2)+
            theme_void()
  
          # precip rate 
          stationPerc<- stationTemp %>%
            group_by(variable) %>%
            mutate(sumPrecip=sum(value, na.rm = TRUE),
                   n=n()) %>%
            group_by(activityCat, add=TRUE) %>%
            summarise(perTotal=sum(value, na.rm = TRUE)/n())
          stationPerc$variable = factor(stationPerc$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
          
          ggplot(stationPerc, aes(variable,perTotal, fill=variable))+
            geom_bar(stat='identity')+
            facet_wrap(~activityCat, nrow=1)+
            ylab("avg precip by class (mm)")+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())+
            ggtitle("Average precipitation by activity class (PRISM)")
          
          
          
#####        
  ### sample day maps showing progression of extent and max precip
          # plot 10 random maps from selected node to assess quality
          at<-c(seq(0.01,50,1),105)
          mapTheme <- rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
          temp<-subLayers[[which(somTime$date>=as.Date("1986-09-22") & somTime$date<=as.Date("1986-09-30"))]]
          temp[temp == 0] <- NA  
          pExamp<-levelplot(temp, contour=FALSE, 
                           margin=FALSE, par.settings=mapTheme, at=at,
                           main="Sample days - PRISM-daily 1981-2019")+
            layer(sp.polygons(aznm, col = 'gray40', lwd=1))
          somTime[which(somTime$date>=as.Date("1983-09-25") & somTime$date<=as.Date("1983-09-30")),c("date","codes","percExtent","percExtent10","maxPrecip","activityCat")]
        
#####        
        
  
##### Predict classification of new events - use getDailyPRISMppt.R
  newPrcp<-stack("~/RProjects/SOMs/monsoonPrecip/SWUS_061522_071722_PRISM_daily_prcp.grd")
  # crop to region
  e <- extent(-115.5,-106,31.3, 37.5) 
  newPrcp <-crop(newPrcp,e)
  # apply mask
  newPrcp <- mask(newPrcp, mask)
  # convert layers to dataframe
    new.layers.df<-(as.data.frame(newPrcp, long=TRUE, xy=TRUE))
    colnames(new.layers.df)<-c("lon","lat","date","value")  
  # long to wide
    new.df.wide<-dcast(new.layers.df, formula = date~lat+lon, value.var = "value")
    new.df.wide[is.na(new.df.wide)] <- 0
  # make prediction based on trained SOM
    som.prediction <- predict(som.gh500, newdata = as.matrix(new.df.wide[,2:ncol(new.df.wide)]),
                              trainX =as.matrix(df.wide[idx,2:ncol(df.wide)]))
    
  # look at dates and mapped units
   predictedUnits<-cbind.data.frame(names(newPrcp),som.prediction$unit.classif)
   colnames(predictedUnits)<-c("date","mapUnit")
   predictedUnits<-merge(predictedUnits,merge(code_grid,activity, by="codes"), by="mapUnit")

   # make calendar plot
   predictedUnits$date<-as.Date(predictedUnits$date, format="X%Y.%m.%d")      
   predictedUnits$wday<-wday(predictedUnits$date, label = T, week_start = 7)
   predictedUnits$week<-epiweek(predictedUnits$date)
   predictedUnits$month<-format(predictedUnits$date,"%m")
   predictedUnits %>%
     ggplot(aes(wday,-week, fill = activityCat)) +
     geom_tile(colour = "white")  + 
     geom_text(aes(label = codes), size = 3) +
     theme(aspect.ratio = 1/2,
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y = element_blank(),
           panel.grid = element_blank(),
           axis.ticks = element_blank(),
           panel.background = element_blank(),
           strip.background = element_blank(),
           strip.text = element_text(face = "bold", size = 15),
           panel.border = element_rect(colour = "black", fill=NA, size=1)) +
     scale_fill_manual(values = c("lightblue", "palegreen", "yellow","violet"))+
     facet_wrap(~month, nrow = 4, ncol = 1, scales = "free") +
     labs(title = "SOM precip classifications - Monsoon 2022")
   # mean regional precip stats
   predictedUnits <- predictedUnits[order(as.Date(predictedUnits$date, format="%Y-%m-%d")),]
   
   predictedUnits$meanPrecip<-apply(new.df.wide[,2:ncol(new.df.wide)], 1, mean, na.rm=TRUE)# mean regional precip
   
   table(predictedUnits$activityCat)
   table(predictedUnits$activityCat2)  
    table(somTime$activityCat)/43
    
         
   #####
   # try with AHPS data ~ adapted from stateMaps.R
   
   # generate dates -- keep with PRISM date
   allDates<-seq(as.Date("2021-07-01"), as.Date("2021-08-30"),1)
   #####
   # download daily netcdf files from NOAA
   ahpsStack<-stack()
   for(i in 1:length(allDates)){
     #build URL
     URL<-paste0("https://water.weather.gov/precip/downloads/",
                 format(allDates[i],"%Y"),"/",format(allDates[i],"%m"),"/",format(allDates[i],"%d"),
                 "/nws_precip_1day_",format(allDates[i],"%Y%m%d"),"_conus.nc")
     download.file(URL, destfile = paste0("/home/crimmins/RProjects/SOMs/temp/",allDates[i],".nc"), method="curl")
     temp<-raster(paste0("/home/crimmins/RProjects/SOMs/temp/",allDates[i],".nc"), varname="observation")
     ahpsStack <- stack(ahpsStack,temp)
     print(allDates[i])
     print(URL)
   }
   # clean up files in temp dir
   
   # use cbRaster for climo/reference grid
   temp<-cbRaster[[1]]
   # reproject AHSP to lat/lon
   proj4string(temp) <- "+proj=longlat +datum=WGS84 +no_defs"
   ahpsStack <- projectRaster(ahpsStack, temp)
   
   # crop to region
   newPrcp<-ahpsStack*25.4
   e <- extent(-115.5,-106,31.3, 37.5) 
   newPrcp <-crop(newPrcp,e)
   # apply mask
   newPrcp <- mask(newPrcp, mask)
   # convert layers to dataframe
   new.layers.df<-(as.data.frame(newPrcp, long=TRUE, xy=TRUE))
   colnames(new.layers.df)<-c("lon","lat","date","value")  
   # long to wide
   new.df.wide<-dcast(new.layers.df, formula = date~lat+lon, value.var = "value")
   new.df.wide[is.na(new.df.wide)] <- 0
   # make prediction based on trained SOM
   som.prediction <- predict(som.gh500, newdata = as.matrix(new.df.wide[,2:ncol(new.df.wide)]),
                             trainX =as.matrix(df.wide[idx,2:ncol(df.wide)]))
   
   # look at dates and mapped units
   predictedUnits<-cbind.data.frame(allDates,som.prediction$unit.classif)
   colnames(predictedUnits)<-c("date","mapUnit")
   predictedUnits<-merge(predictedUnits,merge(code_grid,activity, by="codes"), by="mapUnit")
   
   # make calendar plot
   #predictedUnits$date<-as.Date(predictedUnits$date, format="X%Y.%m.%d")      
   predictedUnits$wday<-wday(predictedUnits$date, label = T, week_start = 7)
   predictedUnits$week<-epiweek(predictedUnits$date)
   predictedUnits$month<-format(predictedUnits$date,"%m")
   predictedUnits %>%
     ggplot(aes(wday,-week, fill = activityCat)) +
     geom_tile(colour = "white")  + 
     geom_text(aes(label = codes), size = 3) +
     theme(aspect.ratio = 1/2,
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y = element_blank(),
           panel.grid = element_blank(),
           axis.ticks = element_blank(),
           panel.background = element_blank(),
           strip.background = element_blank(),
           strip.text = element_text(face = "bold", size = 15),
           panel.border = element_rect(colour = "black", fill=NA, size=1)) +
     scale_fill_manual(values = c("lightblue", "palegreen", "yellow","violet"))+
     facet_wrap(~month, nrow = 3, ncol = 1, scales = "free") +
     labs(title = "SOM precip classifications (AHPS) - Monsoon 2021")
   #####
   