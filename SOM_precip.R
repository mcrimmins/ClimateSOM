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

ptm <- proc.time()

# set rasteroptions
rasterOptions(progress = 'text')

# functions
leap_every_year <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}
perc.rank <- function(x) trunc(rank(x))/length(x)

# map layers
states <- getData('GADM', country='United States', level=1)
  az<-subset(states, NAME_1=="Arizona")
  nm<-subset(states, NAME_1=="New Mexico")
  #aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico")
  aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico"| NAME_1=="California" | NAME_1=="Nevada" | NAME_1=="Utah" | NAME_1=="Colorado")
us <- getData('GADM', country='United States', level=0)
mx <- getData('GADM', country='Mexico', level=0)
#cn <- getData('GADM', country='Canada', level=0)
# HUC4
huc4<-rgdal::readOGR(dsn="~/RProjects/SOMs/monsoonPrecip/shapes", layer="huc4clip")
# station points
stations<-cbind.data.frame(c("Tucson","Phoenix","Flagstaff","Las Vegas","El Paso","Albuquerque"),
                           c(32.1145,33.4373,35.1404,36.0840,31.8053,35.0433),
                           c(-110.9392,-112.0078,-111.6690,-115.1537,-106.3824,-106.6129))
colnames(stations)<-c("stations","latitude","longitude")
coordinates(stations)= ~ longitude+latitude


# load PRISM for SW Region - from dailyDownloadPRISM.R
 prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2019_PRISM_daily_prcp.grd") 
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
  #subLayers<-prcp # FOR PERC GRIDS
  subLayers<-prcp[[which(dates$month %in% mos)]]
    # crop to region
    subLayers<-crop(subLayers,e)
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
    # names
    seasAvgPrecip$anomName<-"normal"
    seasAvgPrecip$anomName[seasAvgPrecip$percRank<=0.33] <- "dry"
    seasAvgPrecip$anomName[seasAvgPrecip$percRank>=0.66] <- "wet"
    # get station extracts
    stationSumSeas<-as.data.frame(t(raster::extract(sumSeas, stations)))
      colnames(stationSumSeas)<-stations$stations
    # get daily station data  
    stationDaily<-as.data.frame(t(raster::extract(subLayers, stations)))
      colnames(stationDaily)<-stations$stations
    # get watershed stats
    huc4SumSeas<-as.data.frame(t(raster::extract(sumSeas, huc4, fun=mean)))
      colnames(huc4SumSeas)<-huc4$NAME
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
      
  # convert layers to dataframe
  #layers.df<-(as.data.frame(subLayers, long=TRUE, xy=TRUE))
  #colnames(layers.df)<-c("lon","lat","date","value")  
  # long to wide
  #df.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
  
  # save raw precip
  #save(df.wide, file="~/RProjects/SOMs/monsoonPrecip/AZNM_PRISM_JAS.RData")
  load("~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS.RData")
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

  
#  SOM screening #####
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
# save(diagnostics, file = "~/RProjects/SOMs/monsoonPrecip/diagnostics_CP2.RData")
######

#####  
# kohonen SOM
  # nrows=3
  # ncols=4
  # ptm <- proc.time()
  #   #som.gh500 <- som(as.matrix(df.wide[,2:ncol(df.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
  #   #set.seed(999) #keep set seed 16, 11 upper left/wet, 9,6 upper right dry/UL wet, 8 UL Wet/LR dry, 7 LR Wet/LL Dry, 5 LR wet/LL dry, 4/100 UL wet/UR dry, 101 wet UR/dry LL, 102/104 dry UL/wet UR, 103 UL Wet/LR Dry
  #   set.seed(996) # 3x4
  #   som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
  #                         #grid = somgrid(ncols, nrows, "rectangular"),
  #                         grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("bubble")),
  #                         #alpha = c(0.05, 0.001), # for online
  #                         radius = c(1,0), # c(3,0.33)
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
  #####
  
  ##### 
  # CP2 - 2 phase SOM training
  # nrows=3
  # ncols=4
  # ptm <- proc.time()
  # #som.gh500 <- som(as.matrix(df.wide[,2:ncol(df.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
  # set.seed(999) #keep set seed 16, 11 upper left/wet, 9,6 upper right dry/UL wet, 8 UL Wet/LR dry, 7 LR Wet/LL Dry, 5 LR wet/LL dry, 4/100 UL wet/UR dry, 101 wet UR/dry LL, 102/104 dry UL/wet UR, 103 UL Wet/LR Dry
  # #set.seed(996) # for 4x4
  # som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
  #                       #grid = somgrid(ncols, nrows, "rectangular"),
  #                       grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
  #                       #alpha = c(0.05, 0.001), # for online
  #                       radius = c(3,1), #(4,1) for 4x5
  #                       mode= "pbatch", # "pbatch" or "online"
  #                       #mode= "online", # "pbatch" or "online"
  #                       #maxNA.fraction = 0.999,
  #                       cores = 7,
  #                       rlen = 500,
  #                       dist.fcts = "sumofsquares")
  # som.gh500.2 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
  #                       #grid = somgrid(ncols, nrows, "rectangular"),
  #                       grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
  #                       #alpha = c(0.05, 0.001), # for online
  #                       radius = c(2,0.33), # c(3,0.33) for 3x5
  #                       mode= "pbatch", # "pbatch" or "online"
  #                       #mode= "online", # "pbatch" or "online"
  #                       #maxNA.fraction = 0.999,
  #                       init= som.gh500$codes,
  #                       cores = 7,
  #                       rlen = 700,
  #                       dist.fcts = "sumofsquares")
  #     ## quantization error:
  #     mean(som.gh500.2$distances)
  #     ## topographical error measures:
  #     source("topo.error.R")
  #     topo.error(som.gh500.2, "nodedist")
  # proc.time() - ptm
  # som.gh500<-som.gh500.2
  #####
  
 # test<-calc(subLayers, sd)
#  test2<-calc(subLayers, mean)
  
  # save SOM output
  #save(som.gh500, file = "~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS_SOM3x4_CP2.RData")
  load("~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS_SOM3x4_CP2.RData") #1- rad(4,1):rlen:5000, #2-rad(3,0.33):rlen:7000
  nrows=3
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
                     round(table(somTime$codes)/length(unique(somTime$year)),1),"dy/yr)")
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
  
  #cbRaster[cbRaster < 0.254] <- NA  
  #cbRaster[cbRaster ==0] <- NA  
  at<-c(seq(0,30,0.3))
  mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  #at<-c(seq(0,10,0.25))
  #mapTheme <- rasterTheme(region=rev(terrain_hcl(12)))
  #at<-c(seq(0,1,0.05))
  #mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
  pPrecip<-levelplot(cbRaster, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,
                     names.attr=col.titles,
            main="Precip Patterns JAS 3x4 SOM - PRISM-daily 1981-2019")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))+
    layer(sp.points(xyCentroid[panel.number()],
                  pch=20, cex=1, col="black"))
  
    #layer(sp.polygons(huc4, col = 'gray20', lwd=0.75))
    #layer(sp.polygons(us, col = 'gray40', lwd=1))+
    #layer(sp.polygons(mx, col = 'gray40', lwd=1))
  # png("/home/crimmins/RProjects/SOMs/precip_SOM_cpc.png", width = 10, height = 6, units = "in", res = 300L)
  # #grid.newpage()
  # print(p, newpage = FALSE)
  # dev.off() 
  
  # spatial pattern correlation
  library(pcaPP)
  somTime$kendall<-NA
  temp<-subLayers
    temp[is.na(temp)] <- 0
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
    somTime$kendall[i]<-cor.fk(values(cbRaster[[somTime$mapUnit[i]]]), values(temp[[i]]))
    
    somTime$rmse[i]<- sqrt(sum((values(subLayers[[i]])-values(cbRaster[[somTime$mapUnit[i]]]))^2, na.rm = TRUE)/length(values(subLayers[[i]])))
    somTime$mae[i]<- sum(abs(values(subLayers[[i]])-values(cbRaster[[somTime$mapUnit[i]]])), na.rm = TRUE)/length(values(subLayers[[i]]))
    print(i)
  }
  rm(temp)
  
  mean(somTime$spearman, na.rm=TRUE)
  median(somTime$spearman, na.rm=TRUE)
  
  mean(somTime$rmse, na.rm=TRUE)
  median(somTime$rmse, na.rm=TRUE)
  
  mean(somTime$pearson, na.rm=TRUE)
  median(somTime$pearson, na.rm=TRUE)
  
  mean(somTime$kendall, na.rm=TRUE)
  median(somTime$kendall, na.rm=TRUE)
  #mean(somTime$pearson, na.rm=TRUE)
  
  ##### 
  # Kendall tau corrs of each day against each map unit
  library(pcaPP)
  temp<-subLayers
  temp[is.na(temp)] <- 0
  
  tauCodeBook <- data.frame(matrix(vector(), 0, length(codeList),
                         dimnames=list(c(), codeList)),
                  stringsAsFactors=F)
  
  for(i in 1:nrow(somTime)){
      for(j in 1:length(codeList)){
        #tauCodeBook[i,j]<-cor.fk(values(cbRaster[[j]]), values(temp[[i]]))
        tauCodeBook[i,j]<-cor(values(cbRaster[[j]]), values(temp[[i]]),
                                 use = "na.or.complete", method="spearman")
        
      }
    print(i)
  }
  somTime$maxTau<-max.col(tauCodeBook)
  somTime$maxTauVal<-apply(tauCodeBook, 1, FUN=max)
  somTime$maxTauCode<-codeList[somTime$maxTau]
  somTime$codeDiff<-somTime$codes==somTime$maxTauCode
     # if tau - NA, replace with 1_1
        somTime$maxTau[is.na(somTime$maxTau)]<-1
  somTime$unitDiff<-somTime$mapUnit-somTime$maxTau      
  
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
                     main="Lowest classification error JAS 3x4 SOM - PRISM-daily 1981-2019")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))

  # max Tau
  highestCorr<-somTime %>% group_by(maxTauCode) %>% slice_max(maxTauVal, n=1)
  corrRaster<-stack()
  for (i in 1:length(codeList)) {
    temp<-subLayers[[which(somTime$date==highestCorr$date[which(highestCorr$maxTau==i)])]]
    corrRaster<-stack(corrRaster, temp)
  }
  corrRaster[corrRaster == 0] <- NA  
  at<-c(seq(0,50,1),100)
  mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  pHighCorr<-levelplot(corrRaster, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,
                       main="Highest tau JAS 3x4 SOM - PRISM-daily 1981-2019")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  
  
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
  somTime$percExtent5<-(rowSums(df.wide[idx,2:ncol(df.wide)]>=5)/(ncol(df.wide)-1))*100 # percent extent >5 mm
  somTime$percExtent10<-(rowSums(df.wide[idx,2:ncol(df.wide)]>=10)/(ncol(df.wide)-1))*100 # percent extent >0
  somTime$maxPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, max) # max value of day
  somTime$meanPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, mean)# mean regional precip
  somTime$medPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, median) # max value of day  
  somTime$percZero<-(rowSums(df.wide[idx,2:ncol(df.wide)]==0)/(ncol(df.wide)-1))*100 # percent 0 precip  
  somTime$sumPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, sum) # max value of day 
  
  #####
  # add in climate indices
  load("~/RProjects/SOMs/monsoonPrecip/climInd.RData")
  climInd<-climInd[,c("date","phase","amplitude","ONI")]
  climInd$date<-climInd$date+1 # align with PRISM doy
    colnames(climInd)[1]<-"date.c"
  somTime<-merge(somTime,climInd, by.x="date",by.y="date.c")
   # ENSO phase
  somTime$ENSO<-"Neutral"
    somTime$ENSO<-ifelse(somTime$ONI <= -0.5, "La Nina", somTime$ENSO)   
    somTime$ENSO<-ifelse(somTime$ONI >= 0.5, "El Nino", somTime$ENSO)   
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
      geom_col() +
      #geom_bar(position="fill", stat="identity") +
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
    
    
  #####    
  # counts of nodes by day of year
    countDOY<-somTime %>% group_by(month,day) %>% count(codes)
      countDOY <- countDOY %>% group_by(month,day) %>% summarize(maxCount=max(n),
                                                                 maxNode= codes[which.max(n)])
    countDOY$date<-as.Date(paste0(countDOY$month,"-",countDOY$day,"-2016"), format="%m-%d-%Y")      
    countDOY$wday<-wday(countDOY$date, label = T, week_start = 7)
    countDOY$week<-epiweek(countDOY$date)
    countDOY %>%
      ggplot(aes(wday,-week, fill = maxCount)) +
      geom_tile(colour = "white")  + 
      geom_text(aes(label = maxNode), size = 3) +
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
      facet_wrap(~month, nrow = 3, ncol = 1, scales = "free") +
      labs(title = "Most frequent nodes by day of year")
    #####  
    
    
  # distribution of daily precip values by node
  extentPerc<-somTime[,c(8,13,14,15)]  
    extentPerc<-melt(extentPerc, id.vars="codes")
    ggplot(extentPerc, aes(x=codes, y=value, fill=variable))+
      geom_boxplot(varwidth = FALSE, position = "dodge2", outlier.alpha = 0.2)+
      ggtitle("Distribution of Daily Extent Precip (%) by nodes")+
      theme(legend.position="bottom")+
      facet_wrap(~codes, scales="free_x")
    
  ggplot(somTime, aes(as.factor(somTime$codes), percExtent))+
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
  
  # plot 10 random maps from selected node to assess quality
  at<-c(seq(0.01,50,1),200)
  mapTheme <- rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
  temp<-subLayers[[idx]]
  temp[temp == 0] <- NA  
  pRand<-levelplot(temp[[sample(which(somTime$codes=="3_1"),10)]], contour=FALSE, at=at,
                     margin=FALSE, par.settings=mapTheme, 
                     main="Precip Patterns from 3_1 - PRISM-daily 1981-2019")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  
  # plot single day
  at<-c(seq(0,50,1),100)
  levelplot(temp[[3420]], contour=FALSE,
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
  at<-c(seq(0,20,0.25),35)
  #at<-c(seq(0,200,5),2000)
  mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), 
                     main="Composite Max Precip JAS 3x4 SOM - PRISM-daily 1981-2019", xlab=NULL, ylab=NULL)+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  # plot median 
  at<-c(seq(0,1,0.05))
  mapTheme<-rasterTheme(region = c("saddlebrown", "sandybrown", "grey","greenyellow","forestgreen"))
  pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at, 
                   main="Composite Median Percentile Precip JAS 3x4 SOM - PRISM-daily 1981-2019", xlab=NULL, ylab=NULL)+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  #####

  #####
  # Composites of seas totals in given years
  yr<-1999
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
    # plot 
    at<-c(seq(0,100,5),250)
    mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
    pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                     main=paste0("Composite Precipitation ",yr), xlab=NULL, ylab=NULL)+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    # percent of seasonal total
    percNode<-(comp/calc(comp, sum, na.rm=TRUE))*100
      names(percNode)<-codeList
    at<-seq(0,100,2.5)
    mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
      pPerc<-levelplot(percNode, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,  
                       main=paste0("Percent of seasonal total by node: ",yr), xlab=NULL, ylab=NULL)+
      layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    # normalized contribution
      # normalize by count?
      counts<-as.vector(table(somTime$codes[which(somTime$year==yr)]))
      norm<-comp/counts
      
      
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
      anomName<-"normal"
      anomYrs<-seasAvgPrecip$`unique(subDates$year)`[which(seasAvgPrecip$anomName==anomName)]
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
      
  ####    
      
      
      
  # Clustering of nodes
  ## Show the U matrix
  Umat <- plot(som.gh500, type="dist.neighbours", main = "SOM neighbour distances")
  ## use hierarchical clustering to cluster the codebook vectors
  som.hc <- cutree(hclust(object.distances(som.gh500, "codes")),12)
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
  # sammon mapping
  library(MASS)
  gh500.codes <- som.gh500$codes
  dis <- dist(as.matrix(som.gh500$codes[[1]]))
  gh500.sam <- sammon(dis)
  plot(gh500.sam$points, type="n")
  text(gh500.sam$points,labels=as.character(1:nrow(code_grid)))
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
    text(gh500.sam$points,labels=as.character(1:nrow(code_grid)))
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
  ggplot(somTime, aes(x=percExtent,y=maxPrecip))+
    geom_point()+
    facet_wrap(~as.factor(codes))+
    ggtitle("Daily Precip Extent vs Max Precip")
  ggplot(somTime, aes(x=meanPrecip,y=spearman))+
    geom_point()+
    facet_wrap(~as.factor(codes))+
    ggtitle("Daily Mean Precip vs Pearson")
    #xlim(0,5)
  
  #
  
  # SOM code vectors
  
  # SOM code vectors distributions
  ggplot(codebook.long, aes(as.factor(codebook.long$codes), value))+
    geom_boxplot(varwidth = TRUE)+
    xlab("Node")+
    ylab("mm")+
    ggtitle("Distribution of Codebook Vectors (precip in mm)")
  
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
  temp<-subset(countsYr, codes=="1_1")
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
                         high = "orange", na.value = NA, name="% of avg",limits=c(0, 2), oob=squish)+
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
    ggtitle("Anom of expected days in each node by MJO (amp>1, *pval<0.05) Phase")+
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
  sumYr2<-somTime %>% 
    group_by(year) %>%
    summarise(sumPrecip=sum(meanPrecip))
  sumYr2$ltAvg<-mean(sumYr2$sumPrecip)
  sumYr2$seasAnom<-sumYr2$sumPrecip-sumYr2$ltAvg
  sumYr2<- merge(sumYr, sumYr2, by='year')
  sumYr2$anomName<-"normal"
    sumYr2$anomName[sumYr2$seasAnom<(-10)] <- "dry"
    sumYr2$anomName[sumYr2$seasAnom>(10)] <- "wet"
  sumYr2 <- sumYr2 %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    sumYr2$row<-as.numeric(sumYr2$row); sumYr2$col<-as.numeric(sumYr2$col); 
  ggplot(sumYr2, aes(x=col, y=-row))+
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
  # counts/anom
    anomCount<-sumYr2 %>% group_by(codes, anomName) %>% summarise(counts=sum(countNodes))
    anomCount <- anomCount %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    anomCount$row<-as.numeric(anomCount$row); anomCount$col<-as.numeric(anomCount$col); 
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=counts))+
      geom_text(aes(label = counts), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 400,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Count of days in each node by Anomaly")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
  # anom of expected counts by anomaly category
    # SOM node counts by anom category
    temp2<-sumYr2 %>% group_by(anomName) %>% summarise(n=sum(countNodes))
      temp2$prop<-temp2$n/sum(temp2$n)
    temp3<-sumYr2 %>% group_by(codes) %>% summarise(n=sum(countNodes))
   
    anomCount<-sumYr2 %>% group_by(codes, anomName) %>% summarise(counts=sum(countNodes))
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
    # ggplot  
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=anomCt))+
      geom_text(aes(label = anomCtLabel), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "lightblue",mid="white", midpoint = 0,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Anom of expected days in each node by Precip Anom (*pval<0.05)")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
    
    
    
    
    
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
  # CREATE COMPOSITES on nodes
  NARR<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2019.grd")
  #NARR<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2019.grd") 
  #NARR<-stack("/scratch/crimmins/NARR/processed/anom/GH500_daily_NARR_WUS_1979_2019_pentMean_anomaly.grd")
  #NARR<-stack("/scratch/crimmins/NARR/processed/anom/PWAT_daily_NARR_WUS_1979_2019_pentMean_anomaly.grd")
  
  # dates - find and remove leap days
  startYr<-1979
  compDates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
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
# RASTERVIS COMPOSITE MAPS
  GH500<-stack("/scratch/crimmins/NARR/processed/00z/GH500_00z_NARR_WUS_1979_2019.grd")
  PWAT<-stack("/scratch/crimmins/NARR/processed/00z/PWAT_00z_NARR_WUS_1979_2019.grd")
  PRCP<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2019.grd")
  # create mean composites
  # dates - find and remove leap days
  startYr<-1979
    compDates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
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
                      main="NARR GH500/Precip Composite Patterns JAS 3x4 SOM (1981-2019)", par.settings=mapTheme)+
              layer(sp.polygons(aznm, col = 'black', lwd=0.5))+
              layer(sp.polygons(countries, col = 'black', lwd=0.5))+
              contourplot(compPWAT,linetype="solid", at=c(25), labels=FALSE, col='darkorange',lwd=0.75)+
              contourplot(compPWAT,linetype="solid", at=c(50), labels=FALSE, col='darkred',lwd=0.75)+
              contourplot(compGH500,linetype="dotdash", at=seq(5700,5850,25), labels=FALSE, col='gray50',lwd=0.5)+
              contourplot(compGH500,linetype="dotted", at=c(seq(5880,6000,5)), labels=FALSE, col='gray28',lwd=0.5)

        png("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/SOM4x3_composites.png", width = 11, height = 8.5, units = "in", res = 300L)
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
 #####       
 
#####
# station metrics
        # bind station datas
        somTime<-cbind.data.frame(somTime,stationDaily)
        
        stationTemp<-melt(somTime[,c(8,9:14)])
        stationTemp$variable = factor(stationTemp$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
        
        ggplot(stationTemp, aes(as.factor(stationTemp$codes), value, color=variable))+
          geom_boxplot(varwidth = FALSE, outlier.shape = 20)+
          ggtitle("Distribution of Daily Station Precip by nodes")
          #ylim(0,40)
           
    # percent of seas total by node
        stationTemp<-melt(somTime[,c(8,9:14)])
          #stationTotal<-stationTemp %>%  group_by(codes, variable) %>% summarise(sumPrecip=sum(value))
          stationPerc<- stationTemp %>%
                  group_by(variable) %>%
                  mutate(sumPrecip=sum(value)) %>%
                  group_by(codes, add=TRUE) %>%
                  summarise(perTotal=sum(value)/min(sumPrecip))
          stationPerc$variable = factor(stationPerc$variable, c("Las Vegas", "Flagstaff", "Phoenix", "Tucson", "Albuquerque", "El Paso"))
          
          
          ggplot(stationPerc, aes(variable,perTotal, fill=variable))+
            geom_bar(stat='identity')+
            facet_wrap(~codes)+
            ylab("avg proportion of seas total")+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank())+
            ggtitle("Average proportion of seasonal total precip by node")
          
          
          
          
  
        
#####        
              
        
#####        
        
  
##### Predict classification of new events - 
  newPrcp<-stack("~/RProjects/SOMs/monsoonPrecip/SWUS_070120_081320_PRISM_daily_prcp.grd")
  # crop to region
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
   predictedUnits<-merge(predictedUnits,code_grid, by="mapUnit")
  