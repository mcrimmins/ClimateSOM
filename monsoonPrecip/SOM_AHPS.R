# SOM for AHPS precipitation 
# adapted from SOM_precip.R
# MAC 9/8/21

library(raster)
library(lubridate)
library(reshape2)
library(kohonen)
# library(tidyr)
# library(ggplot2)
# library(PBSmapping)
# library(dplyr)
# library(cowplot)
# library(scales)
# library(grid)
# library(RColorBrewer)

ptm <- proc.time()

# set rasteroptions
rasterOptions(progress = 'text')

# functions
leap_every_year <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}

# map layers
states <- getData('GADM', country='United States', level=1)
az<-subset(states, NAME_1=="Arizona")
nm<-subset(states, NAME_1=="New Mexico")
#aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico")
aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico"| NAME_1=="California" | NAME_1=="Nevada" | NAME_1=="Utah" | NAME_1=="Colorado"| NAME_1=="Texas")
us <- getData('GADM', country='United States', level=0)
mx <- getData('GADM', country='Mexico', level=0)
#cn <- getData('GADM', country='Canada', level=0)
# ELEVATION grid
elev<-raster("~/RProjects/SOMs/monsoonPrecip/shapes/PRISM_us_dem_4km_asc.asc")

# load PRISM for SW Region - from downloadAHPS_SOMs.R
prcp<- stack("~/RProjects/SOMs/monsoonPrecip/SWUS_010116_123120_AHPS_daily_prcp.grd") 

# dates - find and remove leap days
startYr<-2017
#dates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
dates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2020-12-31"),1))
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
elev<-crop(elev, e)

# convert layers to dataframe
layers.df<-(as.data.frame(subLayers, long=TRUE, xy=TRUE))
colnames(layers.df)<-c("lon","lat","date","value")  
# long to wide
df.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
# save raw precip
#save(df.wide, file="~/RProjects/SOMs/monsoonPrecip/AZwNM_AHPS_JAS_2017_2020.RData")

#load("~/RProjects/SOMs/monsoonPrecip/AZwNM_AHPS_JAS_2017_2020.RData")

# replace NAs with 0's if needed
df.wide[is.na(df.wide)] <- 0

# run SOM with all days or subset based on specific node
idx<-seq(1,nrow(df.wide),1)


  #####
  # CP2 - 2 phase SOM training
  nrows=4
  ncols=4
  ptm <- proc.time()
  #som.gh500 <- som(as.matrix(df.wide[,2:ncol(df.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
  #set.seed(999) #keep set seed 16, 11 upper left/wet, 9,6 upper right dry/UL wet, 8 UL Wet/LR dry, 7 LR Wet/LL Dry, 5 LR wet/LL dry, 4/100 UL wet/UR dry, 101 wet UR/dry LL, 102/104 dry UL/wet UR, 103 UL Wet/LR Dry
  set.seed(123) # 123 for 4x4
  som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
                        #grid = somgrid(ncols, nrows, "rectangular"),
                        grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
                        #alpha = c(0.05, 0.001), # for online
                        radius = c(4,1), #(4,1) for 4x5
                        mode= "pbatch", # "pbatch" or "online"
                        #mode= "online", # "pbatch" or "online"
                        #maxNA.fraction = 0.999,
                        cores = 7,
                        rlen = 500, #5000
                        dist.fcts = "sumofsquares")
  print("completed phase 1, performing phase 2")
  som.gh500.2 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
                        #grid = somgrid(ncols, nrows, "rectangular"),
                        grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
                        #alpha = c(0.05, 0.001), # for online
                        radius = c(2,0.33), # c(3,0.33) for 3x5
                        mode= "pbatch", # "pbatch" or "online"
                        #mode= "online", # "pbatch" or "online"
                        #maxNA.fraction = 0.999,
                        init= som.gh500$codes,
                        cores = 7,
                        rlen = 1500, #7000
                        dist.fcts = "sumofsquares")
      ## quantization error:
      mean(som.gh500.2$distances)
      ## topographical error measures:
      source("topo.error.R")
      topo.error(som.gh500.2, "nodedist")
  proc.time() - ptm
  som.gh500<-som.gh500.2
  #####

 # test<-calc(subLayers, sd)
#  test2<-calc(subLayers, mean)

#   # save SOM output
#    save(som.gh500, file = "~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS_SOM4x4_15K_CP2_1981_2020.RData")
#   # plot(som.gh500, type="changes") # changes, codes, counts, property, quality, mapping

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
  cbRaster[cbRaster ==0] <- NA
  
  
  at<-c(seq(0,20,0.5))
  #at<-c(seq(0,30,2.5))
  mapTheme <- rasterTheme(region = c("lightblue", "blue","green","green4","yellow","red", "red4"))
  pPrecip<-levelplot(cbRaster, contour=FALSE, margin=FALSE, layout=c(ncols,nrows), at=at,
                     par.settings=mapTheme,
                     #par.settings = list(region=c("lightblue", "blue","green","green4","yellow","red", "red4"),
                     #                   axis.line = list(col = textCol[panel.number()])),
                     scales=list(draw=FALSE),
                     main="Precip Patterns JAS 4x4 SOM - AHPS-daily 2017-2020")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    
  
  