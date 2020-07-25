# SOM using Reanalysis data downloaded from ESRL
# adapting code from SOM_reanalysis.R
# MAC 06/25/2020

library(raster)
library(lubridate)
library(reshape2)
library(kohonen)
library(tidyr)
library(ggplot2)
library(PBSmapping)
library(dplyr)


ptm <- proc.time()

# set rasteroptions
rasterOptions(progress = 'text')

# functions
leap_every_year <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}

# map layers
states <- getData('GADM', country='United States', level=1)
us <- getData('GADM', country='United States', level=0)
mx <- getData('GADM', country='Mexico', level=0)
#cn <- getData('GADM', country='Canada', level=0)

# # load pre-processed raster stacks from /ClimPlot/NARR/processNARR.R
gh500<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2019.grd")
# pwat<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2019.grd") 
prcp<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2019.grd") 
# avgCAPE<-stack("/scratch/crimmins/NARR/processed/CAPE_daily_NARR_WUS_1979_2019.grd")
# maxCAPE<-stack("/scratch/crimmins/NARR/processed/maxCAPE_daily_NARR_WUS_1979_2019.grd")
# sfcDP<-stack("/scratch/crimmins/NARR/processed/DPT2m_daily_NARR_WUS_1979_2019.grd")
# uFlux<-stack("/scratch/crimmins/NARR/processed/WVUFLX_daily_NARR_WUS_1979_2019.grd")
# vFlux<-stack("/scratch/crimmins/NARR/processed/WVVFLX_daily_NARR_WUS_1979_2019.grd")

# crop to SW/NMX region
#e <- extent(-120,-98,22, 39)
#gh500 <- (crop(gh500, e)) # can also add in rotate to convert lon to neg

# dates - find and remove leap days
dates<-as.data.frame(seq.Date(as.Date("1979-01-01"),as.Date("2019-12-31"),1))
  colnames(dates)<-"date"
  dates$month<-as.numeric(format(dates$date, "%m"))
  dates$day<-as.numeric(format(dates$date, "%d"))
  dates$year<-as.numeric(format(dates$date, "%Y"))
  dates$doy<-as.numeric(format(dates$date, "%j"))
  dates$doy_ly<-leap_every_year(dates$date) # this day of year without leap day shifts

  # subset layers to months of interest
  mos<-c(6,7,8,9)
  subDates<-dates[which(dates$month %in% mos),]
  subLayers<-gh500[[which(dates$month %in% mos)]]
  subLayers2<-prcp[[which(dates$month %in% mos)]]
  
  # convert layers to dataframe
  layers.df<-(as.data.frame(subLayers, long=TRUE, xy=TRUE))
  colnames(layers.df)<-c("lon","lat","date","value")  
  # long to wide
  df.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
  
  # convert summary layers to dataframe
  layers.df<-(as.data.frame(subLayers2, long=TRUE, xy=TRUE))
  colnames(layers.df)<-c("lon","lat","date","value")  
  # long to wide
  #df2.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
  layers.df$date<-as.Date(layers.df$date, format="X%Y.%m.%d")
  
  
  # kohonen SOM
  nrows=5
  ncols=7
  som.gh500 <- som(as.matrix(df.wide[,2:ncol(df.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
  codebook<-as.data.frame(som.gh500$codes)
  code_grid<-as.data.frame(som.gh500$grid$pts)
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
  codebook.long<-separate(codebook.long, codes, convert = FALSE, remove = FALSE, into = c("xCols", "yRows"), sep="_")
  
  # assign days to nodes
  nodes<-map(som.gh500)
  somTime<-as.data.frame(cbind(df.wide$date, nodes$unit.classif, nodes$distances))
  colnames(somTime)<-c("date","mapUnit","errorDist")
  somTime$date<-as.character(somTime$date)
  somTime$date<-as.Date(somTime$date, format="X%Y.%m.%d")
  somTime$month<-as.numeric(format(somTime$date, format="%m"))
  somTime$year <-as.numeric(format(somTime$date, format="%Y"))
  somTime$day  <-as.numeric(format(somTime$date, format="%d"))
  #somTime<-separate(somTime, as.character(date), convert = TRUE, remove = FALSE, into = c("year","month","day"), sep=".")
  #somTime$date<-as.Date(paste(somTime$year,"-",somTime$day,"-",somTime$month,sep=""), format="%Y-%d-%m")
  somTime$doy<-as.numeric(format(somTime$date, "%j"))
  somTime$mapUnit<-as.integer(somTime$mapUnit)
  somTime$errorDist<-as.numeric(as.character(somTime$errorDist))

  # join nodes/dates to target summary layer
  nodeDates<-somTime[,c("date","mapUnit")]
  layers.df<-merge(layers.df, nodeDates, by="date")
  targetGrid <- layers.df %>% group_by(lon, lat, mapUnit) %>% summarise(medPrec=median(value, na.rm=TRUE))
  targetGrid <- merge(targetGrid, mapunits, by.x="mapUnit", by.y="value")

proc.time() - ptm  
  
  # RASTERVIS mapping of SOM results
  library(rasterVis)
  # create stack of codebook results
  codeList<-unique(codebook.long$codes)
  cbRaster<-stack()
  for (i in 1:length(codeList)) {
    temp<-codebook.long[which(codebook.long$codes==codeList[i]),c(5,4,6)]
    temp<-rasterFromXYZ(temp)
    cbRaster<-stack(cbRaster, temp)
  }
names(cbRaster)<-codeList
  # precip grids
  pptRaster<-stack()
  for (i in 1:length(codeList)) {
    temp<-targetGrid[which(targetGrid$codes==codeList[i]),c(2,3,4)]
    temp<-rasterFromXYZ(temp)
    pptRaster<-stack(pptRaster, temp)
  }
names(pptRaster)<-codeList

# https://gis.stackexchange.com/questions/121306/how-to-reproduce-contour-style-labeling-with-rasterviscontourplot
pptRaster[pptRaster == 0] <- NA
at=c(seq(0.01,8,1),10)
levels<-seq(5500,6000, 25)
p<- levelplot(pptRaster, contour=FALSE, margin=FALSE, par.settings=viridisTheme, at=at,layout=c(ncols,nrows),
              main="GH500 JJAS 5x7 SOM - NARR 1979-2019")+
    contourplot(cbRaster,linetype = "dashed", pretty=TRUE, at=levels,
                labels = list(cex = 0.4),
                label.style = 'align')+
    layer_(sp.polygons(us, col = 'gray40', lwd=1, fill='gray'))+
    layer_(sp.polygons(mx, col = 'gray40', lwd=1, fill='gray'))
png("/home/crimmins/RProjects/SOMs/GH500.png", width = 10, height = 6, units = "in", res = 300L)
#grid.newpage()
print(p, newpage = FALSE)
dev.off() 
  
  
  
    
  ##### MAP RESULTS ----
  # plot map - fix lines http://cameron.bracken.bz/finally-an-easy-way-to-fix-the-horizontal-lines-in-ggplot2-maps
  # plot limits
  xlim = c(-120,-90)
  ylim = c(20,40)
  
  all_states <- map_data("state")
  world<-map_data("world")
  
  colnames(world)<-c("X","Y","PID","POS","region","subregion")
  world = clipPolys(world, xlim=xlim,ylim=ylim, keepExtra=TRUE)
  #colnames(world)<-c("Lon","Lat","PID","POS","region","subregion") 
  
  colnames(all_states)<-c("X","Y","PID","POS","region","subregion")
  all_states = clipPolys(all_states, xlim=xlim,ylim=ylim, keepExtra=TRUE)
  #colnames(all_states)<-c("Lon","Lat","PID","POS","region","subregion")
  
  p <- ggplot()
 p<- p + geom_contour_filled(data=targetGrid, aes(lon, lat, z=medPrec, color=medPrec))+
     geom_polygon( data=world, aes(x=X, y=Y, group = PID),colour="black", fill=NA )+
    scale_x_continuous(breaks = c(-120,-140))+
    stat_contour(data=codebook.long, aes(-lon, lat, z=value, color=..level..), size=1)+
    scale_colour_distiller(palette = "Spectral", name="500mb GPH (m)")+
    #scale_color_continuous(low="blue", high = "red")+ 
    coord_map(xlim = xlim,ylim = ylim)+
    facet_wrap(~codes, nrow = nrows, ncol=ncols)+theme_bw()+
    labs(x="Lon", y="Lat")+
    ggtitle("JJAS 500mb GH Pattern Classification (NARR, 1979-2019)")
  
  ggsave("GH500.png", plot=p, width = 12, height = 10, units = "in")
# ----    