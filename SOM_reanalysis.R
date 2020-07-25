# SOM using Reanalysis data downloaded from ESRL
# adapting code from rncep_som_test.R
# MAC 06/11/2020

# to do: adapt to other near-real time gridded datasets like NARR, CFSR

library(raster)
library(lubridate)
library(reshape2)
library(kohonen)
library(tidyr)
library(ggplot2)
library(PBSmapping)

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

# get grid processed from ~/SWMonsoonTracker/NCEPgrids/createNCEPgrids.R, stored in scratch
gh500<-stack("/scratch/crimmins/esrl/processed/HGT500_NCEP_R1_1948_2019.grd")

# crop to NAME region
e <- extent(360-130, 360-90,15, 50)
gh500 <- (crop(gh500, e)) # can also add in rotate to convert lon to neg

# dates - find and remove leap days
nlayers(gh500) 
  dates<-as.data.frame(seq.Date(as.Date("1948-01-01"),as.Date("2019-12-31"),1))
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
  
# convert layers to dataframe
  layers.df<-(as.data.frame(subLayers, long=TRUE, xy=TRUE))
    colnames(layers.df)<-c("lon","lat","date","value")  
    # long to wide
    df.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
    
# kohonen SOM
  nrows=4
  ncols=6
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
    codebook.long<-separate(codebook.long, variable, convert = TRUE, into = c("lat", "lon"), sep="_")
    codebook.long$lat<-as.numeric(gsub("X", "", codebook.long$lat))
    codebook.long$lon<-codebook.long$lon-360
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
    
    ##### MAP RESULTS
    # plot map - fix lines http://cameron.bracken.bz/finally-an-easy-way-to-fix-the-horizontal-lines-in-ggplot2-maps
    # plot limits
    xlim = c(-130,-90)
    ylim = c(15,50)
    
    all_states <- map_data("state")
    world<-map_data("world")
    
    colnames(world)<-c("X","Y","PID","POS","region","subregion")
    world = clipPolys(world, xlim=xlim,ylim=ylim, keepExtra=TRUE)
    #colnames(world)<-c("Lon","Lat","PID","POS","region","subregion") 
    
    colnames(all_states)<-c("X","Y","PID","POS","region","subregion")
    all_states = clipPolys(all_states, xlim=xlim,ylim=ylim, keepExtra=TRUE)
    #colnames(all_states)<-c("Lon","Lat","PID","POS","region","subregion")
    
    p <- ggplot()
    p + geom_polygon( data=world, aes(x=X, y=Y, group = PID),colour="black", fill=NA )+
      scale_x_continuous(breaks = c(-120,-140))+
      stat_contour(data=codebook.long, aes(lon, lat, z=value, color=..level..), size=1)+
      scale_colour_distiller(palette = "Spectral", name="500mb GPH (m)")+
      #scale_color_continuous(low="blue", high = "red")+ 
      coord_map(xlim = xlim,ylim = ylim)+
      facet_wrap(~codes, nrow = nrows, ncol=ncols)+theme_bw()+
      labs(x="Lon", y="Lat")+
      ggtitle("JJAS 500mb GH Pattern Classification (NCEP R1, 1948-2019)")
    
    
  