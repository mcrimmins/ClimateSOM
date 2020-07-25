# test of RNCEP data download
# 5/8/2017
# copied to IE VM on 10/26/17

library(maps)
library(ggplot2)
library(PBSmapping)
all_states <- map_data("state")
world<-map_data("world")

library(RNCEP)
library(tidyr)
library(reshape2)
library(kohonen)
library(ggplot2)
library(dplyr)
library(rasterVis)


# # download data
# lat.southnorth=c(20,50), lon.westeast=c(210,265) for CoL analysis
#lat.southnorth=c(20,50), lon.westeast=c(235,265) for NAMS
lat1<-20
lat2<-50
lon1<-235
lon2<-265
wx.extent1 <- NCEP.gather(variable='hgt', level=500,
                          months.minmax=c(6,9), years.minmax=c(1979,2019),
                          lat.southnorth=c(lat1,lat2), lon.westeast=c(lon1,lon2),
                          reanalysis2 = TRUE, return.units = TRUE, status.bar = FALSE)

# load sample data
#load("~/RProjects/SOMs/JAS500mb_1979-2017.RData")

# map data
NCEP.vis.area(wx.data=wx.extent1, layer='2019-09-20 00', show.pts=TRUE,
              draw.contours=TRUE, cols=terrain.colors(64), transparency=.6,
              title.args=list(main="Example: select layer by datetime"),
              interp.loess.args=list(span=0.5))

# daily mean values
wx.ag <- NCEP.aggregate(wx.data=wx.extent1, YEARS=TRUE, MONTHS=TRUE,
                        DAYS=TRUE, HOURS=FALSE, fxn='mean')

# convert to dataframe
wx.df <- NCEP.array2df(wx.data=wx.ag, var.names='GH500')
# convert to daily anoms based on moving pentad?

# long to wide
wx.wide<-dcast(wx.df, formula = datetime~latitude+longitude, value.var = "GH500")

# melt back to long, need to separate lat/lon column
# test<-melt(wx.wide)

# split date
#test<-separate(wx.df,datetime,convert=TRUE,into = c("year","month","day","hr"))
# add in doy, look at isd_lite for moving window average


# kohonen SOM
nrows=4
ncols=6
som.gh500 <- som(as.matrix(wx.wide[,2:ncol(wx.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
codebook<-as.data.frame(som.gh500$codes)
code_grid<-as.data.frame(som.gh500$grid$pts)
code_grid$mapUnit<-seq(1,nrow(code_grid))
code_grid<-code_grid %>%
  unite(y,x, col="codes", sep="_") # add mapunit back in if needed
# 
# xCols<-as.data.frame(som.gh500$grid$pts[,1])
# colnames(xCols)<-"X_Cols"
# yRows<-as.data.frame(som.gh500$grid$pts[,2])
# colnames(yRows)<-"Y_Rows"
codebook<-cbind(code_grid,codebook)
codebook.long<-melt(codebook, id.vars = 1)
codebook.long<-separate(codebook.long, variable, convert = TRUE, into = c("lat", "lon"), sep="_")
codebook.long$lat<-as.numeric(gsub("X", "", codebook.long$lat))
codebook.long$lon<-codebook.long$lon-360
codebook.long<-separate(codebook.long, codes, convert = FALSE, remove = FALSE, into = c("xCols", "yRows"), sep="_")

# assign days to nodes
nodes<-map(som.gh500)
somTime<-as.data.frame(cbind(wx.wide$datetime, nodes$unit.classif, nodes$distances))
colnames(somTime)<-c("datetime","mapUnit","errorDist")
somTime<-separate(somTime, datetime, convert = TRUE, remove = FALSE, into = c("year","month","day","hr"), sep="_")
somTime$date<-as.Date(paste(somTime$year,"-",somTime$day,"-",somTime$month,sep=""), format="%Y-%d-%m")
somTime$doy<-as.numeric(format(somTime$date, "%j"))
somTime$mapUnit<-as.integer(somTime$mapUnit)
somTime$errorDist<-as.numeric(as.character(somTime$errorDist))


# plot map - fix lines http://cameron.bracken.bz/finally-an-easy-way-to-fix-the-horizontal-lines-in-ggplot2-maps
# plot limits
xlim = c(-125,-95)
ylim = c(20,50)

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
  ggtitle("Sept 500mb GH Pattern Classification (NCEP R2, 1979-2018)")

# summary plots -- appears to plot opposite up/down from SOM plot
counts <- plot(som.gh500, type="counts", shape = "straight", labels=counts)
codes <- plot(som.gh500, type="codes", shape = "straight")
similarities <- plot(som.gh500, type="quality", palette.name = terrain.colors)
plot(som.gh500, type="dist.neighbours", main = "SOM neighbour distances")

# sammon mapping
library(MASS)
gh500.codes <- som.gh500$codes
dis <- dist(as.matrix(som.gh500$codes[[1]]))
gh500.sam <- sammon(dis)
plot(gh500.sam$points, type="n")
text(gh500.sam$points,labels=as.character(1:nrow(code_grid)))


# join codes to node table
somTime<-left_join(somTime, code_grid)
# plot map units
ggplot(somTime, aes(doy, year)) + 
  geom_tile(aes(fill = mapUnit), colour = "grey") + 
  geom_text(aes(label = codes), size=1.50)+
  scale_fill_gradient2(low = "red", mid = "green",
                       high = "blue", midpoint = 7, space = "Lab",
                       na.value = "grey50", guide = "colourbar")
# plot error
ggplot(somTime, aes(doy, year)) + 
  geom_tile(aes(fill = errorDist), colour = "grey") + 
  scale_fill_gradient2(low = "white", mid = "yellow",
                       high = "red", space = "Lab",
                       na.value = "grey50", guide = "colourbar")


# plot precip with local CPC precip


# # plot precip with rnoaa cpc_prcp
# library(rnoaa)
# i=100
# tmp<-cpc_prcp(date = somTime$date[i])
# tmp$lon<-tmp$lon-180
# 
# for(i in 1:4){
#   tmp<-cpc_prcp(date = somTime$date[i])
#   if (i==1){
#     tmp2 <- tmp
#   }else{
#     tmp2 <-bind_cols(tmp2,tmp[,3]) # brick or stack?
#   }
# }
# 
# tmp<-cpc_prcp(date = "2019-08-16", us=FALSE)
# p <- ggplot(tmp, aes(x=lon, y=lat, fill=precip)) + theme_bw()
# p + geom_tile()+
#   ggtitle("cpc_prcp(date = '2019-08-16', us=FALSE)")
# 
# p <- ggplot(tmp, aes(x=lon, y=lat)) + theme_bw()
# p + geom_raster(aes(fill=precip))

# plot CPC precip from netcdf files
library(ncdf4)
library(raster)

yr1<-1979
yr2<-2018

allPrecip <- stack()
for(i in yr1:yr2){
  cpc.prcp.file <- paste0("/scratch/crimmins/cpc_global/ftp.cdc.noaa.gov/Datasets/cpc_global_precip/precip.",i,".nc")
  #"load" data into R via Raster package
  # cycle through each year, assigning each day to appropriate SOM unit
  prcp.tmp <- brick(cpc.prcp.file, varname="precip",  ncdf=TRUE)
  #plot(prcp.tmp[[1]])

  # grep out month(s) -- set months 
  prcp.tmp <- raster::subset(prcp.tmp, grep('.09.', names(prcp.tmp), value = T, fixed = T)) # use pipe for more months .08.|.09.
  
  # crop extent
  e <- extent(lon1, lon2, lat1, lat2)
  prcp.tmp <- crop(prcp.tmp, e)
  # build stack
  allPrecip <- stack( allPrecip , prcp.tmp)  
  print(i)
}

# get composite means
precipComposite<-stack()
codePatterns<-stack()
codebook.trim<-codebook.long[-c(1:(ncols*nrows)),]

for(i in 1:nrow(code_grid)){
  # build precip composite
  temp<-overlay(subset(allPrecip, which(somTime$mapUnit==i)), fun=function(x){mean(x,na.rm=T)})
  names(temp)<-paste0(code_grid$codes[i]," (n=",length(which(somTime$mapUnit==i)),")")
  precipComposite<-stack(precipComposite,temp)
  # build height pattern raster
  codeTemp<-subset(codebook.trim, codes==code_grid$codes[i])
  codeTemp<-rasterFromXYZ(codeTemp[,c(5,4,6)])
  codePatterns<-stack(codePatterns,codeTemp)
  print(i)
}

precipComposite[precipComposite == 0] <- NA
precipComposite<-raster::rotate(precipComposite)
names(codePatterns)<-seq(1,ncols*nrows,1)

# mapping
statesPoly <- map_data("state")
worldPoly<-map_data("world")
# create color scheme for precip
#mycols <- rasterTheme(region=colorRampPalette(brewer.pal(9,'YlGnBu'))(100))
mycols <- rasterTheme(region=colorRampPalette(rev(brewer.pal(9,'Spectral')))(100))
cols.at <- seq(0, 15, 0.5)
levs.at<-seq(5600,6000, 25)
levelplot(precipComposite, margin=FALSE, layout=c(ncols,nrows),par.settings = mycols,
          at=cols.at, main="Sept 500mb GH Pattern Classification and mean precip (NCEP R2, 1979-2019)")+
          contourplot(codePatterns,at=levs.at,
              lwd = 0.4,
              labels = list(cex = 0.6),
              label.style = 'align')
  #layer(sp.polygons(statesPoly, col = 'gray40', lwd=0.3))

# plot paths of days/transitions on grid for specified time periods
somTime<-separate(data = somTime, col = codes, into = c("codeRow", "codeCol"), sep = "_", convert = TRUE)

ggplot(somTime[somTime$date>="2017-08-01" & somTime$date<="2017-08-31",], aes(codeCol, codeRow))+
  geom_path(aes(codeCol,codeRow, color=doy))+
  scale_y_reverse(lim=c(nrows+1,0))+
  scale_color_distiller(type = "seq", palette = "YlGnBu")+
  xlim(0,ncols+1)+
  theme_bw()

# grab data from ESRL for near real time