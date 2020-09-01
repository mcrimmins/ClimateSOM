# FLASH database data for SOMs
# https://inside.nssl.noaa.gov/flash/database/
# MAC 08/14/2020

library(sp)
library(maptools)
library(raster)

# som extent
e <- extent(-115.5,-106,31.3, 37.5) 

# spatial downloads - nws points
# download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/nws_ff_pts_shp.zip",
#               destfile="monsoonPrecip/shapes/nws_pts.zip")
# unzip ("monsoonPrecip/shapes/nws_pts.zip", exdir = "./monsoonPrecip/shapes/")
library(rgdal)
nwspts <- readOGR( 
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/") , 
  layer="points",
  verbose=FALSE
)
# crop to extent
nwspts<-crop(nwspts,e)
# dates
nwspts$month<-as.numeric(format(as.POSIXct(as.numeric(as.character(nwspts@data$BegUnixTi)), origin="1970-01-01"), "%m"))
nwspts$year<-as.numeric(format(as.POSIXct(as.numeric(as.character(nwspts@data$BegUnixTi)), origin="1970-01-01"), "%Y"))
nwspts$day<-as.numeric(format(as.POSIXct(as.numeric(as.character(nwspts@data$BegUnixTi)), origin="1970-01-01"), "%d"))
nwspts$date<-as.Date(paste0(nwspts$month,"-",nwspts$day,"-",nwspts$year),format="%m-%d-%Y")
# subset
nwspts <- nwspts[nwspts$month %in% c(7,8,9),]
nwspts <- nwspts[nwspts$year %in% c(1981:2019),]


# spatial downloads - nws polygons
# download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/nws_ff_poly_shp.zip",
#               destfile="monsoonPrecip/shapes/nws_poly.zip")
# unzip ("monsoonPrecip/shapes/nws_poly.zip", exdir = "./monsoonPrecip/shapes/")

# nwspoly <- readOGR( 
#   dsn= paste0(getwd(),"/monsoonPrecip/shapes/") , 
#   layer="polygons",
#   verbose=FALSE
# )
# # crop to extent
# nwspoly<-crop(nwspoly,e)


# spatial downloads - usgs points
# download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/usgs_events_shp.zip",
#               destfile="monsoonPrecip/shapes/usgs_pts.zip")
# unzip ("monsoonPrecip/shapes/usgs_pts.zip", exdir = "./monsoonPrecip/shapes/")

huc15 <- readOGR(
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/") ,
  layer="huc15",
  verbose=FALSE
)
huc14 <- readOGR(
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/") ,
  layer="huc14",
  verbose=FALSE
)
huc13 <- readOGR(
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/") ,
  layer="huc13",
  verbose=FALSE
)
hucs<-rbind(rbind(huc13, huc14), huc15)
# crop to extent
hucs<-crop(hucs, e)
# times 
hucs$Start.Time<-as.POSIXct(as.character(hucs$Start.Time))
hucs$End.Time..<-as.POSIXct(as.character(hucs$End.Time..))
hucs$Peak.Time<-as.POSIXct(as.character(hucs$Peak.Time))
# month/date/time
hucs$month<-as.numeric(format(hucs$Peak.Time,"%m"))
hucs$day<-as.numeric(format(hucs$Peak.Time,"%d"))
hucs$year<-as.numeric(format(hucs$Peak.Time,"%Y"))
# subset
hucs <- hucs[hucs$month %in% c(7,8,9),]
hucs <- hucs[hucs$year %in% c(1981:2019),]
hucs$date<-as.Date(format(hucs$Peak.Time,"%Y-%m-%d"))
save(hucs, file="~/RProjects/SOMs/monsoonPrecip/USGS_FFlood.RData")


# SHAVE event data
# download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/shave_shp.zip",
#               destfile="monsoonPrecip/shapes/shave.zip")
# unzip ("monsoonPrecip/shapes/shave.zip", exdir = "./monsoonPrecip/shapes/")
#   unzip ("monsoonPrecip/shapes/shave_shp/shave_hu15.zip", exdir = "./monsoonPrecip/shapes/shp15")
#   #unzip ("monsoonPrecip/shapes/shave_shp/shave_hu14.zip", exdir = "./monsoonPrecip/shapes/shp14")
#   unzip ("monsoonPrecip/shapes/shave_shp/shave_hu13.zip", exdir = "./monsoonPrecip/shapes/shp13")
#   
  
huc15shave <- readOGR( 
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/shp15") , 
  layer="shave_hu15",
  verbose=FALSE
)
huc13shave <- readOGR( 
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/shp13") , 
  layer="shave_hu13",
  verbose=FALSE
)
hucShave<-rbind(huc13shave,huc15shave)
# crop to extent
hucShave<-crop(hucShave, e)
