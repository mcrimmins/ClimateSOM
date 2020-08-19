# FLASH database data for SOMs
# https://inside.nssl.noaa.gov/flash/database/
# MAC 08/14/2020

# spatial downloads - nws points
download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/nws_ff_pts_shp.zip",
              destfile="monsoonPrecip/shapes/nws_pts.zip")
unzip ("monsoonPrecip/shapes/nws_pts.zip", exdir = "./monsoonPrecip/shapes/")
library(rgdal)
nwspts <- readOGR( 
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/") , 
  layer="points",
  verbose=FALSE
)

# spatial downloads - nws polygons
download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/nws_ff_poly_shp.zip",
              destfile="monsoonPrecip/shapes/nws_poly.zip")
unzip ("monsoonPrecip/shapes/nws_poly.zip", exdir = "./monsoonPrecip/shapes/")

nwspoly <- readOGR( 
  dsn= paste0(getwd(),"/monsoonPrecip/shapes/") , 
  layer="polygons",
  verbose=FALSE
)


# spatial downloads - usgs points
download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/usgs_events_shp.zip",
              destfile="monsoonPrecip/shapes/usgs_pts.zip")
unzip ("monsoonPrecip/shapes/usgs_pts.zip", exdir = "./monsoonPrecip/shapes/")

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

# SHAVE event data
download.file("http://www.nssl.noaa.gov/projects/flash/database/2016v1/shave_shp.zip",
              destfile="monsoonPrecip/shapes/shave.zip")
unzip ("monsoonPrecip/shapes/shave.zip", exdir = "./monsoonPrecip/shapes/")
  unzip ("monsoonPrecip/shapes/shave_shp/shave_hu15.zip", exdir = "./monsoonPrecip/shapes/shp15")
  #unzip ("monsoonPrecip/shapes/shave_shp/shave_hu14.zip", exdir = "./monsoonPrecip/shapes/shp14")
  unzip ("monsoonPrecip/shapes/shave_shp/shave_hu13.zip", exdir = "./monsoonPrecip/shapes/shp13")
  
  
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


