# create PRISM climo maps for monsoon SOM study
# MAC 08/05/2020

library(prism)
library(raster)
library(rasterVis)
library(maptools)

# set rasteroptions
rasterOptions(progress = 'text')

# change to data directory
options(prism.path = "~/RProjects/SOMs/monsoonPrecip/climoPRISM/") 

#get_prism_normals(type="ppt",resolution = "4km",mon = 1:12, keepZip=F)

# build stack
names<-ls_prism_data(name=TRUE)
moPPT<-prism_stack(names$files)
percJAS<-calc(moPPT[[7:9]],sum)/calc(moPPT,sum)

plot(moPPT[[7:9]])

# study area box
x_coord <- c(-115.5,-115.5,-106,-106)
y_coord <- c(31.3,37.5,37.5,31.5)
xym <- cbind(x_coord, y_coord)
library(sp)
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))

# rastervis maps
# map layers
states <- getData('GADM', country='United States', level=1)
mlra <- readShapePoly(paste0("/home/crimmins/RProjects/LivnehDrought/shapes/mlra/mlra_v42.shp"))
# read in HUC4
#rgdal::ogrListLayers("~/RProjects/SOMs/monsoonPrecip/shapes/wbdhu4_a_us_september2019/wbdhu4_a_us_september2019.gdb")
huc4<-rgdal::readOGR(dsn = "~/RProjects/SOMs/monsoonPrecip/shapes/wbdhu4_a_us_september2019/wbdhu4_a_us_september2019.gdb", layer="WBDHU4")
# intersect of HUC and study area
huc4clip<-raster::intersect(huc4, sps)
# save clipped shapefile
rgdal::writeOGR(huc4clip,dsn="~/RProjects/SOMs/monsoonPrecip/shapes", layer="huc4clip",driver="ESRI Shapefile")

at <- seq(0, 65, 1)
p0 <- levelplot(percJAS*100, par.settings = YlOrRdTheme, ylab=NULL, xlab=NULL, margin=FALSE,
                at=at, ylim=c(30,40), xlim=c(-120,-100), main="JAS Percent of Annual Total")+ # width=1, height=0.5, row=3, column=1, 
  layer(sp.polygons(states, col = 'gray40', lwd=0.5))+
  layer(sp.polygons(huc4clip, col = 'gray20', lwd=0.5))+
  layer(sp.polygons(sps, col = 'black', lwd=1))


at <- c(seq(0, 150, 5),180)  
p0 <- levelplot(moPPT[[6:10]], par.settings = viridisTheme, ylab=NULL, xlab=NULL, margin=FALSE,
                names.attr=c("June","July","August","September","October"), at=at, layout=c(1,5),
                ylim=c(30,40), xlim=c(-120,-100), main="Monthly Total Precip (mm)")+ # width=1, height=0.5, row=3, column=1, 
  layer(sp.polygons(states, col = 'gray40', lwd=0.5))+
  layer(sp.polygons(sps, col = 'black', lwd=1))





