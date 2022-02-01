# download AHPS for SOMs
# MAC 09/08/21

library(raster)

# load mask
mask<-raster( "~/RProjects/SOMs/monsoonPrecip/SWUS_PRISM_MASK.grd")

# get PRISM ref grid
newPrcp<-stack("~/RProjects/SOMs/monsoonPrecip/SWUS_070121_083021_PRISM_daily_prcp.grd")
# crop to region
e <- extent(-115.5,-106,31.3, 37.5) 
newPrcp <-crop(newPrcp[[1]],e)


# generate dates -- keep with PRISM date
allDates<-seq(as.Date("2017-01-01"), as.Date("2020-12-31"),1)
#####
# download daily netcdf files from NOAA
ahpsStack<-stack()
for(i in 1:length(allDates)){
  #build URL
  URL<-paste0("https://water.weather.gov/precip/downloads/",
              format(allDates[i],"%Y"),"/",format(allDates[i],"%m"),"/",format(allDates[i],"%d"),
              "/nws_precip_1day_",format(allDates[i],"%Y%m%d"),"_conus.nc")
  download.file(URL, destfile = paste0("/home/crimmins/RProjects/SOMs/monsoonPrecip/temp/",allDates[i],".nc"), method="curl")
  temp<-raster(paste0("/home/crimmins/RProjects/SOMs/monsoonPrecip/temp/",allDates[i],".nc"), varname="observation")
  ahpsStack <- stack(ahpsStack,temp)
  print(allDates[i])
  print(URL)
}
# clean up files in temp dir

# use cbRaster for climo/reference grid
#temp<-cbRaster[[1]]
# reproject AHSP to lat/lon
proj4string(newPrcp) <- "+proj=longlat +datum=WGS84 +no_defs"
ahpsStack <- projectRaster(ahpsStack, newPrcp)

# crop to region
prcp<-ahpsStack*25.4
#e <- extent(-115.5,-106,31.3, 37.5) 
prcp <-crop(prcp,newPrcp)
# apply mask
#prcp <- mask(prcp, mask)

# name layers
names(prcp)<-allDates
# write out raw yearly values
writeRaster(prcp,filename="~/RProjects/SOMs/monsoonPrecip/SWUS_010116_123120_AHPS_daily_prcp.grd", overwrite=TRUE)

