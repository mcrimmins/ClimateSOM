# PRISM Daily Download script for SOM - recent precip for classification
# 07/20/20 MAC
# see http://ropensci.github.io/prism/ for more info
# adapted from dailyDownloadPRISM.R

library(prism)
library(raster)

# set rasteroptions
rasterOptions(progress = 'text')

# change to data directory
options(prism.path = "~/RProjects/SOMs/monsoonPrecip/recentPRISM/") 

# download tmean daily data
get_prism_dailys(type="ppt", minDate = "2020-08-01", maxDate = "2020-8-31", keepZip=FALSE)

### NEED TO FIX FILES WHEN THE CHANGE TO FINAL

# Work with PRISM rasters ----
# get date order, make sure in right order
names<-ls_prism_data(name=TRUE)
dates<-as.Date(strptime(substr(names$product_name, 1,11),"%b %d %Y"))
dateOrder<-as.data.frame(sort.int(dates, index.return=TRUE))

# get into stack
precip_stack<-stack()
yr1<-2020
yr2<-2020
for(i in yr1:yr2){
  paste0(i,"-01-01")
  subsetDateX<-dateOrder$ix[which(dateOrder$x >= paste0(i,"-01-01") & dateOrder$x <= paste0(i,"-12-31"))]
  pcpStack <- prism_stack(names$files[subsetDateX])
  
  pcpStack<-crop(pcpStack, extent(-119,-100,28,43))
  
  precip_stack<-stack(precip_stack, pcpStack) 
  
  print(i)
}  

# name layers
names(precip_stack)<-seq.Date(as.Date("2020-07-01"),as.Date("2020-8-31"),1)
# write out raw yearly values
writeRaster(precip_stack,filename="~/RProjects/SOMs/monsoonPrecip/SWUS_070120_083120_PRISM_daily_prcp.grd", overwrite=TRUE)
