# grab monthly PRISM precip from RCC-ACIS, plot monsoon precip
# MAC 11/15/21

# get monthly PRISM from RCC ACIS
library(RCurl)
library(jsonlite)
library(raster)    
# create current date
dateRangeStart="1900-01-01"
dateRangeEnd= "2021-12-31"

# custom functions
perc.rank<-function(x) trunc(rank(x,ties.method = "average"))/length(x)

# generate dates -- keep with PRISM date
allDates<-seq(as.Date(dateRangeStart), as.Date(dateRangeEnd),by="month")

# Set bounding box for PRISM extract
# AZ/NM bbox -115.004883,31.184609,-102.524414,37.387617
ACISbbox<-"-115.5,31.3,-106,37.5" 

# ACIS query in JSON
jsonQuery=paste0('{"bbox":"',ACISbbox,'","sdate":"',dateRangeStart,'","edate":"',dateRangeEnd,'","grid":"21","elems":"mly_pcpn","meta":"ll,elev","output":"json"}') # or uid
#jsonQuery=paste0('{"bbox":"',ACISbbox,'","sdate":"',dateRangeStart,'","edate":"',dateRangeEnd,'","grid":"2","elems":"pcpn","meta":"ll","output":"json"}') # or uid

out<-postForm("http://data.rcc-acis.org/GridData",
              .opts = list(postfields = jsonQuery,
                           httpheader = c('Content-Type' = 'application/json', Accept = 'application/json')))
out<-fromJSON(out)

# convert to list of matrices, flipud with PRISM
matrixList <- vector("list",length(out$data))
for(i in 1:length(out$data)){
  matrixList[[i]]<-apply(t(out$data[[i]][[2]]),1,rev)
}

# read into raster stack
rasterList<-lapply(matrixList, raster)
gridStack<-stack(rasterList)
gridExtent<-extent(min(out$meta$lon), max(out$meta$lon), min(out$meta$lat), max(out$meta$lat))
gridStack<-setExtent(gridStack, gridExtent, keepres=FALSE, snap=FALSE)
names(gridStack)<-allDates
# set 0 and neg to NA
gridStack[gridStack < 0] <- NA
## manage dates
allDates<-as.data.frame(allDates)
allDates$month<-as.numeric(format(allDates$allDates, "%m"))
allDates$year<-as.numeric(format(allDates$allDates, "%Y"))
# set months for monsoon season
idx<-which(allDates$month %in% c(7,8,9)) # grab only summer months
#idx<-which(allDates$month %in% c(7)) # grab only summer months
allDates<-allDates[idx,]
gridStack<-gridStack[[idx]]

# calculate seasonal totals
sumSeas<-stackApply(gridStack, allDates$year, fun = sum)
seasAvgPrecip<-cellStats(sumSeas, 'mean')
seasAvgPrecip<-cbind.data.frame(unique(allDates$year),seasAvgPrecip)
seasAvgPrecip$percRank<-perc.rank(seasAvgPrecip$seasAvgPrecip) 
colnames(seasAvgPrecip)<-c("year","avgPrecip","percRank")
# names for terciles
seasAvgPrecip$anomName<-"normal"
seasAvgPrecip$anomName[seasAvgPrecip$percRank<=0.33] <- "dry"
seasAvgPrecip$anomName[seasAvgPrecip$percRank>=0.66] <- "wet"

library(cowplot)
library(ggplot2)
ggplot(seasAvgPrecip, aes(year,avgPrecip, fill=as.factor(seasAvgPrecip$anomName)) )+
  geom_bar(stat = 'identity')+
  ggtitle("Regional Average Total Precip (July-Aug-Sept)")+
  geom_hline(yintercept=mean(seasAvgPrecip$avgPrecip), color="black")+
  geom_hline(yintercept=median(seasAvgPrecip$avgPrecip), color="red")+
  scale_fill_manual(values = c("saddlebrown", "grey", "forestgreen"), name="tercile")+
  ylab("inches")+
  theme_bw()

write.csv(seasAvgPrecip, file="AZNM_JAS_1900_2021.csv", row.names = FALSE)
