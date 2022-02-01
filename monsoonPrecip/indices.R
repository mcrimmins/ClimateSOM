# get indices for SOM analysis
# MAC 07/29/2020

# get MJO from BOM
# http://www.bom.gov.au/climate/mjo/graphics/rmm.74toRealtime.txt

library(RCurl)
library(reshape2)

# MJO
indexFile<-getURL("http://www.bom.gov.au/climate/mjo/graphics/rmm.74toRealtime.txt")
MJO <- read.table(textConnection(indexFile), skip=2,header=F, sep="")
colnames(MJO)<-c("year", "month", "day", "RMM1", "RMM2", "phase", "amplitude","note")

# ONI from ESRL
# https://psl.noaa.gov/data/correlation/oni.data
indexFile<-getURL("https://psl.noaa.gov/data/correlation/oni.data")
ONI <- read.table(textConnection(indexFile), skip=1, header=F, sep="", nrows=length(seq(1950,2020,1)))
colnames(ONI)<-c("year", seq(1,12,1))

# AO from ESRL
# https://psl.noaa.gov/data/correlation/ao.data
# indexFile<-getURL("https://psl.noaa.gov/data/correlation/ao.data")
# AOI <- read.table(textConnection(indexFile), skip=1, header=F, sep="", nrows=length(seq(1950,2020,1)))
# colnames(ONI)<-c("year", seq(1,12,1))


# combine into time series for SOM analysis 1981-2019
ONI<-melt(ONI, measure.vars = 2:13)
ONI$moDate<-as.Date(paste0(ONI$variable,"-01-",ONI$year),format = "%m-%d-%Y")
MJO$moDate<-as.Date(paste0(MJO$month,"-01-",MJO$year),format = "%m-%d-%Y")
MJO$date<-as.Date(paste0(MJO$month,"-",MJO$day,"-",MJO$year),format = "%m-%d-%Y")

# combine on date
climInd<-merge(MJO,ONI, by="moDate")
colnames(climInd)[13]<-"ONI"

save(climInd, file = "~/RProjects/SOMs/monsoonPrecip/climInd.RData")

# trim to SOM time series
# somInd<-subset(climInd, date>="1981-01-01" & date<="2019-12-31")