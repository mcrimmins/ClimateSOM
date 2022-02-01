# get station data for SOM analysis
# MAC 10/8/20

library(RCurl)
library(jsonlite)
library(lubridate)

# functions
leap_every_year <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}


# Station IDS
# Tucson 028820
# Phoenix 026481
# Flagstaff 023010
# Las Vegas 264436
# El Paso 412797
# Albuquerque 290234

jsonQuery='{"sids":"028820,026481,023010,264436,412797,290234","sdate":"1950-01-01","edate":"2020-12-31","elems":"4"}'

out<-postForm("http://data.rcc-acis.org/MultiStnData", 
              .opts = list(postfields = jsonQuery, 
                           httpheader = c('Content-Type' = 'application/json', Accept = 'application/json')))
out<-fromJSON(out)

# get data from list
ll<-data.frame(matrix(unlist(out$data$meta$ll), nrow=length(out$data$meta$ll), byrow=T))
meta<-out$data$meta
  data<- t(matrix(unlist(out$data$data), nrow=length(out$data$data), byrow=T))
  data[data=="T"]<-0.001
data<-cbind.data.frame(seq(as.Date("1950-01-01",format="%Y-%m-%d"),as.Date("2020-12-31",format="%Y-%m-%d"),1),data)
# convert columns to numeric
unfactorize<-c(2:7)
data[,unfactorize]<-lapply(unfactorize, function(x) as.numeric(as.character(data[,x])))  
colnames(data)<-c("date",meta$name)
# PRISM data correction
data$date<-data$date+1

# subset to July-Aug-Sep 81-19
data$year<-as.numeric(format(data$date,"%Y"))
data$month<-as.numeric(format(data$date,"%m"))
data$doy_ly<-leap_every_year(data$date) 
data<-subset(data, year>=1981)
#data<-data[which(data$month %in% c(7,8,9)),]
data<-data[which(data$doy_ly %in% seq(167,274,1)),]

save(data, file="stationPrecip_thru2020.RData")



  
  