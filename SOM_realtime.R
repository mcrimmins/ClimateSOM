# real time SOM mapping for monsoon precip
# MAC 07/25/2022

library(RCurl)
library(jsonlite)
library(raster)
library(kohonen)
library(tidyr)
library(ggplot2)
library(lubridate)


#####
# get recent PRISM data -- plotMonsoonPRISM.R
# auto date range...start with 6-15 and run on 6-17 to get two days of data, end on 10/1
dateRangeStart="2022-06-15"
dateRangeEnd=as.Date(format(as.POSIXct(Sys.time()),usetz=TRUE, tz="Etc/GMT+7"))-1 # date on local time zone
if(dateRangeEnd<"2022-06-16" | dateRangeEnd>="2022-10-01"){
  stop()
}

# generate dates -- keep with PRISM date
allDates<-seq(as.Date(dateRangeStart), as.Date(dateRangeEnd),1)

# AZ/NM bbox -115.004883,31.184609,-102.524414,37.387617
ACISbbox<-"-116,31,-102,38"

# ACIS query
jsonQuery=paste0('{"bbox":"',ACISbbox,'","sdate":"',dateRangeStart,'","edate":"',dateRangeEnd,'","grid":"21","elems":"pcpn","meta":"ll,elev","output":"json"}') # or uid
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
# grab grid for ts extent
#gridStackTS<-gridStack
# set 0 and neg to NA
gridStack[gridStack <= 0] <- NA
# convert to mm
gridStack<-gridStack*25.4
#####

# load SOM mapping from SOM_precip.R
load("~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JJAS_SOM4x4_15K_CP2_1981_2020.RData") #1- rad(4,1):rlen:5000, #2-rad(3,0.33):rlen:7000
mask<-raster( "~/RProjects/SOMs/monsoonPrecip/SWUS_PRISM_MASK.grd")

##### Predict classification of new events - use getDailyPRISMppt.R
#newPrcp<-stack("~/RProjects/SOMs/monsoonPrecip/SWUS_061522_071722_PRISM_daily_prcp.grd")
# resample to same resolution
newPrcp<-resample(gridStack,mask, method="bilinear")
newPrcp<-crop(newPrcp,mask)
newPrcp<-mask(newPrcp,mask)

# crop to region
#e <- extent(-115.5,-106,31.3, 37.5) # "-115,31,-102,38"
#newPrcp <-crop(gridStack,e)
# apply mask
#newPrcp <- mask(newPrcp, mask)
# convert layers to dataframe
new.layers.df<-(as.data.frame(newPrcp, long=TRUE, xy=TRUE))
colnames(new.layers.df)<-c("lon","lat","date","value")  
# long to wide
new.df.wide<-reshape2::dcast(new.layers.df, formula = date~lat+lon, value.var = "value")
new.df.wide[is.na(new.df.wide)] <- 0
# make prediction based on trained SOM
som.prediction <- predict(som.gh500, newdata = as.matrix(new.df.wide[,2:ncol(new.df.wide)]),
                          trainX =as.matrix(df.wide[idx,2:ncol(df.wide)]))

# look at dates and mapped units
predictedUnits<-cbind.data.frame(names(newPrcp),som.prediction$unit.classif)
colnames(predictedUnits)<-c("date","mapUnit")

# get code_grid
codebook<-as.data.frame(som.gh500$codes)
code_grid<-as.data.frame(round(som.gh500$grid$pts,0))
code_grid$mapUnit<-seq(1,nrow(code_grid))
code_grid<-code_grid %>%
  unite(y,x, col="codes", sep="_") # add mapunit back in if needed
# LARGER ACTIVE HIGH CATEGORY
activity<- rbind.data.frame(cbind(c("1_1"), c("Inactive"),c("Inactive")),
                            cbind(c("1_2","2_1","2_2"),rep("Active-low",3),rep("Active",3)),
                            cbind(c("1_3","3_1","3_2","2_3","3_3","1_4","4_1","4_2","2_4"), rep("Active-high",9),rep("Active",9)),
                            cbind(c("3_4","4_3","4_4"), rep("Widespread"),rep("Active")))
colnames(activity)<-c("codes","activityCat","activityCat2")

# merge code_grid and activity
predictedUnits<-merge(predictedUnits,merge(code_grid,activity, by="codes"), by="mapUnit")

# make calendar plot
predictedUnits$date<-as.Date(predictedUnits$date, format="X%Y.%m.%d")      
predictedUnits$wday<-lubridate::wday(predictedUnits$date, label = T, week_start = 7)
predictedUnits$week<-lubridate::epiweek(predictedUnits$date)
predictedUnits$month<-format(predictedUnits$date,"%m")
predictedUnits$label1<-paste0(predictedUnits$codes,"\n",predictedUnits$date)
predictedUnits$label2<-paste0(predictedUnits$codes,"-",predictedUnits$activityCat,"\n",format(predictedUnits$date, "%m-%d"))
# mean regional precip stats
predictedUnits <- predictedUnits[order(as.Date(predictedUnits$date, format="%Y-%m-%d")),]

# plot calendar
predictedUnits %>%
  ggplot(aes(wday,-week, fill = activityCat)) +
  geom_tile(colour = "white")  + 
  geom_text(aes(label = label1), size = 3) +
  theme(aspect.ratio = 1/2,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_manual(values = c("lightblue", "palegreen", "yellow","violet"))+
  facet_wrap(~month, nrow = 4, ncol = 1, scales = "free") +
  labs(title = "SOM precip classifications - Monsoon 2022")

# stats
#predictedUnits$meanPrecip<-apply(new.df.wide[,2:ncol(new.df.wide)], 1, mean, na.rm=TRUE)# mean regional precip
table(predictedUnits$activityCat)
table(predictedUnits$activityCat2)  
#table(somTime$activityCat)/43


# thumbnail map
# PLOT ALL DAYS THUMBNAILS ----
all_states <- map_data("state")
#extent(newPrcp)
xlim = c(-115,-106)
ylim = c(31,38)
# state boundaries
colnames(all_states)<-c("X","Y","PID","POS","region","subregion")
states= PBSmapping::clipPolys(all_states, xlim=xlim,ylim=ylim, keepExtra=TRUE)

# colorramp for total precip
precipCols<-colorRampPalette(c("lightblue", "dodgerblue3", "palegreen","green4","salmon","orangered3",
                               "lightgoldenrod1","orange2","plum2","purple"))(50)
precBreaks<-seq(0,6,0.5)
precLabs<-as.character(seq(0,6,0.5))
precLabs[13]<-">6"
precLabs[1]<-"0.01"
#precBreaksmin<-seq(1,19,2)

#theme_set(theme_bw())
library(rasterVis)
newPrcpIn<-newPrcp/25.4
names(newPrcpIn)<-predictedUnits$label2
p<-gplot(newPrcpIn) + geom_tile(aes(fill = value)) +
  #scale_fill_gradient2(low = 'white', high = 'blue') +
  #scale_fill_distiller(palette = "Spectral", direction = -1, na.value="burlywood4", 
  #                     name="inches", limits=c(0,20),oob=squish)+
  facet_wrap(~ variable) +
  #sugrrants::facet_calendar(~ variable, nrow = 2)+
  scale_fill_gradientn(colours = precipCols, na.value="burlywood4", 
                       name="inches", limits=c(0,6),oob=scales::squish, breaks=precBreaks, labels=precLabs, expand=NULL)+
  guides(fill= guide_colorbar(barheight=15,nbin = 500, raster = FALSE))+
  coord_equal(xlim = c(-115,-106), ylim = c(31,38), expand = FALSE)+
  xlab("Longitude") + ylab("Latitude") 

# p + geom_text(
#               data    = predictedUnits,
#               mapping = aes(x = -114, y = 37, label = codes)
#               )


p<-p +  geom_polygon( data=states, aes(x=X, y=Y, group = PID),colour="grey", fill=NA, size=0.25 )+
  #scale_x_continuous(breaks = c(-120,-140))+
  #ggtitle("Daily Total Precip (inches) - PRISM")+
  ggtitle(paste0("Daily Precipitation (in.): ",allDates[1]," to ",allDates[length(allDates)]))+
  labs(caption=paste0("Plot created: ",format(Sys.time(), "%Y-%m-%d"),
                      "\nThe University of Arizona\nhttps://cals.arizona.edu/climate/\nData Source: PRISM Climate Group\nRCC-ACIS"))+
  theme(plot.title=element_text(size=14, face = "bold"))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# write out file
library(magick)
png("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/Monsoon2022_SOMs.png",width = 16, height = 10, units = "in", res = 300L)
#grid.newpage()
print(p, newpage = FALSE)
dev.off()

# add logos
# Call back the plot
plot <- image_read("/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/Monsoon2022_SOMs.png")
# And bring in a logo
#logo_raw <- image_read("./logos/UA_CLIMAS_logos.png")
logo_raw <- image_read("/home/crimmins/RProjects/ClimPlot/logos/UA_CSAP_CLIMAS_logos_horiz.png") 
logo <- image_resize(logo_raw, geometry_size_percent(width=120,height = 120))
# Stack them on top of each other
#final_plot <- image_append((c(plot, logo)), stack = TRUE)
#final_plot <- image_mosaic((c(plot, logo)))
final_plot <- image_composite(plot, logo, offset = "+410+2760")
# And overwrite the plot without a logo
image_write(final_plot,"/home/crimmins/RProjects/SOMs/monsoonPrecip/figs/Monsoon2022_SOMs.png")
# END PLOT ALL DAYS THUMBNAILS ----







