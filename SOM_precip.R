# SOM for gridded precipitation 
# adapted from SOM_NARR.R
# MAC 6/26/20

library(raster)
library(lubridate)
library(reshape2)
library(kohonen)
library(tidyr)
library(ggplot2)
library(PBSmapping)
library(dplyr)
library(cowplot)
library(scales)

ptm <- proc.time()

# set rasteroptions
rasterOptions(progress = 'text')

# functions
leap_every_year <- function(x) {
  ifelse(yday(x) > 59 & leap_year(x) == FALSE, yday(x) + 1, yday(x))
}

# map layers
states <- getData('GADM', country='United States', level=1)
  az<-subset(states, NAME_1=="Arizona")
  nm<-subset(states, NAME_1=="New Mexico")
  aznm<-subset(states, NAME_1=="Arizona" | NAME_1=="New Mexico")
us <- getData('GADM', country='United States', level=0)
mx <- getData('GADM', country='Mexico', level=0)
#cn <- getData('GADM', country='Canada', level=0)

# load PRISM for SW Region - from dailyDownloadPRISM.R
 prcp<- stack("/scratch/crimmins/PRISM/processed/SWUS_1981_2019_PRISM_daily_prcp.grd") 
# PRISM percentiles from monsoonPRISM_misc.R
# prcp<-stack("/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2019_PRISM_daily_prcp.grd") 
# load PRISM mask from monsoonPRISM_misc.R
mask<-raster( "~/RProjects/SOMs/monsoonPrecip/SWUS_PRISM_MASK.grd")


# load CPC daily precip from ~/SWMonsoonTracker/NCEPGrids/createNCEPgrids.R
#prcp<- rotate(stack("/scratch/crimmins/cpc_global/processed/CPC_Daily_precip_global_1979_2019_NAClip.grd")) 
#    names(prcp)<-seq.Date(as.Date("1979-01-01"),as.Date("2019-12-31"),1)
    
# # load pre-processed raster stacks from /ClimPlot/NARR/processNARR.R
# gh500<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2019.grd")
# pwat<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2019.grd") 
# prcp<-stack("/scratch/crimmins/NARR/processed/PRCP_daily_NARR_WUS_1979_2019.grd") 
# avgCAPE<-stack("/scratch/crimmins/NARR/processed/CAPE_daily_NARR_WUS_1979_2019.grd")
# maxCAPE<-stack("/scratch/crimmins/NARR/processed/maxCAPE_daily_NARR_WUS_1979_2019.grd")
# sfcDP<-stack("/scratch/crimmins/NARR/processed/DPT2m_daily_NARR_WUS_1979_2019.grd")
# uFlux<-stack("/scratch/crimmins/NARR/processed/WVUFLX_daily_NARR_WUS_1979_2019.grd")
# vFlux<-stack("/scratch/crimmins/NARR/processed/WVVFLX_daily_NARR_WUS_1979_2019.grd")

# crop to SW/NMX region
#e <- extent(-115,-102,25.5, 37.5)
# AZ and western NM
e <- extent(-115.5,-106,31.3, 37.5)
# AZ only
#e <- extent(subset(states, NAME_1=="Arizona"))
#e<-extent(aznm)
#prcp <- (crop(prcp, e)) # can also add in rotate to convert lon to neg

# dates - find and remove leap days
startYr<-1981
dates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
  colnames(dates)<-"date"
  dates$month<-as.numeric(format(dates$date, "%m"))
  dates$day<-as.numeric(format(dates$date, "%d"))
  dates$year<-as.numeric(format(dates$date, "%Y"))
  dates$doy<-as.numeric(format(dates$date, "%j"))
  dates$doy_ly<-leap_every_year(dates$date) # this day of year without leap day shifts
  
# subset layers to months of interest
  #mos<-c(6,7,8,9)
  mos<-c(7,8,9)
  subDates<-dates[which(dates$month %in% mos),]
  #subLayers<-prcp # FOR PERC GRIDS
  subLayers<-prcp[[which(dates$month %in% mos)]]
    # crop to region
    subLayers<-crop(subLayers,e)
    # apply mask
    subLayers <- mask(subLayers, mask)
  
  #subLayers<-prcp # for already subsampled percentiles
  
  # convert layers to dataframe
#  layers.df<-(as.data.frame(subLayers, long=TRUE, xy=TRUE))
#  colnames(layers.df)<-c("lon","lat","date","value")  
  # long to wide
#  df.wide<-dcast(layers.df, formula = date~lat+lon, value.var = "value")
  
  # save raw precip
  #save(df.wide, file="~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS.RData")
  load("~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS.RData")
  # save perc precip
  # save(df.wide, file="~/RProjects/SOMs/monsoonPrecip/perc_AZwNM_PRISM_JAS.RData")
  #  load("~/RProjects/SOMs/monsoonPrecip/perc_AZwNM_PRISM_JAS.RData")
 
  # replace NAs with 0's if needed
  df.wide[is.na(df.wide)] <- 0
  
  # run SOM with all days or subset based on specific node
  idx<-seq(1,nrow(df.wide),1)
  #idx<-which(som.gh500$unit.classif==which(codeList=="1_1"))
  #whichNodes<-c("1_1","1_2","1_3","2_1","2_2","2_3","3_1")
 # idx<-which(som.gh500$unit.classif %in% match(codeList,whichNodes))

  
#  SOM screening #####
# source("topo.error.R")
# ncols=c(2,3,4,3,5,4,6,7,5,4,8,9,6,5,10,11,6,8,5,13,7,14)
# nrows=c(2,2,2,3,2,3,2,2,3,4,2,2,3,4,2,2,4,3,5,2,4,2)
# qe<-c()
# te<-c()
# ptm <- proc.time()
# for(i in 1:length(nrows)){
#   som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]),
#                         #grid = somgrid(ncols[i], nrows[i], "rectangular"),
#                         grid = somgrid(ncols[i], nrows[i], topo="rectangular", neighbourhood.fct = c("gaussian")),
#                         #alpha = c(0.05, 0.001),
#                         #radius = c(5,1),
#                         mode= "pbatch",
#                         cores = 7,
#                         rlen = 500,
#                         dist.fcts = "sumofsquares")
#   ## quantization error:
#   qe[i]<-mean(som.gh500$distances)
#   ## topographical error measures:
#   te[i]<-topo.error(som.gh500, "nodedist")
#   print(paste0(nrows[i],"x",ncols[i]))
# }
# proc.time() - ptm
# 
# diagnostics<-cbind.data.frame(nrows,ncols,qe,te)
# plot(diagnostics$qe, diagnostics$te, xlim=c())
# text(diagnostics$qe, diagnostics$te, labels=paste0(diagnostics$nrows,"-",diagnostics$ncols))
# save(diagnostics, file = "~/RProjects/SOMs/monsoonPrecip/diagnostics.RData")
######

# kohonen SOM
  nrows=4
  ncols=5
  ptm <- proc.time()  
    #som.gh500 <- som(as.matrix(df.wide[,2:ncol(df.wide)]), grid = somgrid(ncols, nrows, "rectangular"))
    set.seed(999) #keep set seed 16, 11 upper left/wet, 9,6 upper right dry/UL wet, 8 UL Wet/LR dry, 7 LR Wet/LL Dry, 5 LR wet/LL dry, 4/100 UL wet/UR dry, 101 wet UR/dry LL, 102/104 dry UL/wet UR, 103 UL Wet/LR Dry  
    som.gh500 <- supersom(as.matrix(df.wide[idx,2:ncol(df.wide)]), 
                          #grid = somgrid(ncols, nrows, "rectangular"),
                          grid = somgrid(ncols, nrows, topo="rectangular", neighbourhood.fct = c("gaussian")),
                          #alpha = c(0.05, 0.001), # for online
                          radius = c(3,1),
                          mode= "pbatch", # "pbatch" or "online"
                          #mode= "online", # "pbatch" or "online"
                          #maxNA.fraction = 0.5,
                          cores = 7,
                          rlen = 500, 
                          dist.fcts = "sumofsquares")
  proc.time() - ptm
  ## quantization error:
  mean(som.gh500$distances)
  ## topographical error measures:
  source("topo.error.R")
  topo.error(som.gh500, "nodedist")
  
  # save SOM output
  #save(som.gh500, file = "~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS_SOM4x5.RData")
  #load("~/RProjects/SOMs/monsoonPrecip/AZwNM_PRISM_JAS_SOM4x5.RData")
  
  #topo.error(som.gh500, "bmu")
  codebook<-as.data.frame(som.gh500$codes)
  code_grid<-as.data.frame(round(som.gh500$grid$pts,0))
  code_grid$mapUnit<-seq(1,nrow(code_grid))
  code_grid<-code_grid %>%
    unite(y,x, col="codes", sep="_") # add mapunit back in if needed
    
  # deal with code book
  codebook<-cbind(code_grid,codebook)
  codebook.long<-melt(codebook, id.vars = 1)
  # separate out mapunits  ----
  mapunits<-codebook.long[1:(nrows*ncols),]
  codebook.long<-codebook.long[-(1:(nrows*ncols)),]
  # ####
  codebook.long<-separate(codebook.long, variable, convert = TRUE, into = c("lat", "lon"), sep="_.")
  codebook.long$lat<-as.numeric(gsub("X", "", codebook.long$lat))
  codebook.long$lon<--codebook.long$lon
  codebook.long<-separate(codebook.long, codes, convert = FALSE, remove = FALSE, into = c("yRows","xCols"), sep="_") # switch position of yRows/xCols
  
  # assign days to nodes
  #nodes<-kohonen::map(som.gh500)
  #somTime<-as.data.frame(cbind(df.wide$date, nodes$unit.classif, nodes$distances))
  somTime<-as.data.frame(cbind(df.wide$date[idx], som.gh500$unit.classif, som.gh500$distances))
  colnames(somTime)<-c("date","mapUnit","errorDist")
  somTime$date<-as.character(somTime$date)
  somTime$date<-as.Date(somTime$date, format="X%Y.%m.%d")
  somTime$month<-as.numeric(format(somTime$date, format="%m"))
  somTime$year <-as.numeric(format(somTime$date, format="%Y"))
  somTime$day  <-as.numeric(format(somTime$date, format="%d"))
  #somTime<-separate(somTime, as.character(date), convert = TRUE, remove = FALSE, into = c("year","month","day"), sep=".")
  #somTime$date<-as.Date(paste(somTime$year,"-",somTime$day,"-",somTime$month,sep=""), format="%Y-%d-%m")
  somTime$doy<-as.numeric(format(somTime$date, "%j"))
  somTime$mapUnit<-as.integer(as.character(somTime$mapUnit))
  somTime$errorDist<-as.numeric(as.character(somTime$errorDist))
  # join codes to node table
  somTime<-left_join(somTime, code_grid)
  # get codeList
  codeList<-unique(codebook.long$codes)
  
  # expand grid to double
  #som.gh500.2 <- expandMap(som.gh500)
  #plot(som.gh500.2, type="counts", shape = "straight", labels=counts)
  summary(som.gh500)
  
    
  # RASTERVIS mapping of SOM results
  library(rasterVis)
  library(colorspace)
  # create stack of codebook results
  cbRaster<-stack()
  for (i in 1:length(codeList)) {
    temp<-codebook.long[which(codebook.long$codes==codeList[i]),c(5,4,6)]
    temp<-rasterFromXYZ(temp)
    cbRaster<-stack(cbRaster, temp)
  }
  names(cbRaster)<-codeList

  cbRaster[cbRaster < 0.254] <- NA  
  at<-c(seq(0,40,1))
  #mapTheme <- rasterTheme(region=rev(terrain_hcl(12)))
  mapTheme <- rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
  pPrecip<-levelplot(cbRaster, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows), at=at,
            main="Precip Patterns JAS 3x3 SOM - PRISM-daily 1981-2019")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
    #layer(sp.polygons(us, col = 'gray40', lwd=1))+
    #layer(sp.polygons(mx, col = 'gray40', lwd=1))
  png("/home/crimmins/RProjects/SOMs/precip_SOM_cpc.png", width = 10, height = 6, units = "in", res = 300L)
  #grid.newpage()
  print(p, newpage = FALSE)
  dev.off() 
  
  # spatial pattern correlation
  #somTime$kendall<-NA
  somTime$spearman<-NA
  somTime$pearson<-NA
  for(i in 1:nrow(somTime)){
     somTime$pearson[i]<-cor(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]),
                              use = "na.or.complete", method="pearson")
    somTime$spearman[i]<-cor(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]),
        use = "na.or.complete", method="spearman")
    # somTime$kendall[i]<-cor(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]),
    #                          use = "na.or.complete", method="kendall")
    print(i)
  }
  
  mean(somTime$spearman, na.rm=TRUE)
  mean(somTime$pearson, na.rm=TRUE)
    # i<-234
  # plot(values(cbRaster[[somTime$mapUnit[i]]]), values(subLayers[[i]]))
  # plot(stack(cbRaster[[somTime$mapUnit[i]]], subLayers[[i]]))
  #   plot(values(cbRaster[[somTime$mapUnit[i]]]), type="l")
  #   lines(values(subLayers[[i]]), col="red")
  
  # PRECIP METRICS FOR DAYS
  # dist of all values
  # quantile(as.matrix(df.wide[2:ncol(df.wide)]), probs = c(.25, .5, .75))
    
  # metrics of spatial autocorr?
  # spatial patterns of precip, metrics for each day
  somTime$percExtent<-(rowSums(df.wide[idx,2:ncol(df.wide)]>0)/(ncol(df.wide)-1))*100 # percent extent >0
  somTime$maxPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, max) # max value of day
  somTime$meanPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, mean)# mean regional precip
  somTime$medPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, median) # max value of day  
  somTime$percZero<-(rowSums(df.wide[idx,2:ncol(df.wide)]==0)/(ncol(df.wide)-1))*100 # percent 0 precip  
  somTime$sumPrecip<-apply(df.wide[idx,2:ncol(df.wide)], 1, sum) # max value of day 
  
  
  # distribution of daily precip values by node
  ggplot(somTime, aes(as.factor(somTime$codes), percExtent))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Daily Extent Precip (%) by nodes")
  ggplot(somTime, aes(as.factor(somTime$codes), maxPrecip))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Max Daily Precip (mm) by nodes")
  ggplot(somTime, aes(as.factor(somTime$codes), meanPrecip))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Mean Daily Precip (mm) by nodes")  
  ggplot(somTime, aes(as.factor(somTime$codes), medPrecip))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Median Daily Precip (mm) by nodes")  
  ggplot(somTime, aes(as.factor(somTime$codes), percZero))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Percent Zero Extent Precip (mm) by nodes")  
  
  # plot 10 random maps from selected node to assess quality
  #at<-c(seq(0.01,45,1),60)
  mapTheme <- rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
  temp<-subLayers[[idx]]
  levelplot(temp[[sample(which(somTime$codes=="1_1"),20)]], contour=FALSE,
                     margin=FALSE, par.settings=mapTheme, 
                     main="Precip Patterns from 1_1 - PRISM-daily 1981-2019")+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  
  # maps on composites of days in nodes ####
  # percentile data for composites
   perc<-stack("/scratch/crimmins/PRISM/processed/JASperRank_SWUS_1981_2019_PRISM_daily_prcp.grd") 
      perc<-crop(perc,e)
      # apply mask
      temp<-mask(perc, mask)
  # composites on raw precip
  #temp<-subLayers[[idx]]
  comp<-stackApply(temp, somTime$mapUnit, fun=IQR)
    tempNames<-as.data.frame(names(comp))
    colnames(tempNames)<-"index"
    tempNames <- tempNames %>% separate(index, c(NA,"mapUnit"), sep = "_", remove=FALSE)
    comp<-subset(comp, order(as.numeric(tempNames$mapUnit)))
  names(comp)<-codeList   
  # plot 
  at<-c(seq(0,5,0.25),45)
  mapTheme<-rasterTheme(region = c("lightblue", "blue", "green","yellow","red"))
  pComp<-levelplot(comp, contour=FALSE, margin=FALSE, par.settings=mapTheme,layout=c(ncols,nrows),  
                     main="Composite IQR Percentile JAS 4x5 SOM - PRISM-daily 1981-2019", xlab=NULL, ylab=NULL)+
    layer(sp.polygons(aznm, col = 'gray40', lwd=1))
  #####

  
  
  # Clustering of nodes
  ## Show the U matrix
  Umat <- plot(som.gh500, type="dist.neighbours", main = "SOM neighbour distances")
  ## use hierarchical clustering to cluster the codebook vectors
  som.hc <- cutree(hclust(object.distances(som.gh500, "codes")),5)
  add.cluster.boundaries(som.gh500, som.hc)
  
  # map clustered nodes ----
  
  #####
  
  # DIAGNOSTICS
  plot(som.gh500, type="changes") # changes, codes, counts, property, quality, mapping
  
  # summary plots -- appears to plot opposite up/down from SOM plot
  counts <- plot(som.gh500, type="counts", shape = "straight", labels=counts)
  codes <- plot(som.gh500, type="codes", shape = "straight")
  similarities <- plot(som.gh500, type="quality", palette.name = terrain.colors)
  plot(som.gh500, type="dist.neighbours", main = "SOM neighbour distances")
  plot(som.gh500)
  # sammon mapping
  library(MASS)
  gh500.codes <- som.gh500$codes
  dis <- dist(as.matrix(som.gh500$codes[[1]]))
  gh500.sam <- sammon(dis)
  plot(gh500.sam$points, type="n")
  text(gh500.sam$points,labels=as.character(1:nrow(code_grid)))
  ##  Polygon version of map ----
  library(sp)
    temp<-list()
    ctr<-1
    for(j in 1:(nrows-1)){
      for(k in 1:(ncols-1)){
        temp[ctr]<-Polygons(list(Polygon(cbind(c(gh500.sam$points[(k+(ncols*j)-ncols),1], gh500.sam$points[(k+(ncols*j)-ncols)+1,1],
                                           gh500.sam$points[k+(ncols*j)+1,1],gh500.sam$points[k+(ncols*j),1]),
                                         c(gh500.sam$points[(k+(ncols*j)-ncols),2], gh500.sam$points[(k+(ncols*j)-ncols)+1,2],
                                           gh500.sam$points[k+(ncols*j)+1,2],gh500.sam$points[k+(ncols*j),2])))),paste0(ctr))
        ctr<-ctr+1
      }
    }  
    sr1<-SpatialPolygons(temp)  
    plot(gh500.sam$points, type="n")
    text(gh500.sam$points,labels=as.character(1:nrow(code_grid)))
    plot(sr1, add=TRUE)
  # ----
    
    
  # GGPLOT versions
  # counts --- plotting incorrectly
  counts<-somTime %>% group_by(codes) %>% count(codes)
    counts <- counts %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    counts$row<-as.numeric(counts$row); counts$col<-as.numeric(counts$col); 
    counts$avg<-(counts$n/length(seq(1981,2019,1)))
pCt<-ggplot(counts, aes(x=col,y=-row))+
      geom_tile(aes(fill = (n/nrow(somTime)*100)))+
      geom_text(aes(label = codes), size=6)+
      scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="% days")+
      #scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 20,
      #                     high = "red", na.value = NA, name="% days")+
      ggtitle("Percent of days/node")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
  # SOM distances
    ndist <- unit.distances(som.gh500$grid)
    cddist <- as.matrix(object.distances(som.gh500, type = "codes"))
    cddist[abs(ndist - 1) > .001] <- NA
    neigh.dists <- colMeans(cddist, na.rm = TRUE)
    som_grid <- som.gh500[[4]]$pts %>%
      as_tibble %>% 
      mutate(id=row_number())
    som_grid$codes<-paste0(som_grid$y,"_",som_grid$x)
    som_grid <- som_grid %>% mutate(dist=neigh.dists)
pDist<-ggplot(som_grid, aes(x=x,y=-y))+
      geom_tile(aes(fill = dist))+
      geom_text(aes(label = codes), size=6)+
      scale_fill_gradient(low = "lightblue", high = "red", na.value = NA, name="distance")+
      ggtitle("Neighborhood Distance")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
# SOM Spearman corr by day
spear<-somTime %>% group_by(codes) %>% summarise(spearR = mean(spearman, na.rm=TRUE))
  spear <- spear %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  spear$row<-as.numeric(spear$row); spear$col<-as.numeric(spear$col);
  pRho<-ggplot(spear, aes(x=col,y=-row))+
    geom_tile(aes(fill = spearR))+
    geom_text(aes(label = round(spear$spearR,2)), size=6)+
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="mean Spearman")+
    ggtitle("Mean Spearman Corr by node")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())  
# boxplot of spearman values
  ggplot(somTime, aes(as.factor(somTime$codes), spearman))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Spearman Rho by nodes")
  ggplot(somTime, aes(as.factor(somTime$codes), pearson))+
    geom_boxplot(varwidth = TRUE)+
    ggtitle("Distribution of Pearson r by nodes") 
  
  # SOM node quality
    sim<-cbind.data.frame(som.gh500$unit.classif,som.gh500$distances)
      colnames(sim)<-c("node","distance")
    sim <- sim %>% group_by(node) %>% summarise(dist = mean(distance))
    sim$code<-codeList
    sim <- sim %>% separate(code, c("row","col"), sep = "_", remove=FALSE)
      sim$row<-as.numeric(sim$row); sim$col<-as.numeric(sim$col); 
  pQ<-ggplot(sim, aes(x=col,y=-row))+
        geom_tile(aes(fill = dist))+
        geom_text(aes(label = code), size=6)+
        scale_fill_gradient(low = "yellow", high = "red", na.value = NA, name="mean error")+
        ggtitle("SOM Node Quality")+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
  plot_grid(pCt,pDist,pQ,pRho, ncol = 1, align = "v")
  # SOM code vectors
  
  # SOM code vectors distributions
  ggplot(codebook.long, aes(as.factor(codebook.long$codes), value))+
    geom_boxplot(varwidth = TRUE)+
    xlab("Node")+
    ylab("mm")+
    ggtitle("Distribution of Codebook Vectors (precip in mm)")
  
  # SOM node counts by year - facet wrap indiv heat maps
  countsYr<-somTime %>% group_by(codes,year) %>% count(codes)
  countsYr <- countsYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  countsYr$row<-as.numeric(countsYr$row); countsYr$col<-as.numeric(countsYr$col); 
    ggplot(countsYr, aes(x=col, y=-row))+
      geom_tile(aes(fill=n))+
      geom_text(aes(label = n), size=4)+
      facet_wrap(~year)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 25,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Count of days in each node by year")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  # time series
  temp<-subset(countsYr, codes=="4_2")
  ggplot(temp, aes(x=year, y=n, color=as.factor(codes)))+
    geom_line()
  
  # SOM node anomaly
  countAnom<-merge(countsYr, counts, by="codes")
  countAnom$anomCT<-countAnom$n.x/countAnom$avg
  ggplot(countAnom, aes(x=col.x, y=-row.x))+
    geom_tile(aes(fill=anomCT))+
    #geom_text(aes(label = anomCT), size=4)+
    facet_wrap(~year)+
    scale_fill_gradient2(low = "purple",mid="white", midpoint = 1,
                         high = "orange", na.value = NA, name="% of avg",limits=c(0, 2), oob=squish)+
    ggtitle("Anom of days in each node by year")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  
  
  # SOM node counts by month
  countsMo<-somTime %>% group_by(codes,month) %>% count(codes)
  countsMo <- countsMo %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  countsMo$row<-as.numeric(countsMo$row); countsMo$col<-as.numeric(countsMo$col); 
    ggplot(countsMo, aes(x=col, y=-row))+
      geom_tile(aes(fill=n))+
      geom_text(aes(label = n), size=4)+
      facet_wrap(~month)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 350,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Count of days in each node by month")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
  
  # classification error by node-year
  errorYr<-somTime %>% 
           group_by(codes,year) %>%
           summarise(error=mean(errorDist))
  errorYr <- errorYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  errorYr$row<-as.numeric(errorYr$row); errorYr$col<-as.numeric(errorYr$col); 
  ggplot(errorYr, aes(x=col, y=-row))+
    geom_tile(aes(fill=round(error,0)))+
    #geom_text(aes(label = error), size=2)+
    facet_wrap(~year)+
    scale_fill_gradient2(low = "blue",mid="yellow", 
                         high = "red", midpoint = 300000, 
                         na.value = NA, name="class. error")+
    ggtitle("Classification error by node-year")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
      # error by year
      errorYr<-somTime %>% 
        group_by(year) %>%
        summarise(error=median(errorDist))
      ggplot(errorYr, aes(x=year, y=error))+
        geom_bar(stat = "identity", fill="red")+
        ggtitle("Median error of all nodes by year")

  # classification error by node-year
    spearYr<-somTime %>% 
        group_by(codes,year) %>%
        summarise(error=mean(spearman))
      spearYr <- spearYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
      spearYr$row<-as.numeric(spearYr$row); spearYr$col<-as.numeric(spearYr$col); 
      ggplot(spearYr, aes(x=col, y=-row))+
        geom_tile(aes(fill=error))+
        #geom_text(aes(label = error), size=2)+
        facet_wrap(~year)+
        scale_fill_gradient2(low = "blue",mid="yellow", 
                             high = "red", midpoint = 0.5, 
                             na.value = NA, name="mean rho")+
        ggtitle("Mean Spearman Corr by node-year")+
        theme_bw()+
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank())
      # error by year
      spearYr<-somTime %>% 
        group_by(year) %>%
        summarise(error=mean(spearman, na.rm=TRUE))
      ggplot(spearYr, aes(x=year, y=error))+
        geom_bar(stat = "identity", fill="red")+
        ggtitle("Mean Spearman Rho of all nodes by year")     
      
        
  # % of seasonal precip by node
  sumYr<-somTime %>% 
    group_by(codes,year) %>%
    summarise(sumPrecip=sum(meanPrecip), countNodes=n())
  sumYr <- sumYr %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
  sumYr$row<-as.numeric(sumYr$row); sumYr$col<-as.numeric(sumYr$col); 
  ggplot(sumYr, aes(x=col, y=-row))+
    geom_tile(aes(fill=round(sumPrecip,0)))+
    #geom_text(aes(label = error), size=2)+
    facet_wrap(~year)+
    scale_fill_gradient2(low = "blue",mid="yellow", 
                         high = "red", midpoint = 20, 
                         na.value = NA, name="Sum Precip")+
    ggtitle("Sum of Precip by node-year")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank())
  
  # more daily/seasonal precip metrics
  sumYr2<-somTime %>% 
    group_by(year) %>%
    summarise(sumPrecip=sum(meanPrecip))
  sumYr2$ltAvg<-mean(sumYr2$sumPrecip)
  sumYr2$seasAnom<-sumYr2$sumPrecip-sumYr2$ltAvg
  sumYr2<- merge(sumYr, sumYr2, by='year')
  sumYr2$anomName<-"normal"
    sumYr2$anomName[sumYr2$seasAnom<(-10)] <- "dry"
    sumYr2$anomName[sumYr2$seasAnom>(10)] <- "wet"
  sumYr2 <- sumYr2 %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    sumYr2$row<-as.numeric(sumYr2$row); sumYr2$col<-as.numeric(sumYr2$col); 
  ggplot(sumYr2, aes(x=col, y=-row))+
      geom_tile(aes(fill=countNodes))+
      #geom_text(aes(label = countNodes), size=2)+
      facet_grid(as.factor(anomName)~year)+
      #facet_wrap(as.factor(anomName)~year)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", 
                           high = "red", midpoint = 25, 
                           na.value = NA, name="Count")+
      ggtitle("Node counts/year - Wet/Normal/Dry")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
  # counts/anom
    anomCount<-sumYr2 %>% group_by(codes, anomName) %>% summarise(counts=sum(countNodes))
    anomCount <- anomCount %>% separate(codes, c("row","col"), sep = "_", remove=FALSE)
    anomCount$row<-as.numeric(anomCount$row); anomCount$col<-as.numeric(anomCount$col); 
    ggplot(anomCount, aes(x=col, y=-row))+
      geom_tile(aes(fill=counts))+
      geom_text(aes(label = counts), size=4)+
      facet_wrap(~anomName)+
      scale_fill_gradient2(low = "lightblue",mid="yellow", midpoint = 400,
                           high = "red", na.value = NA, name="count")+
      ggtitle("Count of days in each node by Anomaly")+
      theme_bw()+
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank())
    
  # trends in nodes plot; plot with precip 
  
  
  # plot map units
  ggplot(somTime, aes(doy, year)) + 
    geom_tile(aes(fill = mapUnit), colour = "grey") + 
    geom_text(aes(label = codes), size=1.50)+
    scale_fill_gradient2(low = "lightblue", mid = "green",
                         high = "orange", midpoint = (ncols*nrows)/2, space = "Lab",
                         na.value = "grey50", guide = "colourbar")+
    ggtitle("Daily Precip Type Classification")
  # plot error
  ggplot(somTime, aes(doy, year)) + 
    geom_tile(aes(fill = errorDist), colour = "grey") + 
    scale_fill_gradient2(low = "white", mid = "yellow",
                         high = "red", space = "Lab",
                         na.value = "grey50", guide = "colourbar")+
    ggtitle("Daily Precip Type Classification Error")
  # plot daily spearman
  ggplot(somTime, aes(doy, year)) + 
    geom_tile(aes(fill = spearman), colour = "grey") + 
    scale_fill_gradient2(low = "lightblue", mid = "yellow",
                         high = "red", space = "Lab",
                         na.value = "grey50", guide = "colourbar")+
    ggtitle("Daily Map Correlation (Rho)")
  
  # create composites on nodes
  NARR<-stack("/scratch/crimmins/NARR/processed/GH500_daily_NARR_WUS_1979_2019.grd")
  #NARR<-stack("/scratch/crimmins/NARR/processed/PWAT_daily_NARR_WUS_1979_2019.grd") 
  
  # dates - find and remove leap days
  startYr<-1979
  compDates<-as.data.frame(seq.Date(as.Date(paste0(startYr,"-01-01")),as.Date("2019-12-31"),1))
    colnames(compDates)<-"date"
  compLayers<-NARR[[match(somTime$date,compDates$date)]]
  compMean<-stackApply(compLayers, somTime$mapUnit, fun=mean)
    names(compMean)<-codeList   
  compSD<-stackApply(compLayers, somTime$mapUnit, fun=sd)
    names(compSD)<-codeList     
  # plot 
    at<-seq(0,60,5)
    mapTheme<-rasterTheme(region=brewer.pal(8,"BrBG")) 
  pComp<-levelplot(compMean, contour=FALSE, margin=FALSE,layout=c(ncols,nrows), 
                 main="PWAT Composite Patterns JAS 2x3 SOM - NARR 1981-2019", par.settings=mapTheme)+
      #contourplot(compSD,linetype = "dashed")+
      layer(sp.polygons(az, col = 'gray40', lwd=1))
  
    library(gridExtra)
    grid.arrange(pPrecip, pComp, ncol=1)
    
  # ggplot to plot rasters of diff resolutions
  # First, to a SpatialPointsDataFrame
  compMean_pts <- rasterToPoints(compMean, spatial = TRUE)
  cbRaster_pts <- rasterToPoints(cbRaster, spatial = TRUE)
  # Then to a 'conventional' dataframe
  compMean_df  <- data.frame(compMean_pts)
    compMean_df<- melt(compMean_df, id.vars=c("x", "y"))
    compMean_df<-subset(compMean_df, variable!="optional")
  cbRaster_df  <- data.frame(cbRaster_pts)
    cbRaster_df<- melt(cbRaster_df, id.vars=c("x", "y"))
    cbRaster_df<-subset(cbRaster_df, variable!="optional")
    rm(compMean_pts)
    rm(cbRaster_pts)
    
  # fortify polygons
    #usPoly <- fortify(us, region="NAME_0")
    library(maps)
    usPoly <- map_data("state")
    mxPoly <- map_data(database = "world", regions = "Mexico")
    
  ggplot() +
        geom_raster(data = cbRaster_df , aes(x = x, y = y, fill = value))+
        geom_polygon(data=usPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
        geom_polygon(data=mxPoly,aes(long, lat, group = group),
                 fill = NA, col = "black", size = 0.2)+
        stat_contour(data = compMean_df , aes(x = x, y = y, z = value, color = ..level..))+
        scale_fill_gradientn(colors = c("lightblue", "blue", "green","yellow","red"))+
        #scale_color_continuous(type = "viridis") +
        scale_color_gradient2(midpoint = 5800, low = "blue", mid = "yellow", high = "red", name = "GH" )+
        #scale_color_gradient2(midpoint = 25, low = "brown", mid = "lightgrey", high = "green", name = "PWAT(mm)" )+
        coord_cartesian(xlim = c(-125, -100), ylim = c(25, 50))+
        facet_wrap(~variable)
  

##### Predict classification of new events
  newPrcp<-stack("~/RProjects/SOMs/monsoonPrecip/SWUS_070120_071520_PRISM_daily_prcp.grd")
  # crop to region
  newPrcp <-crop(newPrcp,e)
  # apply mask
  newPrcp <- mask(newPrcp, mask)
  # convert layers to dataframe
    new.layers.df<-(as.data.frame(newPrcp, long=TRUE, xy=TRUE))
    colnames(new.layers.df)<-c("lon","lat","date","value")  
  # long to wide
    new.df.wide<-dcast(new.layers.df, formula = date~lat+lon, value.var = "value")
    new.df.wide[is.na(new.df.wide)] <- 0
  # make prediction based on trained SOM
    som.prediction <- predict(som.gh500, newdata = as.matrix(new.df.wide[,2:ncol(new.df.wide)]),
                              trainX =as.matrix(df.wide[idx,2:ncol(df.wide)]))
    
  
  