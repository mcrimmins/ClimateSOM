# make map of 500mb days
# code adapted from rncep_som_test.R
# MAC 01/15/19

library(RNCEP)

# # # download data
# lat.southnorth=c(20,50), lon.westeast=c(210,265) # for CoL analysis
# lat.southnorth=c(20,50), lon.westeast=c(235,265) for NAMS
 wx.extent1 <- NCEP.gather(variable='hgt', level=500,
                           months.minmax=c(1,1), years.minmax=c(2019,2019),
                           lat.southnorth=c(20,50), lon.westeast=c(210,265),
                           reanalysis2 = TRUE, return.units = TRUE, status.bar = FALSE)