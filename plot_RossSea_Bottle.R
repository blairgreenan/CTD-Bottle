# plot_RossSea_Bottle_Data.R
# Blair Greenan
# Fisheries and Oceans Canada
# 1 May 2023
#
# load libraries
library(oce)
library(R.matlab)
library(R.oo)
library(tidyverse)
library(lubridate)
library(cmocean)

# load CTD Data from the .mat file
NBP1201_bottle <- readMat("nbp1201_bottle_data_Apr_2017.mat")

# load the water depth at CTD station locations
# this information was compiled from the Cruise Event Log spreadsheet
Water_depth <- read.delim("CTD_STN_WATER_DEPTH.txt", header = FALSE, sep = "\n")

# create a data frame with the CTD data
data <- data.frame(NBP1201_ctd$ctd.data)

# extract the names of the columns in the CTD data and use the trimws
# function to trim the white space at the start and end of the names
# Note: there are a few cases below where symbols were included in the column name
#       and the substr function would not work on those names so these few cases have been renamed.
names(data)[1:9] <- trimws(substr(NBP1201_bottle$CTD.casts.info.columns, 5, 36))
NBP1201_bottle$CTD.BTL.data.columns[3] <- "12 - Density [sigma-theta, Kg/m^3"
NBP1201_bottle$CTD.BTL.data.columns[4] <- "13 - Density2 [sigma-theta, Kg/m^3"
names(data)[10:29] <- trimws(substr(NBP1201_bottle$CTD.BTL.data.columns, 6, 60))
names(data)[30:35] <- trimws(substr(NBP1201_bottle$nutrient.data.columns, 6, 60))
names(data)[36:40] <- trimws(substr(NBP1201_bottle$chlorophyll.data.columns, 6, 60))
names(data)[41] <- trimws(substr(NBP1201_bottle$TM.dissolved.iron.data.columns, 6, 60))
names(data)[42:48] <- trimws(substr(NBP1201_bottle$TM.casts.info.columns, 6, 60))
names(data)[49:51] <- trimws(substr(NBP1201_bottle$particulate.data.columns, 6, 60))
NBP1201_bottle$TM.metals.0.4.mu.columns[1] <- "52 - 0.4 micron filter code"
NBP1201_bottle$TM.metals.0.4.mu.columns[2] <- "53 - 0.4 micron filter number"
names(data)[52:92] <- trimws(substr(NBP1201_bottle$TM.metals.0.4.mu.columns, 6, 60))
names(data)[93:95] <- trimws(substr(NBP1201_bottle$primary.production.columns, 6, 60))
NBP1201_bottle$TM.metals.2mu.columns[1] <- "96 - 2 micron filter"
NBP1201_bottle$TM.metals.2mu.columns[2] <- "96 - 2 micron filter ID"
names(data)[96:136] <- trimws(substr(NBP1201_bottle$TM.metals.2mu.columns, 6, 60))
names(data)[137:177] <- trimws(substr(NBP1201_bottle$HPLC.data.columns, 7, 60))


# CTD metadata (lat, lon, etc.) from the mat file
stn_list <- data.frame(NBP1201_ctd$ctd.list)

# extract the names of the columns in the CTD metadata and use the trimws
# function to trim the white space at the start and end of the names
names(stn_list) <- trimws(substr(NBP1201_ctd$list.variables, 5, 19))

# create a set of unique cast numbers (note that there is a cast 4.2,
# so that means there are 118 locations with 119 CTD casts)
cast <- unique(data[['CTD cast']])

# Edit the longitude so that it matches the convention in the Ross Sea bathymetry data
for (i in 1:length(stn_list[['longitude (deg)']])){
  if(stn_list[['longitude (deg)']][i]>0)
  {stn_list[['longitude (deg)']][i]<-stn_list[['longitude (deg)']][i]-360}
}

# create an empty list for the ctd object
ctd <- list()

# loop to populate the ctd object with data and metadata for each CTD cast
for (i in seq_along(cast)) {
  cat('Cast', cast[i], '...')
  II <- data[['CTD cast']] == cast[i]
  ## create the ctd object
  tmp <- as.ctd(salinity=data$Sal[II],
                temperature=data$Temp[II],
                pressure=data$Pres[II],
                station = data$`CTD cast`[II][1],
                cruise = "NBP1201",
                ship = "R/V Nathaniel B. Palmer") # just use the first element of the station vector in the data - oddity of the MAT file having a station # for each point in the CTD profile
  # add the other fields that were collected by the CTD system on the R/V NBP
  # I am going to ignore the Time vector in the data since it is not particularly useful
  fields <- names(data)
  fields <- fields[-which(fields %in% c('Station', 'Sal', 'Temp', 'Pres', 'Time'))]
  for (f in fields) {
    tmp <- oceSetData(tmp, f, data[[f]][II])
  }
  # add metadata latitude, longitude and start time of the CTD cast to the ctd object
  tmp <- oceSetMetadata(tmp, 'longitude', stn_list[['longitude (deg)']][i])
  tmp <- oceSetMetadata(tmp, 'latitude', stn_list[['latitude (deg)']][i])
  tmp <- oceSetMetadata(tmp, 'waterDepth', Water_depth[i,1])
  tmp <- oceSetMetadata(tmp, 'startTime',
                        as.POSIXct(
                          paste0(stn_list$year[i], '-', stn_list$month[i], '-', stn_list$day[i], ' ', stn_list$hour[i],
                                 ':', stn_list$minute[i], ':', stn_list$second[i]),
                          tz='UTC'))
  ctd[[i]] <- tmp
  cat('\n')
}

# Create a section using the function from the oce package
ctd_section <- as.section(ctd)
# Create a subset for the 2nd occupation of NE-SW section over Ross Bank CTD casts 74-80
RB <- subset(ctd_section, 74 <= stationId & stationId <= 80)
# Grid the data in the subset using the oce SectionGrid function
RBgrid <- sectionGrid(RB, p=seq(0,1000,10))
# Reverse the order of the stations so that it is presented in descending order which presents better as West on left and East on right side of plot
RBgrid2 <- sectionSort(RBgrid, decreasing = TRUE)
# Add topography to the plot
RossSeaBathy <- read.topo(file = "C:/Users/greenanb/Documents/Science Projects/Current/Ross Sea/Data/Ross Sea Bathymetry/topo_198W_174W_78.5S_72.5S_1min.nc")
#dev.new()
#dev.new()
# Print figure to a TIFF file format
tiff("CTD74-80.tiff", width=6, height=6, units='in', res=1200, compression = 'lzw')
par(mfrow=c(3,2))
plot(RBgrid2, which="temperature", ztype = "image", zcol = cmocean('thermal'), zbreaks=seq(-2.5, 0.5, 0.1), showBottom = RossSeaBathy, legend.text = 'A', xlab="", ylim = c(0, 600))
text(5,550,expression("Temperature (\u00B0C)"), adj=0)
plot(RBgrid2, which="Fluor", ztype = "image", zcol = cmocean('algae'), zbreaks=seq(0, 10, 0.1), showBottom = RossSeaBathy, legend.text = 'B', xlab="", ylim = c(0, 600))
text(5,550,expression("Fluorescence (\u03BCg/L)"), adj=0)
plot(RBgrid2, which="salinity", ztype = "image", zcol = cmocean('haline'), zbreaks=seq(34.4, 34.6, 0.01), showBottom = RossSeaBathy, legend.text = 'C', xlab="", ylim = c(0, 600))
text(5,550,"Salinity", adj=0)
plot(RBgrid2, which="Oxy_s", ztype = "image", zcol = cmocean('oxy'), zbreaks=seq(6, 8, 0.1), showBottom = RossSeaBathy, legend.text = 'D', xlab="", ylim = c(0, 600))
text(5,550,"Oxygen (ml/L)", adj=0)
plot(RBgrid2, which="Sig", ztype = "image", zcol = cmocean('dense'), zbreaks=seq(27.6, 27.9, 0.01), showBottom = RossSeaBathy, legend.text = 'E', ylim = c(0, 600))
text(5,550,"Sigmat", adj=0)
plot(RBgrid2, which="Trans", ztype = "image", zcol = cmocean('matter'), zbreaks=seq(90, 100, 0.5), showBottom = RossSeaBathy, legend.text = 'F', ylim = c(0, 600))
text(5,550,"Transmission (%)", adj=0)
# Close of the TIFF image to force the write to file
dev.off()



