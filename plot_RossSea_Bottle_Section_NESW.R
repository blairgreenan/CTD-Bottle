# plot_RossSea_Bottle_Data.R
# Blair Greenan
# Fisheries and Oceans Canada
# 1 May 2023
#
# Description: this script generates a faceted plot with profile of CTD bottle
# that were carried out over the NE-SW section across Survey 1 of Ross Bank.
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

# create a data frame with the CTD data
data <- data.frame(NBP1201_bottle$data)

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


# create a set of unique cast numbers (note that there is a cast 4.2,
# so that means there are 118 locations with 119 CTD casts)
STN <- unique(data[['Station #']])
# create a set of unique TMCTD cast numbers
CTD_cast <- unique(data[['CTD Cast #']])
# Get rid of NAN from the cast # list
CTD_cast <- CTD_cast[-which(is.nan(CTD_cast))]

# Edit the longitude so that it matches the convention in the Ross Sea bathymetry data
for (j in 1:length(data$`Longitude (decimal Deg)`)){
  # if the data is NAN, then skip over this
  if(!is.nan(data$`Longitude (decimal Deg)`[j])){
    if(data$`Longitude (decimal Deg)`[j]>0)
    {data$`Longitude (decimal Deg)`[j]<-data$`Longitude (decimal Deg)`[j]-360}
  }
}


# There seems to be an issue with NA in the cast list that breaks the procedure
# that I used to for the CTD sections, but is not working for the bottle data
# NA_CTD_cast <- !is.na(data$`CTD Cast #`)

# create an empty list for the ctd_bottle object
ctd_bottle <- list()

# loop to populate the ctd object with data and metadata for each CTD cast
for (i in seq_along(CTD_cast)) {
  cat('CTD cast', CTD_cast[i], '...')
#  II <- data[['CTD Cast #']][NA_CTD_cast] == CTD_cast[i]
  II <- which(data[['CTD Cast #']] == CTD_cast[i])
  ## create the ctd_bottle object
  tmp <- as.ctd(salinity=data$`sal00: Salinity, Practical [PSU]`[II],
                temperature = data$`t090C: Temperature [ITS-90, deg C]`[II],
                pressure = data$`prDM: Pressure, Digiquartz [db]`[II],
                station = data$`CTD Cast #`[II][1],
                cruise = "NBP1201",
                ship = "R/V Nathaniel B. Palmer") # just use the first element of the station vector in the data - oddity of the MAT file having a station # for each point in the CTD profile
  # add the other fields that were collected by the CTD system on the R/V NBP
  # I am going to ignore the Time vector in the data since it is not particularly useful
  fields <- names(data)
  fields <- fields[-which(fields %in% c('CTD Cast #', 'sal00: Salinity, Practical [PSU]', 't090C: Temperature [ITS-90, deg C]', 'prDM: Pressure, Digiquartz [db]'))]
  for (f in fields) {
    tmp <- oceSetData(tmp, f, data[[f]][II])
  }
  # add metadata latitude, longitude and start time of the CTD cast to the ctd object
  tmp <- oceSetMetadata(tmp, 'longitude', data$`Longitude (decimal Deg)`[II][1])
  tmp <- oceSetMetadata(tmp, 'latitude', data$`Latitude (decimal Deg)`[II][1])
  tmp <- oceSetMetadata(tmp, 'startTime', ymd(data$`Date  yyyymmdd`[II][1]) + hms(paste(floor(data$`UTC Time hhmm`[II][1]/100), data$`UTC Time hhmm`[II][1]%%100, "00", sep=":")))
  ctd_bottle[[i]] <- tmp
  cat('\n')
  
}

# Create a section using the function from the oce package
ctd_section <- as.section(ctd_bottle)
# Create a subset for the NE-SW section over Ross Bank CTD casts 37-44
# Note that cast 42 is not included because this was a calibration cast at the SeaHorse mooring site and not a part of the CTD secion
RB <- subset(ctd_section, 37 <= stationId & stationId <= 44 & stationId != 42)
# Grid the data in the subset using the oce SectionGrid function
RBgrid <- sectionGrid(RB, p=seq(0,1000,10))
# Reverse the order of the stations so that it is presented in descending order which presents better as West on left and East on right side of plot
RBgrid2 <- sectionSort(RBgrid, decreasing = TRUE)
# Add topography to the plot
RossSeaBathy <- read.topo(file = "C:/Users/greenanb/Documents/Science Projects/Current/Ross Sea/Data/Ross Sea Bathymetry/topo_198W_174W_78.5S_72.5S_1min.nc")
RossSeaBathy@data$z <- -1*RossSeaBathy@data$z
#dev.new()
#dev.new()
# Print figure to a TIFF file format
tiff("CTD37-44_bottle.tiff", width=6, height=6, units='in', res=1200, compression = 'lzw')
par(mfrow=c(3,2))
plot(RBgrid2, which="Phosphate (uM)", ztype = "image", zcol = cmocean('matter'), zbreaks=seq(1, 2.1, 0.1), showBottom = RossSeaBathy, legend.text = 'A', xlab="", ylim = c(600, 0))
text(5,550,expression("Phosphate (\u03BCM)"), adj=0)
plot(RBgrid2, which="Ammonium(uM)", ztype = "image", zcol = cmocean('matter'), zbreaks=seq(-0.2, 1.8, 0.2), showBottom = RossSeaBathy, legend.text = 'B', xlab="", ylim = c(600, 0))
text(5,550,expression("Ammonium (\u03BCM)"), adj=0)
plot(RBgrid2, which="Nitrite (uM)", ztype = "image", zcol = cmocean('matter'), zbreaks=seq(-0.02, 0.15, 0.01), showBottom = RossSeaBathy, legend.text = 'C', xlab="", ylim = c(600, 0))
text(5,550,"Nitrite (\u03BCM)", adj=0)
plot(RBgrid2, which="Silicate (uM)", ztype = "image", zcol = cmocean('matter'), zbreaks=seq(60, 90, 2), showBottom = RossSeaBathy, legend.text = 'D', xlab="", ylim = c(600, 0))
text(5,550,"Silicate (\u03BCM)", adj=0)
plot(RBgrid2, which="N+N (uM)", ztype = "image", zcol = cmocean('matter'), zbreaks=seq(16, 34, 2), showBottom = RossSeaBathy, legend.text = 'E', ylim = c(600, 0))
text(5,550,"Sigmat", adj=0)
plot(RBgrid2, which="Chl a (ug/L)", ztype = "image", zcol = cmocean('algae'), zbreaks=seq(0, 7, 0.5), showBottom = RossSeaBathy, legend.text = 'F', ylim = c(600, 0))
text(5,550,"Chlorophyll a (\u03BCg/L)", adj=0)


# Close of the TIFF image to force the write to file
dev.off()



