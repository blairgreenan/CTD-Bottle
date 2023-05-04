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
cast <- unique(data[['CTD Cast #']])
# Get rid of NAN from the cast # list
cast <- cast[-which(is.nan(cast))]

# Edit the longitude so that it matches the convention in the Ross Sea bathymetry data
for (i in 1:length(data$`Longitude (decimal Deg)`)){
  # if the data is NAN, then skip over this
  if(!is.nan(data$`Longitude (decimal Deg)`[i])){
    if(data$`Longitude (decimal Deg)`[i]>0)
    {data$`Longitude (decimal Deg)`[i]<-data$`Longitude (decimal Deg)`[i]-360}
  }
}

# create an empty list for the ctd_bottle object
ctd_bottle <- list()

# loop to populate the ctd object with data and metadata for each CTD cast
for (i in seq_along(cast)) {
  cat('Cast', cast[i], '...')
  II <- data[['CTD Cast #']] == cast[i]
  ## create the ctd_bottle object
  tmp <- as.ctd(salinity=data$`sal00: Salinity, Practical [PSU]`[II],
                temperature=data$`t090C: Temperature [ITS-90, deg C]`[II],
                pressure=data$`prDM: Pressure, Digiquartz [db]`[II],
                station = data$`CTD Cast #`[II][1],
                cruise = "NBP1201",
                ship = "R/V Nathaniel B. Palmer") # just use the first element of the station vector in the data - oddity of the MAT file having a station # for each point in the CTD profile
  # add the other fields that were collected by the CTD system on the R/V NBP
  # I am going to ignore the Time vector in the data since it is not particularly useful
  fields <- names(data)
  fields <- fields[-which(fields %in% c('Station #', 'sal00: Salinity, Practical [PSU]', 't090C: Temperature [ITS-90, deg C]', 'prDM: Pressure, Digiquartz [db]'))]
  for (f in fields) {
    tmp <- oceSetData(tmp, f, data[[f]][II])
  }
  # add metadata latitude, longitude and start time of the CTD cast to the ctd object
  tmp <- oceSetMetadata(tmp, 'longitude', data$`Longitude (decimal Deg)`[i])
  tmp <- oceSetMetadata(tmp, 'latitude', data$`Latitude (decimal Deg)`[i])
  tmp <- oceSetMetadata(tmp, 'startTime', ymd(data$`Date  yyyymmdd`[i]) + hms(paste(floor(data$`UTC Time hhmm`[i]/100), data$`UTC Time hhmm`[i]%%100, "00", sep=":")))
  ctd_bottle[[i]] <- tmp
  cat('\n')
  
}

# Creat a data frame for the dFe data collected at the 5 stations on Ross Bank
dFe_data <- list()
for (j in c(36, 38, 40, 44, 48)) {
  ctd_bottle[[j]]@data$`dFe (nM)`[which(is.nan(ctd_bottle[[j]]@data$`dFe (nM)`))] <- NA
  ctd_bottle[[j]]@data$pressure[which(is.nan(ctd_bottle[[j]]@data$pressure))] <- NA
  dFe_data <- rbind(dFe_data, data.frame(Pressure=ctd_bottle[[j]]@data$pressure,dFe=ctd_bottle[[j]]@data$`dFe (nM)`,Cast=cast[j], Survey="Survey 1"))
}

# Second occupation of Ross Bank (Station 74 near top of the bank)
ctd_bottle[[74]]@data$`dFe (nM)`[which(is.nan(ctd_bottle[[74]]@data$`dFe (nM)`))] <- NA
ctd_bottle[[74]]@data$pressure[which(is.nan(ctd_bottle[[74]]@data$pressure))] <- NA
dFe_data <- rbind(dFe_data, data.frame(Pressure=ctd_bottle[[74]]@data$pressure,dFe=ctd_bottle[[74]]@data$`dFe (nM)`,Cast=cast[74], Survey="Survey 2"))

# Create pdf file
pdf("dFE.pdf", width=6, height=6)

# Facet plot of the data
ggplot(data=dFe_data) + 
  geom_point(mapping = aes(x=dFe, y=Pressure)) + 
  facet_wrap(~Cast, scales = "fixed") + 
  scale_y_reverse() +
  xlab(expression(paste("Dissolved Iron (",~mu,"M)"))) +
  ylab("Pressure (db)")

dev.off()

# Create TIFF file
tiff("dFE.tiff", width=6, height=6, units='in', res=1200, compression = 'lzw')

# Facet plot of the data
ggplot(data=dFe_data) + 
  geom_point(mapping = aes(x=dFe, y=Pressure)) + 
  facet_wrap(~Cast, scales = "fixed") + 
  scale_y_reverse() +
  xlab(expression(paste("Dissolved Iron (",~mu,"M)"))) +
  ylab("Pressure (db)")

dev.off()
