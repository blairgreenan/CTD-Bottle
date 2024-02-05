# plot_RossSea_Bottle_Data.R
# Blair Greenan
# Fisheries and Oceans Canada
# 1 May 2023
#
# Description: this script generates a faceted plot with profile of dissolved iron
# from the 6 trace metal CTD casts that we carried out over Ross Bank. The first 5 
# casts were on Survey 1 and the final cast was at the Ross Bank summit on the second
# survey
#
# Update 29 Jan 2024: Added the ROMS model dFe profiles at the trace metal cast locations.
# Plotted as red lines on the facet plot with Pete's dFe data.
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
TM_cast <- unique(data[['cast Fe']])
# Get rid of NAN from the cast # list
TM_cast <- TM_cast[-which(is.nan(TM_cast))]

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
for (i in seq_along(TM_cast)) {
  cat('TMCTD cast', TM_cast[i], '...')
  II <- data[['cast Fe']] == TM_cast[i]
  ## create the ctd_bottle object
  tmp <- as.ctd(salinity=data$`sal00: Salinity, Practical [PSU]`[II],
                temperature=data$`t090C: Temperature [ITS-90, deg C]`[II],
                pressure=data$`prDM: Pressure, Digiquartz [db]`[II],
                station = data$`cast Fe`[II][1],
                cruise = "NBP1201",
                ship = "R/V Nathaniel B. Palmer") # just use the first element of the station vector in the data - oddity of the MAT file having a station # for each point in the CTD profile
  # add the other fields that were collected by the CTD system on the R/V NBP
  # I am going to ignore the Time vector in the data since it is not particularly useful
  fields <- names(data)
  fields <- fields[-which(fields %in% c('cast Fe', 'sal00: Salinity, Practical [PSU]', 't090C: Temperature [ITS-90, deg C]', 'prDM: Pressure, Digiquartz [db]'))]
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

# Water depth in metres compiled from nbp1201_ctd_list - BG.xlsx
water_depth = c(674, 172, 302, 297, 594, 170)

# Create a data frame for the dFe data collected at the 5 stations on Ross Bank
dFe_data <- list()
for (j in c(21, 22, 23, 24, 25)) {
#  ctd_bottle[[j]]@data$`dFe (nM)`[which(is.nan(ctd_bottle[[j]]@data$`dFe (nM)`))] <- NA
#  ctd_bottle[[j]]@data$pressure[which(is.nan(ctd_bottle[[j]]@data$pressure))] <- NA
  dFe_omit <- as.numeric(na.omit(ctd_bottle[[j]]@data$`dFe (nM)`))
  Press_omit <- as.numeric(na.omit(ctd_bottle[[j]]@data$`depth Fe`))
  dFe_data <- rbind(dFe_data, data.frame(Pressure=Press_omit,dFe=dFe_omit,Facet_cast=TM_cast[j], WD=water_depth[j-20], Survey="Survey 1"))
}

# Second occupation of Ross Bank (Station 74 near top of the bank)
#ctd_bottle[[32]]@data$`dFe (nM)`[which(is.nan(ctd_bottle[[32]]@data$`dFe (nM)`))] <- NA
#ctd_bottle[[32]]@data$pressure[which(is.nan(ctd_bottle[[32]]@data$pressure))] <- NA
#dFe_data <- rbind(dFe_data, data.frame(Pressure=ctd_bottle[[32]]@data$pressure,dFe=ctd_bottle[[32]]@data$`dFe (nM)`,Facet_cast=TM_cast[32], Survey="Survey 2"))
dFe_omit2 <- as.numeric(na.omit(ctd_bottle[[32]]@data$`dFe (nM)`))
Press_omit2 <- as.numeric(na.omit(ctd_bottle[[32]]@data$`depth Fe`))
dFe_data <- rbind(dFe_data, data.frame(Pressure=Press_omit2,dFe=dFe_omit2,Facet_cast=TM_cast[32], WD=water_depth[6], Survey="Survey 2"))

# create a depth above bottom variable
dFe_data$Above_bottom <- dFe_data$WD-dFe_data$Pressure

# Create pdf file
#pdf("dFE_TMCTD_cast.pdf", width=6, height=6)

############# Add the ROMS model dFe profiles
scale_factor <- 1/2.5 # factor to scale the model results
dFe_model <- list()
model26_data <- read.csv(file = "./ROMS_dFe/Model.TM26.012112.csv", sep = ",", header = FALSE)
model26_depth <- -1*model26_data$V1[model26_data$V4<10]
model26_dFe <- scale_factor*model26_data$V4[model26_data$V4<10]
dFe_model <- rbind(dFe_model, data.frame(Pressure=model26_depth,dFe=model26_dFe,Facet_cast=26))
model27_data <- read.csv(file = "./ROMS_dFe/Model.TM27.012112.csv", sep = ",", header = FALSE)
model27_depth <- -1*model27_data$V1[model27_data$V4<10]
model27_dFe <- scale_factor*model27_data$V4[model27_data$V4<10]
dFe_model <- rbind(dFe_model, data.frame(Pressure=model27_depth,dFe=model27_dFe,Facet_cast=27))
model28_data <- read.csv(file = "./ROMS_dFe/Model.TM28.012112.csv", sep = ",", header = FALSE)
model28_depth <- -1*model28_data$V1[model28_data$V4<10]
model28_dFe <- scale_factor*model28_data$V4[model28_data$V4<10]
dFe_model <- rbind(dFe_model, data.frame(Pressure=model28_depth,dFe=model28_dFe,Facet_cast=28))
model29_data <- read.csv(file = "./ROMS_dFe/Model.TM29.012112.csv", sep = ",", header = FALSE)
model29_depth <- -1*model29_data$V1[model29_data$V4<10]
model29_dFe <- scale_factor*model29_data$V4[model29_data$V4<10]
dFe_model <- rbind(dFe_model, data.frame(Pressure=model29_depth,dFe=model29_dFe,Facet_cast=29))
model30_data <- read.csv(file = "./ROMS_dFe/Model.TM30.012112.csv", sep = ",", header = FALSE)
model30_depth <- -1*model30_data$V1[model30_data$V4<10]
model30_dFe <- scale_factor*model30_data$V4[model30_data$V4<10]
dFe_model <- rbind(dFe_model, data.frame(Pressure=model30_depth,dFe=model30_dFe,Facet_cast=30))
model27_data <- read.csv(file = "./ROMS_dFe/Model.TM27.012112.csv", sep = ",", header = FALSE)
model27_depth <- -1*model27_data$V1[model27_data$V4<10]
model27_dFe <- scale_factor*model27_data$V4[model27_data$V4<10]
dFe_model <- rbind(dFe_model, data.frame(Pressure=model27_depth,dFe=model27_dFe,Facet_cast=37))


dev.new()
# Facet plot of the data
plt <- ggplot(data=dFe_data) + 
  geom_point(mapping = aes(x=dFe, y=Pressure)) + 
  facet_wrap(~Facet_cast, scales = "fixed") + 
  scale_y_reverse() +
#  xlab(expression(paste("Dissolved Iron (",~mu,"M)"))) +  # mistake on my part....units should be nM
  xlab("Dissolved Iron (nM)") +
  ylab("Depth (m)") +
  geom_hline(data=dFe_data, aes(yintercept=WD))

plt <- plt + geom_line(data = dFe_model, aes(x=dFe, y=Pressure), color="red")
  
plt

ggsave("RossBank_dFe.png", width = 10, height = 8, units = c("cm"), dpi = 1200, bg = "white", scale = 1.5)
ggsave("RossBank_dFe.pdf", width = 10, height = 8, units = c("cm"), dpi = 1200, bg = "white", scale = 1.5)

dev.off()

# Depth above bottom plot
dev.new()
ggplot(data=dFe_data) + 
  geom_point(mapping = aes(x=dFe, y=Above_bottom)) + 
  facet_wrap(~Facet_cast, scales = "fixed") + 
  xlab("Dissolved Iron (nM)") +
  ylab("Depth Above Bottom (m)")


ggsave("RossBank_dFe_above_bottom.png", width = 10, height = 8, units = c("cm"), dpi = 1200, bg = "white", scale = 1.5)
ggsave("RossBank_dFe_above_bottom.pdf", width = 10, height = 8, units = c("cm"), dpi = 1200, bg = "white", scale = 1.5)

dev.off()

# Create TIFF file
#tiff("dFE_TMCTD_cast.tiff", width=6, height=6, units='in', res=1200, compression = 'lzw')


