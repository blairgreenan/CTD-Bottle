# RossSea_Bottle_Chl_MVP_FL_calibration.R
# Blair Greenan
# Fisheries and Oceans Canada
# 31 Jan 2024
#
# Description: this script loads the NBP1201 CTD bottle data set and the MVP data.
# The script will loop through all the CTD casts and search for a nearby (in time/space)
# MVP profile. For each cast that meets the "nearby" criteria, the script will extract 
# the bottle depths and chlorophyll data, and will also extract the MVP fluorometer values
# at the same depths as the bottles. This paired chl-FL data set will then be used
# to create a linear fit of the relationship.
#
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

######## Load the MVP Rdata file
load("C:/Science Projects/Ross Sea/Documents/Papers/Ross Bank/Figures/MVP/MVP.RData")

# Adjust the longitudes as they or 0 to 360
neg_lon <- which(MVP_tidy_tibble$lon<0)
MVP_tidy_tibble$lon[neg_lon] <- 360 + MVP_tidy_tibble$lon[neg_lon]

######## End MVP Rdata file

# counter for arrays of chl data
mm <- 1

# create empty vectors
MVP_chl <- c()
CTD_chl <- c()

# Loop through each of the CTD casts
for (jj in 1:length(ctd_bottle)) {
  ctd_time <- ctd_bottle[[jj]]@metadata$startTime
  # Get the index values for the range of times
  MVP_survey_time <- which(MVP_tidy_tibble$daten>=(ctd_time-hours(2)) & MVP_tidy_tibble$daten<=(ctd_time+hours(2)))
  if (length(MVP_survey_time)>0) {
    MVP_survey_depth_match <- which(MVP_tidy_tibble$pres[MVP_survey_time] %in% ctd_bottle[[jj]]@data$`Nominal Depth (m)`)
    for (kk in 1:length(MVP_survey_depth_match)) {
      MVP_chl[mm] <- MVP_tidy_tibble$chl[MVP_survey_time[MVP_survey_depth_match[kk]]]
      # get index in the bottle data based on the depth
      bottle_index <- which(ctd_bottle[[jj]]@data$`Nominal Depth (m)` == MVP_tidy_tibble$pres[MVP_survey_time[MVP_survey_depth_match[kk]]])
      CTD_chl[mm] <- ctd_bottle[[jj]]@data$`Chl a (ug/L)`[bottle_index]
      mm <- mm + 1
    }
  }
}

# remove any instances of NaN that occur in the CTD chl bottle data
CTD_chl_rm_nan <- CTD_chl[-which(is.nan(CTD_chl))]
MVP_chl_rm_nan <- MVP_chl[-which(is.nan(CTD_chl))]

# create a tibble
chl_tibble <- tibble(CTD_chl_rm_nan, MVP_chl_rm_nan)
# Rename the 'original_name' variable to 'new_name'
chl_tibble <- chl_tibble %>%
  rename(MVP_Fluorometer = MVP_chl_rm_nan)
chl_tibble <- chl_tibble %>%
  rename(CTD_Chla = CTD_chl_rm_nan)

# Create a linear model
#chl_model <- lm(chl_tibble$MVP_Fluorometer ~ chl_tibble$CTD_Chla)
chl_model <- lm(chl_tibble$CTD_Chla ~ chl_tibble$MVP_Fluorometer)

# Print a summary of the linear regression
summary(chl_model)


# plot the data with ggplot
plt <- ggplot(data = chl_tibble, aes(x=MVP_Fluorometer, y=CTD_Chla)) + 
  geom_point() + ylab("Chlorophyll (\u03BCg/L)") +
  xlab("MVP Fluorometer") + xlim(c(-1, 25)) + ylim(c(-1, 25))
# geom_point() + xlim(c(-1, 25)) + ylim(c(-1, 25)) + xlab("Chlorophyll (\mug/L") +
  
# Overlay the linear fit using geom_smooth
plt <- plt + geom_smooth(method = "lm", se = FALSE, color = "red")

plt

ggsave("MVP_FL_calibration.png")


