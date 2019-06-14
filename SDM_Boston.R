# Author: Patrick Sogno
# E-mail: patrick.sogno@stud-mail.uni-wuerzburg.de
# Date: 05 June 2019
# RStudio version: 1.1.456
# R version: 3.5.1

#--------------------------------------------------#
# Background: This code is part of a university
# project on spatial modeling and prediction.
# The objective to spatially model user types for
# rental bikes in Boston, Massachusetts, USA has
# been developed by Ronja Lappe and me.
# This code and Ronja's are part of the same work.
# For any questions regarding her code, please
# contact her directly:
# ronja.lappe@stud-mail.uni-wuerzburg.de
#--------------------------------------------------#

# => Hypothesis: Based on the bike data of Boston for [DATE], [TIMESPAN] we can:
# 1. Find spatial distribution at that point in time for different user types
# 2. Build an occurence model that is transferable to other times.
# 3. We will find that bike usage is not correlated to seasons.
# 4. We will find that bike usage is not correlated to specific demographic 
# 5. user types
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#
# In this code:
#
# Spatial distribution modelling.
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


######------------------------ 0. Set environment ------------------------######

## Paths:

ipath ="D:/Patrick/Documents/Dokumente/SoSe2019/MET/SpatialModeling/BikeData_bo"  
opath = getwd()

bikeDir <- paste0(opath, "/", "BikeData_bo")

## Packages:

library('dplyr')
library('rgdal')
library('raster')
library('sdm')
library('data.table')

## Data names:

distanceStack <- "distanceStack.tif"
distanceStation <- "distanceStation.tif"
NDVIraster <- "NDVI.tif"
tripdata <- "201705-hubway-tripdata.csv"
stations <- "previous_Hubway_Stations_as_of_July_2017.csv"



######------------------------- End of section 0 -------------------------######



######------------------------ 1. Read occurences ------------------------######


## Subset data to only one day:

# read file
trips <- read.csv(paste0(ipath, "/", tripdata))
head(trips)
class(trips$starttime)

# subset df: Only trips that were on May 1, 2017:

trips$starttime <- as.character(trips$starttime) # factor to character
temp <- strsplit(trips$starttime, " ") # split string
mat <- matrix(unlist(temp), ncol=2, byrow=TRUE) # build matrix from strings
tempdf <- as.data.frame(mat) # build df from matrix
colnames(tempdf) <- c("Date", "Time")
head(tempdf)

trips$date <- tempdf$Date # make new column with only start date
trips$time <- tempdf$Time # make new column with only start time

mytrips <- trips[which(trips$date == "2017-05-01"),] # select trips of May, 1

mytrips$time <- as.character(mytrips$time) # same as above but for time...
temp <- strsplit(mytrips$time, ":")
mat <- matrix(unlist(temp), ncol=3, byrow=TRUE)
tempdf <- as.data.frame(mat)
colnames(tempdf) <- c("h", "m", "s")
mytrips$time <- as.numeric(tempdf$h)

mytrips <- mytrips[which(mytrips$time == 9),] # select trips that occured at 9h
mytrips

mysubs <- mytrips[which(mytrips$usertype == "Subscriber"),] # split by usertype
mycusS <- mycusE <- mytrips[which(mytrips$usertype == "Customer"),] # ...
mycusS



## Spatial prediction


# read predictors
preds <- brick(paste0(bikeDir, "/", distanceStack))
plot(preds[[2]])
mycrs <- crs(preds)
mycrs
# Make 2 spdfs (start points, end points):
coordinates(mycusS) <- c("start.station.longitude", "start.station.latitude")
proj4string(mycusS) <- CRS(as.character(mycrs))
class(mycusS)
plot(mycusS, col = "red")
#
coordinates(mycusE) <- c("end.station.longitude", "end.station.latitude")
proj4string(mycusE) <- CRS(as.character(mycrs))
class(mycusE)
plot(mycusE, add = T, col = "blue")
#

# read occurences
occS <- mycusS
occE <- mycusE


# transform CRS of spdfs
plot(preds[[1]])
occS <- spTransform(occS, preds@crs)
plot(occS, col = "red", add = T)
occE <- spTransform(occE, preds@crs)
plot(occE, col = "blue", add = T)

# read stations
mystations <- read.csv(paste0(bikeDir, "/", stations))
coordinates(mystations) <- c("Longitude", "Latitude")
proj4string(mystations) <- CRS(as.character(mycrs))
plot(mystations, add = T)

# create absence data
'%!in%' <- function(x,y)!('%in%'(x,y)) # build "not in" function
occS$start.station.name
stationNames <- mystations$Station
occSS <- mystations[which(mystations$Station %in% occS$start.station.name), ]
plot(occSS)
absSS <- mystations[which(mystations$Station %!in% occS$start.station.name), ]
plot(absSS)
plot(occS, col = "red", add = T)
#
startSpecies <- mystations
startSpecies$occ <- NA
for (i in 1:length(startSpecies)) {
  if (startSpecies$Station[i] %!in% occS$start.station.name) {
    startSpecies$occ[i] = 0
  }else {
    startSpecies$occ[i] = 1
  }
}
#
plot(startSpecies[startSpecies$occ == 1,], col = "green")
plot(startSpecies[startSpecies$occ == 0,], col = "red", add = T)

#

## Train SDM:
names(occS)
d <- sdmData(formula = occ~layer, train = startSpecies, predictors = preds) # produce data for spatial prediction
d

m1 <- sdm(students~., data = d, methods = c('glm', 'svm')) # fit species distribution model

p1 <- predict(m1, newdata = preds, filename = paste0(opath, "/", "sdm_preds_1.grd"), overwrite = T) # predict

plot(p1)

