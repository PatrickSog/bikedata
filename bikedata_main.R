####################### Spatial prediction of rented bikes in Bosten ######################

# Author: Patrick Sogno and Ronja Lappe
# Email: patrick.sogno@stud-mail.uni-wuerzburg.de/ ronja.lappe@posteo.de

#----------------------------------------Background---------------------------------------#

# Background: This code is part of a university project on spatial modeling and prediction.

#-------------------------------------Problem formulation---------------------------------#

# => Problem formulation: Where and when do we find a high density of rented & 
# returned bikes? Which land use types seem to detemermine this distribution? 
# can we predict the occurance of rented and returned bikes for other times or cities? 

# => Hypothesis: 
# 1. Bike rentals and returns are generally highest at public transport stations (pts)
# 2. Bike rentals at touristic places are highest in the afternoons and on weekends
# 3. Bike returns in the mornings are highest in business districts
# 4. Bike rentals in the mornings are highest at pts
# 5. Bike returns in the evenings are highest at pts
# 6. Bike rentals in the evenings are highest at business districts… ?
# 7. ...? 

#----------------------------------------- Data sets -------------------------------------#

# rental stations 
# bike rentals & returns (station, time, trip duration, gender, age, and more)
# OSM data (public transport stations (railway, bus), tourisim, leisure, buisness?)
# evtl. additional climate data...

#------------------------------------ Content of this code -------------------------------#

# 1. data download
# 2. data preparation 
# 3. species distribution modelling 
# 4. visualization 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#+++++++++++#

######------------------------------- 0. Set environment ----------------------------######

## paths
dirs <- new.env()
dirs$main_dir <- "/Users/Ronjamac/Documents/Studium_Geographie/Master_EAGLE/Spatial_modelling_prediction/Spatial_modeling_prediction/bikedata"

## packages 
loadandinstall <- function(mypkg) {if (!is.element(mypkg, installed.packages()[,1])){install.packages(mypkg)}; library(mypkg, character.only=TRUE)  }
loadandinstall('bikedata')
loadandinstall('dplyr')
loadandinstall('rgdal')
loadandinstall('raster')
loadandinstall('sdm')
loadandinstall('data.table')
loadandinstall('dodgr')
loadandinstall('tidyverse')
loadandinstall('osmdata')
loadandinstall('osmplotr')
loadandinstall('OpenStreetMap')
loadandinstall('ellipse')
loadandinstall('lubridate')

## parameters 
city <- 'bo' # Bosten(bo), Chicago(ch), Washington, D.C.(dc), Los Angeles(la), London(lo), 
             # Minnesota(mn), New York City(ny), Philadelphia(ph), San Fransico(sf)
city_osm <- 'boston' 
dates <- 201705 # only May 2017, timespan also possible, e.g. 201705:201708

mycrs <- "+proj=longlat +datum=WGS84"

######------------------------- 1. Find, download & load data -----------------------######

## bike data (if more than one month... use list and write as loop)
bikedata_dl_zip <- dl_bikedata(city=city,data_dir = dirs$main_dir, 
                               dates = dates)                         # download
bikedata_dl <- unzip(bikedata_dl)                                     # unzip
bikedata <- read.csv(file.path(dirs$main_dir,basename(bikedata_dl)))  # read into R

### extract stations 
store_bikedata(bikedb = 'bikedb', city = city, data_dir = dirs$main_dir, dates = dates) # create database of downloaded bikedata
bike_stations <- bike_stations(bikedb = 'bikedb', city = city)       # extract stations for selected city from database
file_bike_stations <- paste(dates,city,"bike_stations.csv",sep = "_")# create filename
write.csv(stations_dl,file = file_bike_stations)                     # safe on harddisk
stations <- read.csv(file_bike_stations)                             # read into R
coordinates(stations) <- c("longitude", "latitude")                  # convert to spatial
proj4string(stations) <- CRS("+proj=longlat +datum=WGS84")           # add coordinate system  
plot(stations,col="blue",add=T)

## OSM land use (could also be simplified I guess)
  # railway stations
  preds <- new.env()                                                     # create new environment for predictors 
  preds$c1 <- opq(bbox = city_osm) %>%                                   # query from OSM
      add_osm_feature(key = "public_transport") %>%
      add_osm_feature(key = "railway") %>%
      osmdata_xml(file = "OSMcommuter.osm", quiet = F)                   # save on harddisk
  preds$c1 <- sf::st_read('OSMcommuter.osm', layer = 'points', quiet = F)# read into R
  preds$c1 <- as(preds$c1, 'Spatial')                                    # convert to spatial
  plot(preds$c1,col="#009900")
  points(preds$c1,col="#009900",pch=20, add=T)

  # bus stops
  preds$c2<- opq(bbox = city_osm) %>%                     
     add_osm_feature(key = "public_transport") %>%
     osmdata_xml(file = "OSMcommuter2.osm", quiet = F)
  preds$c2 <- sf::st_read('OSMcommuter2.osm', layer = 'points', quiet = F)   
  preds$c2 <- as(preds$c2, 'Spatial')
  plot(preds$c2, col = "#FFCC00", add = T)
  points(preds$c2, col = "#FFCC00", pch=20, add = T)
  
  # touristic places
  preds$t <- opq(bbox = city_osm) %>%                     
     add_osm_feature(key = "tourism") %>%
     osmdata_xml(file = "OSMtourist.osm", quiet = F)
  preds$t <- sf::st_read('OSMtourist.osm', layer = 'points', quiet = F)   
  preds$t <- as(preds$t, 'Spatial')
  plot(preds$t, add=T)

  # Leisure-related places
  preds$l <- opq(bbox = city_osm) %>%                     
     add_osm_feature(key = "leisure") %>%
     osmdata_xml(file = "OSMleisure.osm", quiet = F)
  preds$l <- sf::st_read('OSMleisure.osm', layer = 'points', quiet = F)
  preds$l <- as(preds$l, 'Spatial')
  plot(preds$l, col="green",add=T)

  # Offices
  preds$o <- opq(bbox = city_osm) %>%                     
     add_osm_feature(key = "office") %>%
      osmdata_xml(file = "OSMoffice.osm", quiet = F)
  preds$o <- sf::st_read('OSMoffice.osm', layer = 'points', quiet = F)
  preds$o <- as(preds$o, 'Spatial')
  plot(preds$o, col="red",add=T)
  points(preds$o,col="blue",add=T,pch=20)


######------------------------------- 2. Data preparation ---------------------------######

#--------------------------------- 2.1 Rasterize predictors ------------------------------#

# (could probably also be done more efficiently using lists and for-loops ;))

r <- new.env()                                # create new environment for distance rasters
e <- extent(stations)                         # use spatial extent of projected bike stations
  # railway stations 
  r$c1 <- raster(e, ncols = 100, nrows = 100) # For commuter land use 1
  proj4string(r$c1) <- CRS("+proj=longlat +datum=WGS84")
  r$c1 <- distanceFromPoints(object = r$c1, xy = preds$c1)

  # bus stops 
  r$c2 <- raster(e, ncols = 100, nrows = 100) # For commuter land use 2
  proj4string(r$c2) <- CRS("+proj=longlat +datum=WGS84")
  r$c2 <- distanceFromPoints(object = r$c2, xy = preds$c2)
  plot(r$c2)
  
  # touristic places 
  r$t <- raster(e, ncols = 100, nrows = 100) # For tourist land use
  proj4string(r$t) <- CRS("+proj=longlat +datum=WGS84")
  r$t <- distanceFromPoints(object = r$t, xy = preds$t)
  plot(r$t)
  
  # leisure places
  r$l <- raster(e, ncols = 100, nrows = 100) # For leisure land use
  proj4string(r$l) <- CRS("+proj=longlat +datum=WGS84")
  r$l <- distanceFromPoints(object = r$l, xy = preds$l)
  plot(r$l)
  
  # offices
  r$o <- raster(e, ncols = 100, nrows = 100) # For leisure land use
  proj4string(r$o) <- CRS("+proj=longlat +datum=WGS84")
  r$o <- distanceFromPoints(object = r$o, xy = preds$o)
  plot(r$o)
  
  # create raster stack 
  r_stack <- stack(r$c1,r$c2,r$l,r$s,r$t,r$o)
  names(r_stack) <- c("c1","c2","l","s","t","o")
  writeRaster(r_stack,"distanceStack.tif",overwrite=T)
  
## Rasterize bike stations 
r$s <- raster(e, ncols = 100, nrows = 100) # For distance to bike stations
proj4string(r$s) <- CRS("+proj=longlat +datum=WGS84")
r$s <- distanceFromPoints(object = r$c2, xy = stations)
plot(r$s)
writeRaster(r$s,"distanceStation.tif")

#-------------------------------- 2.2 Manipulate bike dataset ----------------------------#

# split date-time column
occ <- bikedata 
occ$starttime <- ymd_hms(occ$starttime)
occ$stoptime <- ymd_hms(occ$stoptime)
occ$start_year <- year(occ$starttime)
occ$start_month <- month(occ$starttime)
occ$start_day <- day(occ$starttime)
occ$start_wday <- wday(occ$starttime)
occ$start_hour <- hour(occ$starttime)
occ$stop_year <- year(occ$stoptime)
occ$stop_month <- month(occ$stoptime)
occ$stop_day <- day(occ$stoptime)
occ$stop_wday <- wday(occ$stoptime)
occ$stop_hour <- hour(occ$stoptime)

# create 2 spatial dataframe with trip_start and trip_end points 
occ_s <- occ
coordinates(occ_s) <- c("start.station.longitude", "start.station.latitude")
proj4string(occ_s) <- CRS(mycrs)

occ_e <- occ
coordinates(occ_e) <- c("end.station.longitude", "end.station.latitude")
proj4string(occ_e) <- CRS(mycrs)

#------------------------------ 2.3 Subset bike dataset by time --------------------------#

# FOR EXAMPLE: 1st of May at 8 am and 4 pm
occ_s_8am <- occ_s[which(occ_s$start_hour==8 & occ_s$start_day==1),]
occ_e_8am <- occ_e[which(occ_e$start_hour==8 & occ_e$start_day==1),]
occ_4pm <- occ[which(occ$start_hour==16 & occ$start_day==1),]

# OR: rented bikes on mondays at 8 am 
occ_8am_mon <- occ[which(occ$start_hour==8 & occ$start_wday==1),]

# saturday afternoons
occ_s_4pm_sat <- occ_s[which(occ_s$start_hour==16 & occ_s$start_wday==6),]

# and so on... according to what's interesting for our analysis 


#----------------------- 2.4 Create presence and absence at stations ---------------------#

# we have to work on this... right now, we only know, if one or several bikes are present
# at a certain station, but we don't now anything about the density (which would be 
# essential for our analysis, I guess). 
# it would further be nice to create an occurance data.frame which can directly be 
# subsetted by time (wday,hour, etc.)
# any ideas?? 

'%!in%' <- function(x,y)!('%in%'(x,y)) # build "not in" function

# all start trip data 
pres_all <- stations
for (i in 1:length(pres)) {
  if (pres_all$name[i] %!in% occ_s$start.station.name) {
    pres_all$occ_s[i] = 0
  }else {
    pres_all$occ_s[i] = 1
  }
}
head(pres_all,n=50) # lots of absence points even though all trip starts in may are considered... how is that possible? 

# start trip at 8 on 1st May
pres <- stations
pres$occ_s_8am <- NA
for (i in 1:length(pres)) {
  if (pres$name[i] %!in% occ_s_8am$start.station.name) {
    pres$occ_s_8am[i] = 0
  }else {
    pres$occ_s_8am[i] = 1
  }
}
# end trip at 8 on 1st of May
pres$occ_e_8am <- NA
for (i in 1:length(pres)) {
  if (pres$name[i] %!in% occ_e_8am$end.station.name) {
    pres$occ_e_8am[i] = 0
  }else {
    pres$occ_e_8am[i] = 1
  }
}

# trip start saturday afternoon
pres$occ_s_4pm_sat <- NA
for (i in 1:length(pres)) {
  if (pres$name[i] %!in% occ_s_4pm_sat$start.station.name) {
    pres$occ_s_4pm_sat[i] = 0
  }else {
    pres$occ_s_4pm_sat[i] = 1
  }
}

pres

plot(pres[pres_all$occ_s == 1,], col = "black")
plot(pres[pres_all$occ_s== 0,], col = "red", add = T)

# just try some plotting (not relevant for the analysis)
e_l <- as.list(e) # save station extent as list to use it as input for the openmap 
map <- openmap(upperLeft = e_l[c(4,1)], lowerRight = e_l[c(3,2)], type = "opencyclemap", mergeTiles = T)
map <- openproj(map,projection = mycrs)
plot(map)
points(pres[pres$occ_s_8am ==1,], pch=20, col="black")
points(pres[pres$occ_s_8am ==0,], pch=20, col="red")

#--------------------------------- 2.4 Check for collinearity ----------------------------#
# visual inspection of collinearity
cm <- cor(getValues(r_stack), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))

#----------------------------- 3. Species distribution model -----------------------------#

# the following models show a general occurance probability of trip_starts and trip_ends
# at the rentals stations based on the chosen times and predictors...

# subsetting of times should be improved and as said before, densitiy information might
# contain more relevant information... 

# all start trip data (only regarding train stations and offices)
d <- sdmData(formula = occ_s ~ c1+o, train = pres_all, predictors = r_stack)
m1 <- sdm(occ_s~., data = d, methods = c('glm', 'svm')) # fit species distribution model
p1 <- predict(m1, newdata = r_stack, overwrite = T) # predict
plot(p1$id_1.sp_1.m_glm)

# start trip 1st of may 8 am 
d_1 <- sdmData(formula = occ_8am ~ c1+o, train = pres, predictors = r_stack)
m1_1 <- sdm(occ_s_8am~., data = d_1, methods = c('glm', 'svm')) # fit species distribution model
p1_1 <- predict(m1_1, newdata = r_stack, overwrite = T) # predict
plot(p1_1$id_1.sp_1.m_glm)

# end trip 1st of may 8 am 
d_2 <- sdmData(formula = occ_e_8am ~ c1+o, train = pres, predictors = r_stack)
m1_2 <- sdm(occ_e_8am~., data = d_2, methods = c('glm', 'svm')) # fit species distribution model
p1_2 <- predict(m1_2, newdata = r_stack, overwrite = T) # predict
plot(p1_2$id_1.sp_1.m_glm)

# trip start saturday afternoons
d_3 <- sdmData(formula = occ_s_4pm_sat ~ c1+o+t, train = pres, predictors = r_stack)
m1_3 <- sdm(occ_s_4pm_sat~., data = d_3, methods = c('glm', 'svm')) # fit species distribution model
p1_3 <- predict(m1_3, newdata = r_stack, overwrite = T) # predict
plot(p1_3$id_1.sp_1.m_glm)
