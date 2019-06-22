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
loadandinstall('plotrix')
loadandinstall('scales')

## parameters 
city <- 'bo' # Bosten(bo), Chicago(ch), Washington, D.C.(dc), Los Angeles(la), London(lo), 
# Minnesota(mn), New York City(ny), Philadelphia(ph), San Fransico(sf)
city_osm <- 'boston' 
dates <- 201701:201712 # only May 2017, timespan also possible, e.g. 201705:201708

mycrs <- "+proj=longlat +datum=WGS84"

######------------------------- 1. Find, download & load data -----------------------######

## bike data (if more than one month... use list and write as loop)
bikedata_dl_zip <- as.list(dl_bikedata(city=city,data_dir = dirs$main_dir, 
                                       dates = dates))                                       # download
bikedata_dl <- lapply(bikedata_dl_zip, function(x) unzip(x,exdir = dirs$main_dir))                                 # unzip
bikedata <- lapply(bikedata_dl, function(x) read.csv(file.path(dirs$main_dir,basename(x))))  # read into R

### extract stations 
store_bikedata(bikedb = 'bikedb', city = city, data_dir = dirs$main_dir, dates = dates) # create database of downloaded bikedata
stations <- bike_stations(bikedb = 'bikedb', city = city)       # extract stations for selected city from database
stations <- as.data.frame(stations)

## OSM land use (could also be simplified I guess)
# railway stations
preds <- new.env()                                                     # create new environment for predictors 
preds$c1 <- opq(bbox = city_osm) %>%                                   # query from OSM
  add_osm_feature(key = "public_transport") %>%
  add_osm_feature(key = "railway") %>%
  osmdata_xml(file = "OSMcommuter.osm", quiet = F)                     # save on harddisk
preds$c1 <- sf::st_read('OSMcommuter.osm', layer = 'points', quiet = F)# read into R
preds$c1 <- as(preds$c1, 'Spatial')                                    # convert to spatial
plot(preds$c1,col="#009900")
points(preds$c1,col="#009900",pch=20, add=T)

# touristic places
preds$t <- opq(bbox = city_osm) %>%                     
  add_osm_feature(key = "tourism") %>%
  osmdata_xml(file = "OSMtourist.osm", quiet = F)
preds$t <- sf::st_read('OSMtourist.osm', layer = 'points', quiet = F)   
preds$t <- as(preds$t, 'Spatial')
points(preds$t, col="red",pch=20,add=T)

# Offices
preds$o <- opq(bbox = city_osm) %>%                     
  add_osm_feature(key = "office") %>%
  osmdata_xml(file = "OSMoffice.osm", quiet = F)
preds$o <- sf::st_read('OSMoffice.osm', layer = 'points', quiet = F)
preds$o <- as(preds$o, 'Spatial')
plot(preds$o, col="red",add=T)
points(preds$o,col="blue",add=T,pch=20)

# university/college buildings 
preds$u <- opq(bbox = city_osm) %>%                     
  add_osm_feature(key = "amenity", value = c("university","college","research_institute")) %>%
  osmdata_xml(file = "OSMuniversity.osm", quiet = F)
preds$u <- sf::st_read('OSMuniversity.osm', layer='points',quiet = F)
preds$u <- as(preds$u, 'Spatial')
points(preds$u,col="purple",add=T,pch=20)

available_features()

# possible further predictors? 

######------------------------------- 2. Data preparation ---------------------------######

#--------------------------------- 2.1 Rasterize predictors ------------------------------#

# (could probably also be done more efficiently using lists and for-loops ;))

r <- new.env()                                # create new environment for distance rasters
stations_sp <- stations
coordinates(stations_sp) <- c("longitude","latitude")
proj4string(stations_sp) <- CRS(mycrs)
e <- extent(stations_sp)                      # use spatial extent of projected bike stations

# railway stations 
r$c1 <- raster(e, ncols = 100, nrows = 100) # For commuter land use 1
proj4string(r$c1) <- CRS("+proj=longlat +datum=WGS84")
r$c1 <- distanceFromPoints(object = r$c1, xy = preds$c1)

# touristic places 
r$t <- raster(e, ncols = 100, nrows = 100) # For tourist land use
proj4string(r$t) <- CRS("+proj=longlat +datum=WGS84")
r$t <- distanceFromPoints(object = r$t, xy = preds$t)
plot(r$t)

# offices
r$o <- raster(e, ncols = 100, nrows = 100) # For leisure land use
proj4string(r$o) <- CRS("+proj=longlat +datum=WGS84")
r$o <- distanceFromPoints(object = r$o, xy = preds$o)
plot(r$o)

# universities
r$u <- raster(e, ncols = 100, nrows = 100) # For leisure land use
proj4string(r$u) <- CRS("+proj=longlat +datum=WGS84")
r$u <- distanceFromPoints(object = r$u, xy = preds$u)
plot(r$u)

# create raster stack 
r_l <- as.list(r)
r_stack <- stack(r_l)
writeRaster(r_stack,file.path(dirs$main_dir,"distanceStack.tif"),overwrite=T)

#-------------------------------- 2.2 Manipulate bike dataset ----------------------------#

# create one occurence dataframe from list of bikedata frames 
for (i in 1:length(bikedata)) {
  if (i == 1) {
    occ <- bikedata[[i]]
  }else{
    occ <- rbind(occ, bikedata[[i]])
  }
}

# split date-time column
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

# trip start at 8am, 4 pm and 8pm during the week... for trips less than 1 hour...
my_occ_s <- list()
my_occ_s[[1]] <- occ_s[which(occ_s$start_hour == 8 & occ_s$start_day <=5 &
                               occ_s$tripduration < 3600 & occ_s$start_month == 5),]
my_occ_s[[2]] <- occ_s[which(occ_s$start_hour == 16 & occ_s$start_day <=5 &
                               occ_s$tripduration < 3600& occ_s$start_month == 5),]
my_occ_s[[3]] <- occ_s[which(occ_s$start_hour == 20 & occ_s$start_day <=5 &
                               occ_s$tripduration < 3600& occ_s$start_month == 5),]

# trip ends 
my_occ_e <- list()
names(occ_e)
my_occ_e[[1]] <- occ_e[which(occ_e$stop_hour == 8 & occ_e$start_day <=5 &
                               occ_e$tripduration < 3600& occ_s$start_month == 5),]
my_occ_e[[2]] <- occ_e[which(occ_e$stop_hour == 16 & occ_e$start_day <=5 &
                               occ_e$tripduration < 3600& occ_s$start_month == 5),]
my_occ_e[[3]] <- occ_e[which(occ_e$stop_hour == 20 & occ_e$start_day <=5 &
                               occ_e$tripduration < 3600& occ_s$start_month == 5),]

#----------------------- 2.4 Create presence and absence at stations ---------------------#
# trip starts
pres_s <- list()
for(i in 1:length(my_occ_s)){
  frq <- table(my_occ_s[[i]]$start.station.name)
  frq <- as.data.frame(frq)
  names(frq) <- c("name","freq")
  st_frq <- base::merge(stations,frq,all=T)
  st_frq <- na.omit(st_frq)
  coordinates(st_frq) <- c("longitude","latitude")
  proj4string(st_frq) <- CRS(mycrs)
  pres_s[[i]] <- st_frq
}

# trip ends 
pres_e <- list()
for(i in 1:length(my_occ_e)){
  frq <- table(my_occ_e[[i]]$end.station.name)
  frq <- as.data.frame(frq)
  names(frq) <- c("name","freq")
  st_frq <- base::merge(stations,frq,all=T)
  st_frq <- na.omit(st_frq)
  coordinates(st_frq) <- c("longitude","latitude")
  proj4string(st_frq) <- CRS(mycrs)
  pres_e[[i]] <- st_frq
}


#--------------------------------- 2.4 Check for collinearity ----------------------------#
# visual inspection of collinearity
cm <- cor(getValues(r_stack), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))

#----------------------------- 3. Species distribution model -----------------------------#

# trip starts
p_s <- list()
for(i in 1:length(pres_s)){
  d <- sdmData(formula = freq ~ c1+o+t+u, train = pres_s[[i]], predictors = r_stack)
  m <- sdm(freq~., data = d, methods = c('rf','svm')) # fit species distribution model
  p_s[[i]] <- predict(m, newdata = r_stack, overwrite = T) # predict
}

# trip ends 
p_e <- list()
for(i in 1:length(pres_e)){
  d <- sdmData(formula = freq ~ c1+o+t+u, train = pres_e[[i]], predictors = r_stack)
  m <- sdm(freq~., data = d, methods = c('rf','svm')) # fit species distribution model
  p_e[[i]] <- predict(m, newdata = r_stack, overwrite = T) # predict
}

#--------------------------------------- 4. Visualization --------------------------------------#

# base map 
e_l <- as.list(e) # save station extent as list to use it as input for the openmap 
map <- openmap(upperLeft = e_l[c(4,1)], lowerRight = e_l[c(3,2)], type = "opencyclemap", mergeTiles = T)
map <- openproj(map,projection = mycrs)
plot(map,add=F)

# frequency per station 
plot(map,add=F)
points(pres_s[[1]],cex=pres_s[[1]]$freq/10,col=alpha("blue",0.4),pch=20) # start wdays 8 am
plot(map,add=F)
points(pres_e[[1]],cex=pres_e[[1]]$freq/10,col=alpha("blue",0.4),pch=20) # end wdays 8 am
points(pres_s[[2]],cex=pres_s[[2]]$freq/10,col=alpha("blue",0.4),pch=20) # start wdays 4 pm
points(pres_s[[3]],cex=pres_s[[3]]$freq/10,col=alpha("blue",0.4),pch=20) # start wdays 8 pm


# predicted frequency 
plot(map,add=F)
plot(p_s[[1]]$id_2.sp_1.m_svm, alpha=0.5,add=T)
plot(map,add=F)
plot(p_e[[1]]$id_2.sp_1.m_svm, alpha=0.5,add=T)

plot(p_s[[2]]$id_2.sp_1.m_svm)
plot(p_s[[3]]$id_2.sp_1.m_svm)
plot(p_s[[4]]$id_2.sp_1.m_svm,alpha=0.5,add=T)

p_s_stack <- stack(p_s[[1]],p_s[[2]],p_s[[3]])
plotRGB(p_s_stack,1,3,5,stretch="lin",alpha=200,add=T)