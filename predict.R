#------------------------------- 5. Predict for another city ---------------------------------#
## NEW PREDICTORS
city_osm_1 <- 'washington'
e_1 <- getbb('washington')

## OSM land use (could also be simplified I guess)
# railway stations
preds_1 <- new.env()                                                     # create new environment for predictors 
preds_1$c1 <- opq(bbox = city_osm_1) %>%                                   # query from OSM
  add_osm_feature(key = "public_transport") %>%
  add_osm_feature(key = "railway") %>%
  osmdata_xml(file = "OSMcommuter_1.osm", quiet = F)                     # save on harddisk
preds_1$c1 <- sf::st_read('OSMcommuter_1.osm', layer = 'points', quiet = F)# read into R
preds_1$c1 <- as(preds_1$c1, 'Spatial')                                    # convert to spatial
plot(preds_1$c1,col="#009900")
points(preds_1$c1,col="#009900",pch=20, add=T)

# touristic places
preds_1$t <- opq(bbox = city_osm_1) %>%                     
  add_osm_feature(key = "tourism") %>%
  osmdata_xml(file = "OSMtourist_1.osm", quiet = F)
preds_1$t <- sf::st_read('OSMtourist_1.osm', layer = 'points', quiet = F)   
preds_1$t <- as(preds_1$t, 'Spatial')
points(preds_1$t, col="red",pch=20,add=T)

# Offices
preds_1$o <- opq(bbox = city_osm_1) %>%                     
  add_osm_feature(key = "office") %>%
  osmdata_xml(file = "OSMoffice_1.osm", quiet = F)
preds_1$o <- sf::st_read('OSMoffice_1.osm', layer = 'points', quiet = F)
preds_1$o <- as(preds_1$o, 'Spatial')
plot(preds_1$o, col="red",add=T)
points(preds_1$o,col="blue",add=T,pch=20)

# university/college buildings 
preds_1$u <- opq(bbox = city_osm_1) %>%                     
  add_osm_feature(key = "amenity", value = c("university","college","research_institute")) %>%
  osmdata_xml(file = "OSMuniversity_1.osm", quiet = F)
preds_1$u <- sf::st_read('OSMuniversity_1.osm', layer='points',quiet = F)
preds_1$u <- as(preds_1$u, 'Spatial')
points(preds_1$u,col="purple",add=T,pch=20)

#--------------------------------- 5.1 Rasterize predictors ------------------------------#

# (could probably also be done more efficiently using lists and for-loops ;))

r_1 <- new.env()                                # create new environment for distance rasters
e_1 <- extent(e_1)                     # use spatial extent of projected bike stations
# railway stations 
r_1$c1 <- raster(e_1, ncols = 100, nrows = 100) 
proj4string(r_1$c1) <- CRS("+proj=longlat +datum=WGS84")
r_1$c1 <- distanceFromPoints(object = r_1$c1, xy = preds_1$c1)
plot(r_1$c1)

# touristic places 
r_1$t <- raster(e_1, ncols = 100, nrows = 100) 
proj4string(r_1$t) <- CRS("+proj=longlat +datum=WGS84")
r_1$t <- distanceFromPoints(object = r_1$t, xy = preds_1$t)
plot(r_1$t)

# offices
r_1$o <- raster(e_1, ncols = 100, nrows = 100) 
proj4string(r_1$o) <- CRS("+proj=longlat +datum=WGS84")
r_1$o <- distanceFromPoints(object = r_1$o, xy = preds_1$o)
plot(r_1$o)

# universities
r_1$u <- raster(e_1, ncols = 100, nrows = 100) 
proj4string(r_1$u) <- CRS("+proj=longlat +datum=WGS84")
r_1$u <- distanceFromPoints(object = r_1$u, xy = preds_1$u)
plot(r_1$u)

# create raster stack 
r_l_1 <- as.list(r_1)
r_stack_1 <- stack(r_l_1)
writeRaster(r_stack,file.path(dirs$main_dir,"distanceStack_1.tif"),overwrite=T)


#----------------------------- 5.2 Species distribution model -----------------------------#

# trip starts
p_s1 <- list()
for(i in 1:length(pres_s)){
  d1 <- sdmData(formula = freq ~ c1+o+t+u, train = pres_s[[i]], predictors = r_stack)
  m1 <- sdm(freq~., data = d1, methods = c('rf','svm')) # fit species distribution model
  p_s1[[i]] <- predict(m1, newdata = r_stack_1, overwrite = T) # predict
}

#----------------------------------- 5.3 Visualization -----------------------------------#
# first plots
plot(p_s1[[1]]$id_2.sp_1.m_svm) # weekdays
plot(p_s1[[2]]$id_2.sp_1.m_svm) # weekends
plot(p_s1[[3]]$id_2.sp_1.m_svm) # weekdays at 8 am
plot(p_s1[[4]]$id_2.sp_1.m_svm) # weekdays at 4 pm
plot(p_s1[[5]]$id_2.sp_1.m_svm) # weekend at 8 am
plot(p_s1[[6]]$id_2.sp_1.m_svm) # weekend at 4 am


# base map (needs to be replaced with your functions etc. :))
e_l1 <- as.list(e_1) # save station extent as list to use it as input for the openmap 
map_1 <- openmap(upperLeft = e_l1[c(4,1)], lowerRight = e_l1[c(3,2)], type = "opencyclemap", mergeTiles = T)
map_1 <- openproj(map_1,projection = mycrs)
plot(map_1,add=F)
plot(p_s1[[1]]$id_2.sp_1.m_svm,alpha=0.7,add=T) # weekdays
plot(map_1,add=F)
plot(p_s1[[2]]$id_2.sp_1.m_svm,alpha=0.5,add=T) # weekends
plot(map_1,add=F)
plot(p_s1[[3]]$id_2.sp_1.m_svm,alpha=0.5,add=T)
plot(map_1,add=F)
plot(p_s1[[4]]$id_2.sp_1.m_svm,alpha=0.5,add=T)
plot(map_1,add=F)
plot(p_s1[[5]]$id_2.sp_1.m_svm,alpha=0.5,add=T)

