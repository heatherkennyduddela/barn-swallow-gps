
# calculate kernel density estimations using adehabitatHR

# Heather Kenny-Duddela
# Sep 11, 2023

# libraries
library(tidyverse)
library(ggbreak)
library(stringr)
library(adehabitatHR)
# sp package used by adehabitatHR was retired in 2023, need sf instead
library(sf)
library(ggmap) # for mapping
library(RColorBrewer) # for color palletes
# use ggspatial instead for scale bar and north arrow
# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
library(ggspatial)

# set working directory
setwd("~/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/process-spatial-data")

 
# Read in GPS points
allsites3 <- read.csv("input-data/GPS points all sites.csv")
allsites3 <- allsites3[,-1] # remove extra index column

# read in location of barns
barns <- read.csv("input-data/Barns where tags were deployed 2021.csv")

# new table with only bird and lat long (required format for KDE), keep days with few points
allsites_gps <- allsites3[,c(17, 9, 10)]

#indicate which columns give coordinates
coordinates(allsites_gps) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(allsites_gps) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# convert to sf object, and specify lat long columns
gpsEN <- st_as_sfc(allsites_gps, coords=c("Longitude", "Latitude"))
st_crs(gpsEN)

## Probably don't need this code chunk anymore (11/14/24)

# p4s <- "+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"
# # transform points to the proper projection
# st_transform(gpsEN, crs=p4s)

# convert to "SpatialPoints" to be recognized by adehabitatHR
gpsEN <- as(gpsEN, "Spatial")
# the converstion changed the crs, so reset to correct
gpsEN <- spTransform(allsites_gps, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))



#-------------------------------------------------------------------------------
# Calculate kdes
#-------------------------------------------------------------------------------

# calculate UD for all birds using the href method to calculate the smoothing parameter
all_kde <- kernelUD(gpsEN, h="href")
# visualize results
image(all_kde)

# check h-value smoothing parameter, calculated for each animal
all_kde[[1]]@h #175.3883, bc09

all_kde[[2]]@h #108.2257, bc14

all_kde[[3]]@h #128.4294, cao3

all_kde[[4]]@h #92.25929, co23

all_kde[[5]]@h #112.4564, co27

all_kde[[6]]@h #434.4182, co31

all_kde[[7]]@h #93.91769, mb26

all_kde[[8]]@h #97.5634, mb69

all_kde[[9]]@h #175.9738, sc108

all_kde[[10]]@h #268.9254, sc80

all_kde[[11]]@h #101.1918, st20


# calculate average h-val
havg <- mean(c(175.3883, 108.2257, 128.4294, 92.25929 ,112.4564 ,434.4182, 
               93.91769,
             97.5634, 175.9738, 268.9254 ,101.1918 ))

# calculate median h to reduce emphasis on co31
hmed <- median(c(175.3883, 108.2257, 128.4294, 92.25929 ,112.4564 ,
                 434.4182, 93.91769,
                 97.5634, 175.9738, 268.9254 ,101.1918 ))

kde_med <- kernelUD(gpsEN, h=hmed)
# image(kde_med)

kde_med2 <- kernelUD(gpsEN, h=hmed, grid=500, same4all=T)
# image(kde_med2)

kde_med3 <- kernelUD(gpsEN, h=hmed, grid=200, same4all=T, extent=0.02)
# image(kde_med3)

kde_med4 <- kernelUD(gpsEN, h=hmed, grid=500, same4all=F)
# image(kde_med4)


## calculate 90% utilization area

# this calculates the 90% home range in vector form (rather than raster)
range.indiv <- getverticeshr(all_kde, percent=90)



# Try plotting using ggspatial and sf------------------------------------------



## Make spatial polygon maps for each bird overall


# You will need needed to register with Stadia Map and get an API key
register_stadiamaps("[insert API key here]", write=F)


basemap <- get_stadiamap(bbox = c(left = min(allsites_gps@coords[,1])-0.01, 
                                bottom = min(allsites_gps@coords[,2])-0.01, 
                                right = max(allsites_gps@coords[,1])+0.01, 
                                top = max(allsites_gps@coords[,2])+0.01),
                         maptype = "stamen_terrain_background",
                         zoom=12)

base2 <- get_stadiamap(bbox = c(left = min(allsites_gps@coords[,1])-0.01, 
                                bottom = min(allsites_gps@coords[,2])-0.01, 
                                right = max(allsites_gps@coords[,1])+0.01, 
                                top = max(allsites_gps@coords[,2])+0.01),
                       maptype = "stamen_terrain",
                       zoom=12)

# problem is that map is in long lat and KDE polygons are in utm
# convert polygons to long lat

range.indiv.longlat <- spTransform(range.indiv, CRS('+proj=longlat +datum=WGS84 +no_defs'))


########### Calculate 50% KDEs for all ################

range.indiv.50 <- getverticeshr(all_kde, percent=50)
plot(range.indiv.50)

range.indiv.50.longlat <- spTransform(range.indiv.50, CRS('+proj=longlat +datum=WGS84 +no_defs'))

# custom colors
my.color = c("#FFCC33", "#006600","#FF0000","#660033","#FF6600","#0000CC","#330066","#993300","#CC00FF","#666666","#00CCCC")

# convert polygons and points to sf
range.indiv.longlat.sf <- sf::st_as_sf(range.indiv.longlat)
range.indiv.50.sf <- sf::st_as_sf(range.indiv.50.longlat)
allsites_gps.sf <- sf::st_as_sf(allsites_gps)
# remove Hepp and Mary Ann barns, only keep original nesting sites
barns.sf <- sf::st_as_sf(subset(barns, barns$Site!="Hepp" &
                                  barns$Site!="Mary Ann"), 
                         coords=c("Long", "Lat"))
# add crs for barns
sf::st_crs(range.indiv.50.sf)
st_crs(barns.sf) <- 4326

# Manuscript Figure 2
# both 90% and 50% together
# note that the key to get ggmap and sf to work together is to set
# inherit.aes = F for all of the sf layers
ggmap(base2, extent="device") +
  geom_sf(data=range.indiv.longlat.sf, 
          aes(colour = id, fill = id), alpha = 0.3, 
          inherit.aes = F) +
  geom_sf(data=range.indiv.50.sf, 
          aes(color=id, fill=id), alpha=0.1, linetype="dashed", 
          inherit.aes = F) +
  geom_sf(data=allsites_gps.sf, aes(color=nest), shape=1,
          inherit.aes = F) +
  geom_sf(data=barns.sf, fill="black", color="gray", shape=24, size=2,
          inherit.aes = F) +
  annotation_scale(location="br", width_hint=0.3) +
  scale_fill_manual(labels=c("Blue Cloud 09", "Blue Cloud 14", "Cathy's 03" , 
                             "Cook 23", "Cook 27", "Cook 31", "Make Believe 26",
                             "Make Believe 69", "Schaaps 108", "Schaaps 80", 
                             "Struthers 20"), values=my.color) +
  scale_color_manual(labels=c("Blue Cloud 09", "Blue Cloud 14", "Cathy's 03" ,
                              "Cook 23", "Cook 27", "Cook 31", 
                              "Make Believe 26", "Make Believe 69", 
                              "Schaaps 108", "Schaaps 80", "Struthers 20"), 
                     values=my.color) + 
  labs(fill="Bird ID") + guides(color="none") +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Movement areas and GPS points for \ntracked females 2021") +
  scale_x_continuous(breaks=c(-105.26, -105.20, -105.14), limits=c(-105.262, -105.138)) +
  theme_classic()

ggsave("map of kde areas scalebar .png", h=7.5, w=6)
ggsave("map of kde areas scalebar arrow.pdf", h=7.5, w=6)


################################################################################
# Note: this section of code no longer runs because sp has been depreciated
# Plots can be created using ggspatial and sf as done above for Figure 2,
# however, I have not implemented the updated code here. These plots were made
# when sp was still functional. 
################################################################################



##### zoomed in map for middle sites

# get extent to points
middle.sites <- subset(allsites3, allsites3$nest=="bc09" |
                         allsites3$nest=="bc14" |
                         allsites3$nest=="st20" |
                         allsites3$nest=="mb26" |
                         allsites3$nest=="mb69")

middle.sites_gps <- middle.sites[,c(17, 9, 10)]

#indicate which columns give coordinates
coordinates(middle.sites_gps) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(middle.sites_gps) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
middle.gpsEN <- spTransform(middle.sites_gps, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))

middle.base <- get_stadiamap(bbox = c(left = min(middle.sites_gps@coords[,1])-0.005, 
                                  bottom = min(middle.sites_gps@coords[,2])-0.005, 
                                  right = max(middle.sites_gps@coords[,1])+0.005, 
                                  top = max(middle.sites_gps@coords[,2])+0.005),
                         maptype = "stamen_terrain",
                         zoom=12)

my.color2 = c("#FFCC33", "#006600","#330066","#993300","#00CCCC")

ggmap(middle.base) + 
  geom_polygon(data = fortify(range.indiv.longlat[c(1,2,7,8,11),]), 
               aes(x=long, y=lat, colour = id, fill = id, group=group),
               alpha = 0.3) +
  geom_polygon(data=range.indiv.50.longlat[c(1,2,7,8,11),],
               aes(x=long, y=lat, color=id, group=group, fill=id),
               alpha=0.1, linetype="dashed") +
  scale_fill_manual(labels=c("Blue Cloud 09", "Blue Cloud 14", 
                             "Make Believe 26", "Make Believe 69", 
                             "Struthers 20"), values=my.color2) +
  scale_color_manual(labels=c("Blue Cloud 09", "Blue Cloud 14",  
                              "Make Believe 26", "Make Believe 69", 
                              "Struthers 20"), values=my.color2) + 
  labs(fill="Bird ID") + guides(colour="none") +
  geom_point(data = middle.sites, 
             aes(x = Longitude, y = Latitude, colour = nest),
             shape=1) +
  geom_point(data=barns[c(2,3,6),], aes(x=Long, y=Lat), fill="black", 
             color="gray", shape=24, size=2) +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Movement areas and GPS points for tracked females, middle sites") +
  scalebar(x.min=-105.20, x.max=-105.178, y.min=40.12, y.max=40.125, dist=1, 
           location="bottomright", dist_unit="km", transform=T, model="WGS84", 
           height=0.2, st.dist=0.3)

# Figure S5 b
ggsave("map of kde areas for middle sites.png", h=6, w=8.5)
ggsave("map of kde areas for middle sites.pdf", h=6, w=8.5)


### zoomed in map of Cooks -----------------------------------------------------

# get extent to points
cooks.detail <- subset(allsites3, allsites3$nest=="co23" |
                         allsites3$nest=="co27")

cooks.detail_gps <- cooks.detail[,c(17, 9, 10)]

#indicate which columns give coordinates
coordinates(cooks.detail_gps) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(cooks.detail_gps) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
cooks.gpsEN <- spTransform(cooks.detail_gps, 
                           CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))

cooks.base <- get_stadiamap(bbox = c(left = min(cooks.detail_gps@coords[,1])-0.005, 
                                      bottom = min(cooks.detail_gps@coords[,2])-0.005, 
                                      right = max(cooks.detail_gps@coords[,1])+0.005, 
                                      top = max(cooks.detail_gps@coords[,2])+0.005),
                             maptype = "stamen_terrain",
                             zoom=15)

my.color3 = c("#660033","#FF6600")

ggmap(cooks.base) + 
  geom_polygon(data = fortify(range.indiv.longlat[c(4,5),]), 
               aes(x=long, y=lat, colour = id, fill = id, group=group),
               alpha = 0.3) +
  geom_polygon(data=range.indiv.50.longlat[c(4,5),],
               aes(x=long, y=lat, color=id, group=group, fill=id),
               alpha=0.1, linetype="dashed") +
  scale_fill_manual(labels=c("Cook 23", "Cook 27"), values=my.color3) +
  scale_color_manual(labels=c("Cook 23", "Cook 27"), values=my.color3) + 
  labs(fill="Bird ID") + guides(colour="none") +
  geom_point(data = cooks.detail, 
             aes(x = Longitude, y = Latitude, colour = nest),
             shape=1) +
  geom_point(data=barns[c(4),], aes(x=Long, y=Lat), fill="black", color="gray", 
             shape=24, size=2) +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Movement areas and GPS points for tracked females, Cooks detail") +
  scalebar(x.min=-105.185, x.max=-105.168, y.min=40.176, y.max=40.178, dist=1, 
           location="bottomright", dist_unit="km", transform=T, model="WGS84", 
           height=0.3, st.dist=0.3)

# Figure S5 a
ggsave("map of kde areas for Cooks detail.png", h=6, w=8.5)
ggsave("map of kde areas for Cooks detail.pdf", h=6, w=8.5)


### Zoomed in map for Schaaps---------------------------------------------------

# get extent to points
schaaps.detail <- subset(allsites3, allsites3$nest=="sc80" |
                         allsites3$nest=="sc108")

schaaps.detail_gps <- schaaps.detail[,c(17, 9, 10)]

#indicate which columns give coordinates
coordinates(schaaps.detail_gps) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(schaaps.detail_gps) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
schaaps.gpsEN <- spTransform(schaaps.detail_gps, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))

schaaps.base <- get_stadiamap(bbox = c(left = min(schaaps.detail_gps@coords[,1])-0.005, 
                                     bottom = min(schaaps.detail_gps@coords[,2])-0.005, 
                                     right = max(schaaps.detail_gps@coords[,1])+0.005, 
                                     top = max(schaaps.detail_gps@coords[,2])+0.005),
                            maptype = "stamen_terrain",
                            zoom=13)

my.color4 = c("#CC00FF","#666666")

ggmap(schaaps.base) + 
  geom_polygon(data = fortify(range.indiv.longlat[c(9, 10),]), 
               aes(x=long, y=lat, colour = id, fill = id, group=group),
               alpha = 0.3) +
  geom_polygon(data=range.indiv.50.longlat[c(9, 10),],
               aes(x=long, y=lat, color=id, group=group, fill=id),
               alpha=0.1, linetype="dashed") +
  scale_fill_manual(labels=c("Schaaps 108", "Schaaps 80"), values=my.color4) +
  scale_color_manual(labels=c("Schaaps 108", "Schaaps 80"), values=my.color4) + 
  labs(fill="Bird ID") + guides(colour="none") +
  geom_point(data = schaaps.detail, 
             aes(x = Longitude, y = Latitude, colour = nest),
             shape=1) +
  geom_point(data=barns[c(1),], aes(x=Long, y=Lat), fill="black", color="gray", shape=24, size=2) +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Movement areas and GPS points for \ntracked females, Schaaps detail") +
  scalebar(x.min=-105.265, x.max=-105.245, y.min=40.065, y.max=40.08, dist=1, 
           location="bottomright", dist_unit="km", transform=T, model="WGS84", 
           height=0.4, st.dist=0.05, st.size=3)

# Figure S6
ggsave("map of kde areas for Schaaps detail.png", h=6, w=5)
ggsave("map of kde areas for Schaaps detail.pdf", h=6, w=5)




# view areas from the vector mode
as.data.frame(range.indiv)


# compare areas to those calculated with the fixed bandwidth
range.fixed <- getverticeshr(kde_med)
as.data.frame(range.fixed)

################################################################################
# This section of code uses ggspatial and sf, rather than sp
# Plot polygons of daily movement to show variation within females across days
################################################################################

#-------------------------------------------------------------------------------
# Plot daily polygons
#-------------------------------------------------------------------------------

# BlueCloud-09 -----------------------------------------------------------------

## bc09
bc09_nest <- filter(allsites3, grepl("bc09", trackday))
bc09_nest.poly <- bc09_nest[,c(9, 10, 21)]
bc09_barn <- filter(barns, grepl("BlueCloud", Site))

#indicate which columns give coordinates
coordinates(bc09_nest.poly) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(bc09_nest.poly) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
nest.polyEN_bc09 <- spTransform(bc09_nest.poly, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))
# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
bc09_nest_EN100 <- mcp(bc09_nest.poly, percent = 100, unout = "m2")


# get base map
basemap_bc09 <- get_stadiamap(bbox = c(left = min(bc09_nest.poly@coords[,1])-0.005, 
                                       bottom = min(bc09_nest.poly@coords[,2])-0.005, 
                                       right = max(bc09_nest.poly@coords[,1])+0.005, 
                                       top = max(bc09_nest.poly@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
bc09_EN100.sf <- sf::st_as_sf(bc09_nest_EN100)
bc09.points.sf <- sf::st_as_sf(bc09_nest[,c(9,10,21)], 
                               coords=c("Longitude","Latitude"))
bc09.barn.sf <- st_as_sf(bc09_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(bc09.points.sf) <- 4326
st_crs(bc09.barn.sf) <- 4326

#Create the map with foraging range polygons! 
mymap.bc09 <- ggmap(basemap_bc09) + 
  geom_sf(data = bc09_EN100.sf,
               aes(colour = id, fill = id),
               alpha = 0.3, inherit.aes = F) +
  geom_sf(data = bc09.points.sf,
             aes(color = trackday), inherit.aes = F) +
  geom_sf(data=bc09.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for BlueCloud-09") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  annotation_north_arrow(location="br", width=unit(1, "cm"), 
                         pad_y=unit(1,"cm")) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

mymap.bc09

# Figure 1b
ggsave("generated-datatables/bc09_daily_polygons.png", h=4, w=5)


# BlueCloud-14 -----------------------------------------------------------------

## bc14
bc14_nest <- filter(allsites3, grepl("bc14", trackday))
bc14_nest.poly <- bc14_nest[,c(9, 10, 21)]
bc14_barn <- filter(barns, grepl("BlueCloud", Site))

#indicate which columns give coordinates
coordinates(bc14_nest.poly) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(bc14_nest.poly) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
nest.polyEN_bc14 <- spTransform(bc14_nest.poly, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))
# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
bc14_nest_EN100 <- mcp(bc14_nest.poly, percent = 100, unout = "m2")


# get base map
basemap_bc14 <- get_stadiamap(bbox = c(left = min(bc14_nest.poly@coords[,1])-0.005, 
                                       bottom = min(bc14_nest.poly@coords[,2])-0.005, 
                                       right = max(bc14_nest.poly@coords[,1])+0.005, 
                                       top = max(bc14_nest.poly@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
bc14_EN100.sf <- sf::st_as_sf(bc14_nest_EN100)
bc14.points.sf <- sf::st_as_sf(bc14_nest[,c(9,10,21)], 
                               coords=c("Longitude","Latitude"))
bc14.barn.sf <- st_as_sf(bc14_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(bc14.points.sf) <- 4326
st_crs(bc14.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_bc14) + 
  geom_sf(data = bc14_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = bc14.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=bc14.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for BlueCloud-14") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure 1a
ggsave("generated-datatables/bc14_daily_polygons.png", h=4, w=5)


# Cathys-03 -----------------------------------------------------------------

## ca03
ca03_nest <- filter(allsites3, grepl("ca03", trackday))
ca03_nest.poly <- ca03_nest[,c(9, 10, 21)]
ca03_barn <- filter(barns, grepl("Cathys", Site))

#indicate which columns give coordinates
coordinates(ca03_nest.poly) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(ca03_nest.poly) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
nest.polyEN_ca03 <- spTransform(ca03_nest.poly, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))
# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
ca03_nest_EN100 <- mcp(ca03_nest.poly, percent = 100, unout = "m2")


# get base map
basemap_ca03 <- get_stadiamap(bbox = c(left = min(ca03_nest.poly@coords[,1])-0.005, 
                                       bottom = min(ca03_nest.poly@coords[,2])-0.005, 
                                       right = max(ca03_nest.poly@coords[,1])+0.005, 
                                       top = max(ca03_nest.poly@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
ca03_EN100.sf <- sf::st_as_sf(ca03_nest_EN100)
ca03.points.sf <- sf::st_as_sf(ca03_nest[,c(9,10,21)], 
                               coords=c("Longitude","Latitude"))
ca03.barn.sf <- st_as_sf(ca03_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(ca03.points.sf) <- 4326
st_crs(ca03.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_ca03) + 
  geom_sf(data = ca03_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = ca03.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=ca03.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for Cathys-03") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure S4a
ggsave("generated-datatables/ca03_daily_polygons.png", h=4, w=5)


# Cooks-23 -----------------------------------------------------------------

## co23
co23_nest <- filter(allsites3, grepl("co23", trackday))
co23_nest.poly <- co23_nest[,c(9, 10, 21)]
# the last tracking day has only 3 points so remove those 
co23_nest.poly2 <- subset(co23_nest.poly, co23_nest.poly$trackday!="co23-06-11")
co23_barn <- filter(barns, grepl("Cooks", Site))

#indicate which columns give coordinates
coordinates(co23_nest.poly2) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(co23_nest.poly2) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
co23_nest_EN100 <- mcp(co23_nest.poly2, percent = 100, unout = "m2")


# get base map
basemap_co23 <- get_stadiamap(bbox = c(left = min(co23_nest.poly2@coords[,1])-0.005, 
                                       bottom = min(co23_nest.poly2@coords[,2])-0.005, 
                                       right = max(co23_nest.poly2@coords[,1])+0.005, 
                                       top = max(co23_nest.poly2@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
co23_EN100.sf <- sf::st_as_sf(co23_nest_EN100)
co23.points.sf <- sf::st_as_sf(subset(co23_nest[,c(9,10,21)], 
                                      co23_nest$trackday!="co23-06-11"), 
                               coords=c("Longitude","Latitude"))
co23.barn.sf <- st_as_sf(co23_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(co23.points.sf) <- 4326
st_crs(co23.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_co23) + 
  geom_sf(data = co23_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = co23.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=co23.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for Cooks-23") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure 1c
ggsave("generated-datatables/co23_daily_polygons.png", h=4, w=5)


# Cooks-31 -----------------------------------------------------------------

## co31
co31_nest <- filter(allsites3, grepl("co31", trackday))
co31_nest.poly <- co31_nest[,c(9, 10, 21)]
# the last tracking day has only 3 points so remove those 
co31_nest.poly2 <- subset(co31_nest.poly, co31_nest.poly$trackday!="co31-06-22")
co31_barn <- filter(barns, grepl("Cooks", Site))

#indicate which columns give coordinates
coordinates(co31_nest.poly2) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(co31_nest.poly2) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
co31_nest_EN100 <- mcp(co31_nest.poly2, percent = 100, unout = "m2")


# get base map
basemap_co31 <- get_stadiamap(bbox = c(left = min(co31_nest.poly2@coords[,1])-0.005, 
                                       bottom = min(co31_nest.poly2@coords[,2])-0.005, 
                                       right = max(co31_nest.poly2@coords[,1])+0.005, 
                                       top = max(co31_nest.poly2@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
co31_EN100.sf <- sf::st_as_sf(co31_nest_EN100)
co31.points.sf <- sf::st_as_sf(subset(co31_nest[,c(9,10,21)], 
                                      co31_nest$trackday!="co31-06-22"), 
                               coords=c("Longitude","Latitude"))
co31.barn.sf <- st_as_sf(co31_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(co31.points.sf) <- 4326
st_crs(co31.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_co31) + 
  geom_sf(data = co31_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = co31.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=co31.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for Cooks-31") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  annotation_north_arrow(location="br", width=unit(1, "cm"), 
                         pad_y=unit(1,"cm")) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure 1d
ggsave("generated-datatables/co31_daily_polygons.png", h=4, w=5)


# MakeBelieve-26 ---------------------------------------------------------------

## mb26
mb26_nest <- filter(allsites3, grepl("mb26", trackday))
mb26_nest.poly <- mb26_nest[,c(9, 10, 21)]
mb26_barn <- filter(barns, grepl("MakeBelieve", Site))

#indicate which columns give coordinates
coordinates(mb26_nest.poly) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(mb26_nest.poly) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
mb26_nest_EN100 <- mcp(mb26_nest.poly, percent = 100, unout = "m2")


# get base map
basemap_mb26 <- get_stadiamap(bbox = c(left = min(mb26_nest.poly@coords[,1])-0.005, 
                                       bottom = min(mb26_nest.poly@coords[,2])-0.005, 
                                       right = max(mb26_nest.poly@coords[,1])+0.005, 
                                       top = max(mb26_nest.poly@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
mb26_EN100.sf <- sf::st_as_sf(mb26_nest_EN100)
mb26.points.sf <- sf::st_as_sf(mb26_nest[,c(9,10,21)], 
                               coords=c("Longitude","Latitude"))
mb26.barn.sf <- st_as_sf(mb26_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(mb26.points.sf) <- 4326
st_crs(mb26.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_mb26) + 
  geom_sf(data = mb26_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = mb26.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=mb26.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for MakeBelieve-26") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure S4d
ggsave("generated-datatables/mb26_daily_polygons.png", h=4, w=5)


# MakeBelieve-69 ---------------------------------------------------------------

## mb69
mb69_nest <- filter(allsites3, grepl("mb69", trackday))
mb69_nest.poly <- mb69_nest[,c(9, 10, 21)]
mb69_barn <- filter(barns, grepl("MakeBelieve", Site))

#indicate which columns give coordinates
coordinates(mb69_nest.poly) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(mb69_nest.poly) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
mb69_nest_EN100 <- mcp(mb69_nest.poly, percent = 100, unout = "m2")


# get base map
basemap_mb69 <- get_stadiamap(bbox = c(left = min(mb69_nest.poly@coords[,1])-0.005, 
                                       bottom = min(mb69_nest.poly@coords[,2])-0.005, 
                                       right = max(mb69_nest.poly@coords[,1])+0.005, 
                                       top = max(mb69_nest.poly@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
mb69_EN100.sf <- sf::st_as_sf(mb69_nest_EN100)
mb69.points.sf <- sf::st_as_sf(mb69_nest[,c(9,10,21)], 
                               coords=c("Longitude","Latitude"))
mb69.barn.sf <- st_as_sf(mb69_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(mb69.points.sf) <- 4326
st_crs(mb69.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_mb69) + 
  geom_sf(data = mb69_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = mb69.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=mb69.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for MakeBelieve-69") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure S4c
ggsave("generated-datatables/mb69_daily_polygons.png", h=4, w=5)


# Schaaps 108 ---------------------------------------------------------------

## sc108
sc108_nest <- filter(allsites3, grepl("sc108", trackday))
sc108_nest.poly <- sc108_nest[,c(9, 10, 21)]
sc108_barn <- filter(barns, grepl("Schaaps", Site))

#indicate which columns give coordinates
coordinates(sc108_nest.poly) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(sc108_nest.poly) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
sc108_nest_EN100 <- mcp(sc108_nest.poly, percent = 100, unout = "m2")


# get base map
basemap_sc108 <- get_stadiamap(bbox = c(left = min(sc108_nest.poly@coords[,1])-0.005, 
                                       bottom = min(sc108_nest.poly@coords[,2])-0.005, 
                                       right = max(sc108_nest.poly@coords[,1])+0.005, 
                                       top = max(sc108_nest.poly@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
sc108_EN100.sf <- sf::st_as_sf(sc108_nest_EN100)
sc108.points.sf <- sf::st_as_sf(sc108_nest[,c(9,10,21)], 
                               coords=c("Longitude","Latitude"))
sc108.barn.sf <- st_as_sf(sc108_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(sc108.points.sf) <- 4326
st_crs(sc108.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_sc108) + 
  geom_sf(data = sc108_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = sc108.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=sc108.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for Schaaps-108") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure S4e
ggsave("generated-datatables/sc108_daily_polygons.png", h=4, w=5)


# Schaaps 80 ---------------------------------------------------------------

## sc80
sc80_nest <- filter(allsites3, grepl("sc80", trackday))
sc80_nest.poly <- sc80_nest[,c(9, 10, 21)]
# last trackday has only 4 points, so remove these
sc80_nest.poly2 <- subset(sc80_nest.poly, sc80_nest.poly$trackday!="sc80-06-16")
sc80_barn <- filter(barns, grepl("Schaaps", Site))

#indicate which columns give coordinates
coordinates(sc80_nest.poly2) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(sc80_nest.poly2) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
sc80_nest_EN100 <- mcp(sc80_nest.poly2, percent = 100, unout = "m2")


# get base map
basemap_sc80 <- get_stadiamap(bbox = c(left = min(sc80_nest.poly2@coords[,1])-0.005, 
                                       bottom = min(sc80_nest.poly2@coords[,2])-0.005, 
                                       right = max(sc80_nest.poly2@coords[,1])+0.005, 
                                       top = max(sc80_nest.poly2@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
sc80_EN100.sf <- sf::st_as_sf(sc80_nest_EN100)
sc80.points.sf <- sf::st_as_sf(subset(sc80_nest[,c(9,10,21)], 
                                      sc80_nest$trackday!="sc80-06-16"), 
                               coords=c("Longitude","Latitude"))
sc80.barn.sf <- st_as_sf(sc80_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(sc80.points.sf) <- 4326
st_crs(sc80.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_sc80) + 
  geom_sf(data = sc80_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = sc80.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=sc80.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for Schaaps-80") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure S4f
ggsave("generated-datatables/sc80_daily_polygons.png", h=4, w=5)


# Struthers 20 ---------------------------------------------------------------

## st20
st20_nest <- filter(allsites3, grepl("st20", trackday))
st20_nest.poly <- st20_nest[,c(9, 10, 21)]
# last trackday has only 3 points, so remove these
st20_nest.poly2 <- subset(st20_nest.poly, st20_nest.poly$trackday!="st20-06-11")
st20_barn <- filter(barns, grepl("Struthers", Site))

#indicate which columns give coordinates
coordinates(st20_nest.poly2) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(st20_nest.poly2) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")

# create the polygons
#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
st20_nest_EN100 <- mcp(st20_nest.poly2, percent = 100, unout = "m2")


# get base map
basemap_st20 <- get_stadiamap(bbox = c(left = min(st20_nest.poly2@coords[,1])-0.005, 
                                       bottom = min(st20_nest.poly2@coords[,2])-0.005, 
                                       right = max(st20_nest.poly2@coords[,1])+0.005, 
                                       top = max(st20_nest.poly2@coords[,2])+0.005), 
                              maptype="stamen_terrain",
                              zoom = 12)

# convert polygons and points to sf
st20_EN100.sf <- sf::st_as_sf(st20_nest_EN100)
st20.points.sf <- sf::st_as_sf(subset(st20_nest[,c(9,10,21)], 
                                      st20_nest$trackday!="st20-06-11"), 
                               coords=c("Longitude","Latitude"))
st20.barn.sf <- st_as_sf(st20_barn, coords=c("Long","Lat"))
# add crs for points
st_crs(st20.points.sf) <- 4326
st_crs(st20.barn.sf) <- 4326

#Create the map with foraging range polygons! 
ggmap(basemap_st20) + 
  geom_sf(data = st20_EN100.sf,
          aes(colour = id, fill = id),
          alpha = 0.3, inherit.aes = F) +
  geom_sf(data = st20.points.sf,
          aes(color = trackday), inherit.aes = F) +
  geom_sf(data=st20.barn.sf, fill="black", color="gray", shape=24, size=3,
          inherit.aes = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Daily polygons for Struthers-20") +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location="br", width_hint=0.3) +
  guides(color=guide_legend(title="Tracking day"),
         fill=guide_legend(title="Tracking day"))

# Figure S4b
ggsave("generated-datatables/st20_daily_polygons.png", h=4, w=5)


# Caculate daily polygon size for all sites together --------------------------

# new table with only trackday and lat long (required format for MPC), days with few points removed
allsites_good <- subset(allsites3, allsites3$fewpoints==FALSE)
allsites_gps <- allsites_good[,c(21, 9, 10)]

#indicate which columns give coordinates
coordinates(allsites_gps) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(allsites_gps) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
gpsEN <- spTransform(allsites_gps, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))


#Percent defines what percentage of points you want to keep; unout defines what
#units you want the output to be in. I've selected meters squared
allsites_MPC_EN100 <- mcp(gpsEN, percent = 100, unout = "m2")


# put areas into dataframe
areas_trackday <- as.data.frame(cbind(allsites_MPC_EN100$id, allsites_MPC_EN100$area))
colnames(areas_trackday) <- c("trackday","MPC100")
# convert to km2 instead of m2
areas_trackday$km2 <- as.numeric(areas_trackday$MPC100)/1000000

# add site and date info to areas df
areas_trackday$nest <- str_sub(areas_trackday$trackday,1,4)
areas_trackday$date <- as.Date(paste("2021-",str_sub(areas_trackday$trackday, 6,10)))
# fix info for sc108
areas_trackday$date[69:75] <- as.Date(
  paste("2021-",str_sub(areas_trackday$trackday[69:75], 7,11)))
areas_trackday$nest[69:75] <- str_sub(areas_trackday$trackday[69:75],1,5)

# all sites together by date
ggplot(areas_trackday, aes(x=date, y=km2, color=nest)) + geom_point() + 
  geom_line() +
  scale_color_manual(labels=c("Blue Cloud 09", "Blue Cloud 14", "Cathy's 03" ,
                                          "Cook 23", "Cook 27", "Cook 31", 
                                          "Make Believe 26", "Make Believe 69", 
                                          "Schaaps 108", "Schaaps 80", "Struthers 20"), 
                                 values=my.color) +
  scale_y_break(c(2.5,7), scales=0.25, ticklabels = c(7, 8)) +
  ylim(0,8) + xlab("Day of tracking") +
  ylab("Minimum convex polygon area\n(square km)") +
  ggtitle("Daily movement range over time for each female") +
  guides(color=guide_legend(title="Bird ID"))

# Figure 1E
ggsave("generated-datatables/movement range over time all birds.png", h=4, w=10)


#-------------------------------------------------------------------------------
# Calculate KDE at multiple probability levels
#-------------------------------------------------------------------------------


### to calculate the range at multiple probability levels, we need to use raster mode

vud.indiv <- getvolumeUD(all_kde)

kernel.probs <- kernel.area(vud.indiv, percent=seq(50, 95, by=5), unout="ha")

vud.fixed <- getvolumeUD(kde_med)

kernel.probs.fixed <- kernel.area(vud.fixed, percent=seq(50, 95, by=5))

plot(kernel.probs)

# try converting to dataframe to compare different calculations

kernel.probs.t <- t(kernel.probs)
kernel.probs.df <- as.data.frame(kernel.probs.t[,c(9,10)])
colnames(kernel.probs.df) <- c("indiv.hval.90", "indiv.hval.95")

kernel.fixed.t <- t(kernel.probs.fixed)
kernel.fixed.df <- as.data.frame(kernel.fixed.t[,c(9,10)])
colnames(kernel.fixed.df) <- c("fixed.hval.90", "fixed.hval.95")

both.df <- cbind(kernel.probs.df, kernel.fixed.df)
both.df$bird <- rownames(both.df)

# reset graphing settings
dev.off()

library(ggplot2)

ggplot(both.df, aes(x=indiv.hval.90, y=fixed.hval.90)) + geom_line() + 
  geom_point(aes(color=bird), size=3) + 
  geom_abline(slope=1, intercept=0, linetype="dashed")+
  ggtitle("Comparison of 90% ranges") +
  xlab("Area with individual bandwidth") +
  ylab("Area with fixed bandwidth (113)")

ggplot(both.df, aes(x=indiv.hval.95, y=fixed.hval.95)) + geom_line() + 
  geom_point(aes(color=bird), size=3) + 
  geom_abline(slope=1, intercept=0, linetype="dashed")+
  ggtitle("Comparison of 95% ranges")+
  xlab("Area with individual bandwidth") +
  ylab("Area with fixed bandwidth (113)")


# visualize "core" areas (50 and 75 UD)

kernel.core.df <- as.data.frame(kernel.probs.t[,c(1,6)])
colnames(kernel.core.df) <- c("indiv.hval.50", "indiv.hval.75")

kernel.fixed.core <- as.data.frame(kernel.fixed.t[,c(1,6)])
colnames(kernel.fixed.core) <- c("fixed.hval.50", "fixed.hval.75")

both.core <- cbind(kernel.core.df, kernel.fixed.core)
both.core$bird <- rownames(both.core)

ggplot(both.core, aes(x=indiv.hval.50, y=fixed.hval.50)) + geom_line() + 
  geom_point(aes(color=bird), size=3) + 
  geom_abline(slope=1, intercept=0, linetype="dashed")+
  ggtitle("Comparison of 50% ranges") +
  xlab("Area with individual bandwidth") +
  ylab("Area with fixed bandwidth (113)")

ggplot(both.core, aes(x=indiv.hval.75, y=fixed.hval.75)) + geom_line() + 
  geom_point(aes(color=bird), size=3) + 
  geom_abline(slope=1, intercept=0, linetype="dashed")+
  ggtitle("Comparison of 75% ranges") +
  xlab("Area with individual bandwidth") +
  ylab("Area with fixed bandwidth (113)")


# check correlation between 50% UD and 90% UD

areas_50_95 <- cbind(both.df, both.core[,-5])

ggplot(areas_50_95, aes(x=indiv.hval.50, indiv.hval.90)) +
  geom_line() + geom_point(aes(color=bird), size=3) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ggtitle("Comparison of 50% and 90% UD, individual bandwidth"  ) +
  xlab("Area of 50% UD") +
  ylab("Area of 90% UD")

ggplot(areas_50_95, aes(x=fixed.hval.50, fixed.hval.90)) +
  geom_line() + geom_point(aes(color=bird), size=3) +
  geom_abline(slope=1, intercept=50, linetype="dashed") +
  ggtitle("Comparison of 50% and 90% UD, fixed bandwidth"  ) +
  xlab("Area of 50% UD") +
  ylab("Area of 90% UD")


########### Calculate UDs using indiv bandwidths with minimum cutoff #######

bandwidths <- c(175.3883, 108.2257, 128.4294, 92.25929 ,112.4564 ,434.4182, 93.91769,
                97.5634, 175.9738, 268.9254 ,101.1918 )

ggplot(as.data.frame(bandwidths), aes(x=bandwidths)) + geom_histogram()

# use median bandwidth for all birds under the median (112.4564)
# use individual bandwidths for birds above the median

vud.large <- vud.indiv[c(1,3,6,9,10,11)]

vud.small <- vud.fixed[c(2,4,5,7,8)]

# can't run kernel.area function on subset of the vud.indiv because that changes
# the object class to list instead of "estUDm". Instead calculate kernels and 
# then subset the indiv and fixed versions at the end

probs.indiv <- kernel.area(vud.indiv, percent=seq(50, 95, by=5))
probs.large <- probs.indiv[c(1,3,6,9,10,11)]

probs.fixed <- kernel.area(vud.fixed, percent=seq(50, 95, by=5))
probs.small <- probs.fixed[c(2,4,5,7,8)]

probs.all <- cbind(probs.large, probs.small)
probs.all$probabilities <- row.names(probs.all)

# convert probabilitiles table to long format
probs.all.long <- as.data.frame(t(probs.all))
colnames(probs.all.long) <- paste("prob",colnames(probs.all.long),sep="_")
probs.all.long <- probs.all.long[-12,] # remove extra row
probs.all.long$bird <- row.names(probs.all.long)


write.csv(probs.all.long, "KDE kernels with cutoff value.csv", row.names=F)


################################################################################

### Calculate proportion of points within 50m of the barn


# For barn coordinates, note that Hepp and Mary Ann are not primary tagging sites. 
# One female from Cooks moved to nest at Hepp for her replacement clutch. 
# One female from Cathy's moved to Mary Ann for her replacement clutch 


# add barn coordinates to GPS point table
colnames(barns)[1] <- "site"
barn_bird <- left_join(allsites3, barns, by="site")
colnames(barn_bird)[23:24] <- c("barn_lat","barn_long")


### Calculate distances between bird points and the home barn

library(raster)

barn_bird$distance <- pointDistance(barn_bird[, c(10,9)], # bird locations, order long lat
                                    barn_bird[,24:23], # barn locations, order long lat
                                    longlat=T) 


### calculate proportion points within 50m and 100m

barn_bird$close50 <- barn_bird$distance <= 50
barn_bird$close100 <- barn_bird$distance <=100

# save file with raw distances
write.csv(barn_bird, "point.distances.hkd.csv")

# summarize proportion close for each bird
barn_bird_prop <- barn_bird %>%
  group_by(nest) %>%
  summarise(total.points = n(),
            number.close50 = sum(close50),
            number.close100 = sum(close100),
            prop.close50 = number.close50/total.points,
            prop.close100 = number.close100/total.points)

# add proportion close to table with areas
colnames(barn_bird_prop)[1] <- "bird"
probs.all.close <- left_join(probs.all.long, barn_bird_prop, by="bird")


write.csv(probs.all.close, "KDE probs_updated.csv", row.names=F)







