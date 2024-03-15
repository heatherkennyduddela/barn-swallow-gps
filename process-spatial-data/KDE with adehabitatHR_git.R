
# calculate kernel density estimations using adehabitatHR

# Heather Kenny-Duddela
# Sep 11, 2023

# libraries
library(tidyverse)
library(stringr)
library(adehabitatHR)

# read in data
# note that this table is created by the script "Movement range 2021 allsites.R" 
# in the Swift GPS Data Files folder. Copied here for convenience
allsites3 <- read.csv("GPS points all sites.csv")
allsites3 <- allsites3[,-1] # remove extra index column

barns <- read.csv("Barns where tags were deployed 2021.csv")

# new table with only bird and lat long (required format for KDE), keep days with few points
allsites_gps <- allsites3[,c(17, 9, 10)]

#indicate which columns give coordinates
coordinates(allsites_gps) <- c("Longitude", "Latitude")
#Indicate what format the coordinates are in
proj4string(allsites_gps) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84
+datum=WGS84")
# transform to easting and northing
gpsEN <- spTransform(allsites_gps, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))


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
havg <- mean(c(175.3883, 108.2257, 128.4294, 92.25929 ,112.4564 ,434.4182, 93.91769,
             97.5634, 175.9738, 268.9254 ,101.1918 ))

# calculate median h to reduce emphasis on co31
hmed <- median(c(175.3883, 108.2257, 128.4294, 92.25929 ,112.4564 ,434.4182, 93.91769,
                 97.5634, 175.9738, 268.9254 ,101.1918 ))

kde_med <- kernelUD(gpsEN, h=hmed)
image(kde_med)

kde_med2 <- kernelUD(gpsEN, h=hmed, grid=500, same4all=T)
image(kde_med2)

kde_med3 <- kernelUD(gpsEN, h=hmed, grid=200, same4all=T, extent=0.02)
image(kde_med3)

kde_med4 <- kernelUD(gpsEN, h=hmed, grid=500, same4all=F)
image(kde_med4)


## calculate 95% utilization area

# this calculates the 90% home range in vector form (rather than raster)
range.indiv <- getverticeshr(all_kde, percent=90)

plot(range.indiv)

plot(range.indiv, col=1:11)


#### Try plotting areas using ggplot
df <- fortify(range.indiv)
ggplot(df, aes(x=long, y=lat, fill=id, group=group)) +
  geom_polygon(alpha=0.4)

ggplot() +
  geom_polygon(data=df, aes(x=long, y=lat, fill=id, group=group), alpha=0.4)



#Make spatial polygon maps for each bird overall
library(ggmap)

# You will need needed to register with Stadia Map and get an API key
register_stadiamaps("[insert API key here]", write=F)


basemap <- get_stadiamap(bbox = c(left = min(allsites_gps@coords[,1])-0.005, 
                                bottom = min(allsites_gps@coords[,2])-0.005, 
                                right = max(allsites_gps@coords[,1])+0.005, 
                                top = max(allsites_gps@coords[,2])+0.005),
                         maptype = "stamen_terrain_background",
                         zoom=12)

base2 <- get_stadiamap(bbox = c(left = min(allsites_gps@coords[,1])-0.005, 
                                bottom = min(allsites_gps@coords[,2])-0.005, 
                                right = max(allsites_gps@coords[,1])+0.005, 
                                top = max(allsites_gps@coords[,2])+0.005),
                       maptype = "stamen_terrain",
                       zoom=12)

# problem is that map is in long lat and KDE polygons are in utm
# convert polygons to long lat

range.indiv.longlat <- spTransform(range.indiv, CRS('+proj=longlat +datum=WGS84 +no_defs'))

ggplot(fortify(range.indiv.longlat), aes(x=long, y=lat, fill=id, group=group)) +
  geom_polygon(alpha=0.4)


#Turn the spatial data frame of points into just an R data frame for plotting in ggmap
allsites_gpsdf <- data.frame(allsites_gps@coords, 
                              id = allsites_gps@data$nest)

ggplot() +
  geom_polygon(data=fortify(range.indiv.longlat), aes(x=long, y=lat, fill=id, group=group), alpha=0.4) +
  geom_point(data=allsites_gpsdf,aes(x = Longitude, y = Latitude, colour = id),
             shape=1)


#Create the map with foraging range polygons! 
ggmap(basemap) + 
  geom_polygon(data = fortify(range.indiv.longlat), 
               aes(x=long, y=lat, colour = id, fill = id, group=group),
               alpha = 0.3) + 
  geom_point(data = allsites_gpsdf, 
             aes(x = Longitude, y = Latitude, colour = id),
             shape=1) 


########### Try plotting 50% KDEs for all ################

range.indiv.50 <- getverticeshr(all_kde, percent=50)
plot(range.indiv.50)

range.indiv.50.longlat <- spTransform(range.indiv.50, CRS('+proj=longlat +datum=WGS84 +no_defs'))

ggmap(basemap) + 
  geom_polygon(data = fortify(range.indiv.50.longlat), 
               aes(x=long, y=lat, colour = id, fill = id, group=group),
               alpha = 0.3) + 
  geom_point(data = allsites_gpsdf, 
             aes(x = Longitude, y = Latitude, colour = id),
             shape=1) 

library(RColorBrewer)
library(ggsn)

my.color = c("#FFCC33", "#006600","#FF0000","#660033","#FF6600","#0000CC","#330066","#993300","#CC00FF","#666666","#00CCCC")

# both 90 and 50 together
fullmap <- ggmap(base2) + 
  geom_polygon(data = fortify(range.indiv.longlat), 
               aes(x=long, y=lat, colour = id, fill = id, group=group),
               alpha = 0.3) +
  geom_polygon(data=range.indiv.50.longlat,
               aes(x=long, y=lat, color=id, group=group, fill=id),
               alpha=0.1, linetype="dashed") +
  scale_fill_manual(labels=c("Blue Cloud 09", "Blue Cloud 14", "Cathy's 03" , "Cook 23", "Cook 27", "Cook 31", "Make Believe 26", "Make Believe 69", "Schaaps 108", "Schaaps 80", "Struthers 20"), values=my.color) +
  scale_color_manual(labels=c("Blue Cloud 09", "Blue Cloud 14", "Cathy's 03" , "Cook 23", "Cook 27", "Cook 31", "Make Believe 26", "Make Believe 69", "Schaaps 108", "Schaaps 80", "Struthers 20"), values=my.color) + 
  labs(fill="Bird ID") + guides(color="none") +
  geom_point(data = allsites_gpsdf, 
             aes(x = Longitude, y = Latitude, colour = id),
             shape=1) +
  geom_point(data=barns, aes(x=Long, y=Lat), fill="black", color="gray", shape=24, size=2) +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Movement areas and GPS points for \ntracked females 2021") +
  scalebar(x.min=-105.15, x.max=-105.19, y.min=40.08, y.max=40.09, dist=2, 
           location="bottomright", dist_unit="km", transform=T, model="WGS84", 
           height=0.3, st.dist=0.5)
fullmap

map.with.arrow <- north2(fullmap, x=0.25, y=0.8, symbol=12, scale=0.07)


ggsave("map of kde areas scalebar arrow.png", h=7.5, w=6)
ggsave("map of kde areas scalebar arrow.pdf", h=7.5, w=6)


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
  scale_fill_manual(labels=c("Blue Cloud 09", "Blue Cloud 14", "Make Believe 26", "Make Believe 69", "Struthers 20"), values=my.color2) +
  scale_color_manual(labels=c("Blue Cloud 09", "Blue Cloud 14",  "Make Believe 26", "Make Believe 69", "Struthers 20"), values=my.color2) + 
  labs(fill="Bird ID") + guides(colour="none") +
  geom_point(data = middle.sites, 
             aes(x = Longitude, y = Latitude, colour = nest),
             shape=1) +
  geom_point(data=barns[c(2,3,6),], aes(x=Long, y=Lat), fill="black", color="gray", shape=24, size=2) +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Movement areas and GPS points for tracked females, middle sites") +
  scalebar(x.min=-105.20, x.max=-105.178, y.min=40.12, y.max=40.125, dist=1, 
           location="bottomright", dist_unit="km", transform=T, model="WGS84", 
           height=0.2, st.dist=0.3)

ggsave("map of kde areas for middle sites.png", h=6, w=8.5)
ggsave("map of kde areas for middle sites.pdf", h=6, w=8.5)


### zoomed in map of Cooks

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
cooks.gpsEN <- spTransform(cooks.detail_gps, CRS("+init=epsg:32613 +proj=utm +zone=13 +units=m ellps=WGS84"))

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
  geom_point(data=barns[c(4),], aes(x=Long, y=Lat), fill="black", color="gray", shape=24, size=2) +
  ylab("Latitude") + xlab("Longitude") +
  ggtitle("Movement areas and GPS points for tracked females, Cooks detail") +
  scalebar(x.min=-105.185, x.max=-105.168, y.min=40.176, y.max=40.178, dist=1, 
           location="bottomright", dist_unit="km", transform=T, model="WGS84", 
           height=0.3, st.dist=0.3)

ggsave("map of kde areas for Cooks detail.png", h=6, w=8.5)
ggsave("map of kde areas for Cooks detail.pdf", h=6, w=8.5)


### Zoomed in map for Schaaps

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

ggsave("map of kde areas for Schaaps detail.png", h=6, w=5)
ggsave("map of kde areas for Schaaps detail.pdf", h=6, w=5)




# view areas from the vector mode
as.data.frame(range.indiv)


# compare areas to those calculated with the fixed bandwidth
range.fixed <- getverticeshr(kde_med)
as.data.frame(range.fixed)

################################################################################

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


# For barn coordinates, note that Hepp and Mary Ann are not primary tagging sites. One female from
# Cooks moved to nest at Hepp for her replacement clutch. 
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







