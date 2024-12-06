
################################################################################
# Calculate daily movement metrics for each female
# Table S5
# Heather Kenny-Duddela
# Nov 18, 2024
################################################################################

# set working directory
setwd("~/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/process-spatial-data")

# load data
dist <- read.csv("generated-datatables/point.distances.hkd.csv")

# load libraries
library(dplyr)


# Calculate daily prop50m and prop100m -----------------------------------------

# summarise by trackday
dist.daily <- dist %>%
  group_by(nest, trackday) %>%
  summarise(num.points = n(),
            prop.past.50m = 1 - (sum(close50)/num.points),
            prop.past.100m = 1 - (sum(close100)/num.points),
            max.dist = max(distance),
            mean.dist = mean(distance))

# summarise daily measures by bird
dist.var.bird <- dist.daily %>%
  group_by(nest) %>%
  summarise(track.days = n(),
            mean.points = mean(num.points),
            sd.points = sd(num.points),
            mean.past.50m = mean(prop.past.50m),
            sd.past.50m = sd(prop.past.50m),
            mean.past.100m = mean(prop.past.100m),
            sd.past.100m = sd(prop.past.100m),
            mean.max.dist = mean(max.dist),
            sd.max.dist = sd(max.dist),
            daily.mean.dist = mean(mean.dist),
            daily.mean.dist.sd = sd(mean.dist))

write.csv(dist.var.bird, "generated-datatables/daily movement variation.csv", 
          row.names = F)

# summarise across all birds
dist.all <- dist %>%
  group_by(nest) %>%
  summarise(num.points = n(),
            prop.past.50m = 1 - (sum(close50)/num.points),
            prop.past.100m = 1 - (sum(close100)/num.points),
            max.dist = max(distance),
            mean.dist = mean(distance))

dist.all2 <- dist.all %>%
  summarise(mean.num.points = mean(num.points),
            sd.num.points = sd(num.points),
            mean.past.50m = mean(prop.past.50m),
            sd.past.50m = sd(prop.past.50m),
            mean.past.100m = mean(prop.past.100m),
            sd.past.100m = sd(prop.past.100m),
            mean.dist.tot = mean(mean.dist),
            sd.dist = sd(mean.dist, na.rm=T))

