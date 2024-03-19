
################################################################################
# Script for resampling GPS points to estimate max distance traveled
# Heather Kenny-Duddela
# Nov 7, 2023
################################################################################

# Note: because this script uses a random sampling process, running it and generating
# a new final "estimated max distance from resampling.csv" file may result in
# values that differ slightly from the original estimated max distance values used
# in publication. 

# libraries
library(tidyverse)
library(ggplot2)

# Load GPS data
points <- read.csv("generated-datatables/point.distances.hkd.csv")

# summarise number of points per bird and day
summary <- points %>%
  group_by(nest, FIX.date2) %>%
  summarise(num.points = n())


# now calculate total number of points for each bird across all days
summary2 <- points %>%
  group_by(nest) %>%
  summarise(tot.points = n())


# calculate number of tracking days
tracking.days <- points %>%
  group_by(nest) %>%
  summarise(track.days = length(unique(FIX.date2)))

# combine total points and number of tracking days
days.tot.points <- left_join(summary2, tracking.days, by="nest")

# write.csv(days.tot.points, "total points and number tracking days by birds_v2.csv")

# smallest number of points is 61, so down sample to 60

# rename data frame and columns to minimize changes to following code
points.clean <- points
colnames(points.clean)[18] <- "Site"

# set up storage, resample each bird 100 times
# columns are birds, rows are resampled points

storage <- matrix(NA, 100,11)
colnames(storage) <- summary2$nest


# Loop to resample points for each bird
for (i in 1:length(summary2$nest)) {
   bird <- summary2$nest[i]
   rows <- summary2$tot.points[i]
   bird.points <- subset(points.clean, points.clean$Site == summary2$nest[i])
  # loop to repeat sampling 100 times
  for (j in 1:100) {
    sample <- sample(rows, 60, replace=F)
    sample.points <- bird.points$distance[sample]
    max <- max(sample.points)
    
    storage[j,i] <- max
   
  }
  
}

# convert storage to dataframe
storage.df <- as.data.frame(storage)

# visualize histograms

ggplot(storage.df, aes(x=bc09)) + geom_histogram()

ggplot(storage.df, aes(x=bc14)) + geom_histogram()

ggplot(storage.df, aes(x=ca03)) + geom_histogram()

ggplot(storage.df, aes(x=co23)) + geom_histogram()


# add range, mean, and sd to summary2

summary2$min <- NA
summary2$max <- NA
summary2$mean <- NA
summary2$sd <- NA

for (i in 1:length(summary2$nest)) {
  summary2$min[i] <- min(storage.df[,i])
  summary2$max[i] <- max(storage.df[,i])
  summary2$mean[i] <- mean(storage.df[,i])
  summary2$sd[i] <- sd(storage.df[,i])
}

summary2$diff.max.mean <- summary2$max - summary2$mean

# write.csv(summary2, "generated-datatables/estimated max distance from resampling.csv")


                         