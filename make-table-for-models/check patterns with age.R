
################################################################################
# Check relationships with female age
# Heather Kenny-Duddela
# Nov 17, 2024
################################################################################

# set working directory
setwd("~/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/make-table-for-models")

# load data
birds <- read.csv("output-files/table for movement and mating models_BES_R.csv")
age <- read.csv("input-files/table of estimated ages 2021.csv")
mass <- read.csv("input-files/plumage-color-mass.csv")

# libraries
library(dplyr)

## add ages to bird table

# add bird id from mass table
colnames(mass)[1:2] <- c("bird", "band")

age2 <- left_join( mass[ ,1:2], age, by="band")

# add ages to bird table
birds2 <- left_join(birds, age2, by="bird")

# add age category
birds2$age.cat <- ifelse(birds2$est.age>1, "ASY", "SY")

# remove Cooks-27
birds3 <- subset(birds2, birds2$FamilyID_mom!="Cooks-27")

#-------------------------------------------------------------------------------

# age category with female plumage
ggplot(birds3, aes(x=age.cat, y=B_avg.bright.female)) +
  geom_boxplot() +
  geom_point(aes(color=certainty), 
             position=position_jitter(h=0, w=0.1))+
  ggtitle("Female color by age category") +
  xlab("Age category") + ylab("Female belly brightness")

# Figure S3
ggsave("output-files/female color by age category.png", h=3, w=4)

# Wilcoxon test to compare groups
wilcox.test(B_avg.bright.female ~ age.cat, data=birds3, alternative="two.sided")
# Wilcoxon rank sum exact test
# 
# data:  B_avg.bright.female by age.cat
# W = 12, p-value = 0.8333
# alternative hypothesis: true location shift is not equal to 0


#-------------------------------------------------------------------------------
# Check correlation of female mass with group size
#-------------------------------------------------------------------------------


ggplot(birds, aes(x=mass, y=group.size)) +
  geom_point()

cor.test(birds$mass, birds$group.size, method="spearman")
# data:  birds$mass and birds$group.size
# S = 165.27, p-value = 0.4607
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.2487608
