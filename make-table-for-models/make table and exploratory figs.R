
# GPS tagging 2021 - analysis for hypothesis testing
# Version for BES submission (June 2024)
# Heather Kenny-Duddela

# libraries
library(tidyverse)
library(ggplot2)

### Set working directory
setwd("C:/Users/heath/Documents/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/make-table-for-models")

### load data

# paternity data for each clutch
clutch <- read.csv("input-files/fert_2021_by_clutchID.csv")

# update fertilization type for the Hepp clutch, where the dad is from Cooks
# This should actually be classified as ep_same site since Cooks was the 
# original tagging site for this female
clutch$fert_type[which(clutch$clutch_id_ind1=="Hepp (near Cook)_1_1" &
                         clutch$Site_dad=="Cooks")] <- "ep_same"

# change factor order to match previous plots
clutch$fert_type <- factor(clutch$fert_type, 
                                           levels=c("wp","ep_same","mom_unk",
                                                    "dad_unk","ep_diff"))

# reorder brood factors
clutch$Brood_ind1 <- factor(clutch$Brood_ind1, levels=c("collected","1","2"))


# plot by fertilization type and facet by brood number
ggplot(clutch, aes(x=FamilyID_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Family ID") + ylab("Number of\n fertilizations") + 
  facet_grid(Brood_ind1~.) +
  ggtitle("Fertilization types across broods 2021")


# calculate proportion within-pair for all clutches
clutch.wp <- clutch %>%
  group_by(clutch_id_ind1, FamilyID_mom, Brood_ind1, Site_ind1) %>%
  summarise(num.wp = sum(subset(fert, fert_type=="wp")),
            num.ep = sum(subset(fert, fert_type!="wp")),
            clutch.size = mean(clutch_size),
            prop.wp = num.wp/clutch.size,
            prop.ep = num.ep/clutch.size)


### Subset tagged females only

# KDE estimates which are in hectares
areas <- read.csv("input-files/KDE probs_updated.csv")


# add column for family IDs to areas data
areas$FamilyID_mom <- NA

areas$FamilyID_mom[which(areas$bird=="bc09")] <- "BlueCloud-09"
areas$FamilyID_mom[which(areas$bird=="bc14")] <- "BlueCloud-14"
areas$FamilyID_mom[which(areas$bird=="ca03")] <- "Cathys-03"
areas$FamilyID_mom[which(areas$bird=="co23")] <- "Cooks-23"
areas$FamilyID_mom[which(areas$bird=="co27")] <- "Cooks-27"
areas$FamilyID_mom[which(areas$bird=="co31")] <- "Cooks-31"
areas$FamilyID_mom[which(areas$bird=="mb26")] <- "MakeBelieve-26"
areas$FamilyID_mom[which(areas$bird=="mb69")] <- "MakeBelieve-69"
areas$FamilyID_mom[which(areas$bird=="sc108")] <- "Schaaps-108"
areas$FamilyID_mom[which(areas$bird=="sc80")] <- "Schaaps-80"
areas$FamilyID_mom[which(areas$bird=="st20")] <- "Struthers-20"



# subset clutch data
clutch.move <- subset(clutch, clutch$FamilyID_mom %in% areas$FamilyID_mom)

# subset wp data
clutch.wp.move <- subset(clutch.wp, clutch.wp$FamilyID_mom %in% areas$FamilyID_mom)


# plot paternity proportions with just GPS females
# Manuscript Figure S1

ggplot(clutch.move, aes(x=FamilyID_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Female ID") + ylab("Number of\n fertilizations") + 
  facet_grid(Brood_ind1~.) +
  ggtitle("Fertilization types across broods for 2021 tagged females")

ggsave("output-files/FigS1_fert by brood 2021 tagged females.png", w=5.5, h=4)
ggsave("output-files/FigS1_fert by brood 2021 tagged females.pdf", w=5.5, h=4)



# calculate proportion wp in collected and replacement broods with SD

# just collected brood
clutch.wp.move.col <- subset(clutch.wp.move, clutch.wp.move$Brood_ind1=="collected")
mean(clutch.wp.move.col$prop.wp) # average wp in collected 0.636 +- 0.466
sd(clutch.wp.move.col$prop.wp)

# just replacement brood
clutch.wp.move.rep <- subset(clutch.wp.move, clutch.wp.move$Brood_ind1=="1")
mean(clutch.wp.move.rep$prop.wp) # average wp in collected 0.516 +- 0.392
sd(clutch.wp.move.rep$prop.wp)

# also calculate proportion ep in collected and replacement
# collected
mean(clutch.wp.move.col$prop.ep) # average ep in collected 0.363 +- 0.466
sd(clutch.wp.move.col$prop.ep)
# replacement
mean(clutch.wp.move.rep$prop.ep) # average ep in replacement 0.483 +- 0.392
sd(clutch.wp.move.rep$prop.ep)

###############################################################################

# Make the dataframe for the models

# take only replacement clutches (brood 1) from the tagged females
dat <- subset(clutch.wp.move, clutch.wp.move$Brood_ind1==1)

# add back in missing Cooks-27 whose replacement clutch failed
dat <- rbind(dat, clutch.wp.move[7,])
dat[11,c(1,3,5:9)] <- NA

# collected clutches only
dat1.2 <- subset(clutch.wp.move, clutch.wp.move$Brood_ind1=="collected")
# change labels to indicate collected
colnames(dat1.2)[5:9] <- c("num.wp.collected", "num.ep.collected", 
                           "clutch.size.collected", "prop.wp.collected", 
                           "prop.ep.collected")

# add collected clutches to main data
dat1.3 <- left_join(dat, dat1.2[,c(2,5:9)], by="FamilyID_mom")

# add areas data
dat2 <- left_join(dat1.3, areas, by="FamilyID_mom")

# change column names 
colnames(dat2)[29:30] <- c("proportion.close.50","proportion.close.100")

# correct Cooks site for the Hepp bird
dat2[dat2$FamilyID_mom=="Cooks-31",]$Site_ind1 <- "Cooks"

# load data with number of mates
mates.replace <- read.csv("input-files/num mates replacement clutch.csv")
mates.collect <- read.csv("input-files/num mates collected clutch.csv")

colnames(mates.replace)[c(5,7)] <- c("mates.replacement", "chicks.replacement")
colnames(mates.collect)[c(5,7)] <- c("mates.collected", "chicks.collected")

# use belly avg brightness
color <- read.csv("input-files/plumage-color-mass.csv")
# pull out only nestID and female and male color, and female mass
color2 <- color[,c(1,5,7,8)]
colnames(color2)[1] <- "bird"

# add estimated max distance from barn
est.max <- read.csv("input-files/estimated max distance from resampling_original.csv")
colnames(est.max)[c(2,6)] <- c("bird","est.max.dist")

dat2.2 <- left_join(dat2, color2, by="bird")

dat2.4 <- left_join(dat2.2, est.max[,c(2,6)], by="bird")

dat2.5 <- left_join(dat2.4, mates.collect[,c(2,5,7)], by="FamilyID_mom")

dat3 <- left_join(dat2.5, mates.replace[,c(2,5,7)], by="FamilyID_mom")

# add columns for EP mates, rather than total including WP male

# for collected, Schaaps-80 and BlueCloud-09 did not mate with social male
# subtract 1 from all rows except these two
dat3$mates.collected.ep <- NA
dat3$mates.collected.ep[-c(1,9)] <- dat3$mates.collected[-c(1,9)] -1
dat3$mates.collected.ep[c(1,9)] <- dat3$mates.collected[c(1,9)]
# for replacement, MakeBelieve-69 and Cooks-31 did not mate with their social male
# subtract 1 from all rows except these two
dat3$mates.replacement.ep <- NA
dat3$mates.replacement.ep[-c(4,5)] <- dat3$mates.replacement[-c(4,5)] -1
dat3$mates.replacement.ep[c(4,5)] <- dat3$mates.replacement[c(4,5)]


colnames(dat3)[c(31,32)] <- c("B_avg.bright.female", "B_avg.bright.male")

# try log-transforming area metrics to make them less skewed

dat3$log.kde90 <- log(dat3$prob_90)
hist(dat3$prob_90)
hist(dat3$log.kde90)

dat3$log.kde50 <- log(dat3$prob_50)
hist(dat3$prob_50)
hist(dat3$log.kde50)

dat3$log.max.dist <- log(dat3$est.max.dist)
hist(dat3$est.max.dist)
hist(dat3$log.max.dist)

hist(dat3$num.ep)

dat4 <- dat3

# add breeding group sizes
dat4$group.size <- c(5,5,5,5,5,5,2,10,10,1,5)

# calculate change in proportion EP and number EP young between collected and 
# replacement
dat4$diff.ep.prop <- dat4$prop.ep - dat4$prop.ep.collected

# change in number EP chicks
dat4$diff.num.ep <- dat4$num.ep - dat4$num.ep.collected

# change in number of EP sires
dat4$diff.ep.sires <- dat4$mates.replacement.ep - dat4$mates.collected.ep

# save final table
write.csv(dat4, "output-files/table for movement and mating models_BES.csv", row.names = F) 

################################################################################
# Correlation plot among all variables------------------------------------------

library(psych)

# correlations among predictor variables: log kde50, log kde90, prop.close.50,
# prop.close.100, est max dist, female plumage, male plumage, group size

# note that the MASS package masks the select function from dplyr, so you
# can't run the below line if MASS is currently loaded

vis <- select(dat4,log.kde50, log.kde90, log.max.dist, proportion.close.50,
              proportion.close.100, B_avg.bright.female, B_avg.bright.male, 
              mass, group.size)

pairs.panels(vis[,-c(1:3)],
             method="spearman",
             density=T,
             ellipses=F,
             smooth=F,
             stars=T)

# Manuscript Fig S3
# correlations among just the movement variables
vis2 <-  select(dat4,log.kde50, log.kde90, log.max.dist, proportion.close.50,
                proportion.close.100)
png("output-files/FigS3_correlation plot movement vars.png", pointsize=10, width=4500, height=3000, res=600)
pairs.panels(vis2[,4:8],
             method="spearman",
             density=T,
             ellipses=F,
             smooth=F,
             stars=T)
dev.off()



# correlation between collected and replacement prop ep
ggplot(dat4, aes(x=prop.ep, y=prop.ep.collected)) + geom_point(size=3) + 
  xlab("Proportion EP young in replacement") +
  ylab("Proportion EP young in collected clutch") +
  ggtitle("Mating patterns of females across broods")

cor.test(dat4$prop.ep, dat4$prop.ep.collected, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$prop.ep and dat4$prop.ep.collected
# S = 203.6, p-value = 0.5154
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2339114 

# correlation between collected and replacement for number EP offspring
cor.test(dat4$num.ep, dat4$num.ep.collected, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$num.ep and dat4$num.ep.collected
# S = 215.6, p-value = 0.3888
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.3066423 


# correlation between collected and replacement for number EP mates
cor.test(dat4$mates.collected.ep, dat4$mates.replacement.ep, method="spearman")

# Spearman's rank correlation rho
# 
# data:  dat4$mates.collected.ep and dat4$mates.replacement.ep
# S = 196.83, p-value = 0.5934
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1928792 

ggplot(dat4, aes(x=mates.collected.ep, y=mates.replacement.ep)) + 
  geom_point(size=3, position=position_jitter(height=0.08, width=0.04),
             alpha=0.3) + 
  xlab("Number mates in collected clutch") +
  ylab("Number mates in replacement clutch") +
  ggtitle("Mating patterns of females across broods")


# histogram of diff between collected and replacement
ggplot(dat4, aes(x=diff.ep.prop)) + 
  geom_histogram(binwidth=0.3, fill="lightblue", color="black") +
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  ggtitle("Change in proportion of EP offspring") +
  xlab("Difference in proportion of EP offspring from \ncollected to replacement")

ggsave("diff in propEP histogram.png")

# histogram of diff in number EP offspring
ggplot(dat4, aes(x=diff.num.ep)) + 
  geom_histogram(fill="lightblue", color="black", binwidth=1) +
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  ggtitle("Change in number of EP offspring") +
  xlab("Difference in number of EP offspring from \ncollected to replacement")

ggsave("diff in number EP young histogram.png")

# histogram of diff in number EP sires
ggplot(dat4, aes(x=diff.ep.sires)) + 
  geom_histogram(fill="lightblue", color="black", binwidth=1) +
  geom_vline(xintercept=0, color="red", linetype="dashed")+
  ggtitle("Change in number of EP sires") +
  xlab("Difference in number of EP sires from \ncollected to replacement")

ggsave("diff in number EP sires histogram.png")

#-------------------------------------------------------------------------------
### Reaction norm plots for change in clutch size, EP mating--------------------

# make table to show change in paternity
pat.change1 <- dat4 %>% ungroup() %>%
  select(FamilyID_mom, num.ep, clutch.size, prop.ep, 
                     mates.replacement.ep)

pat.change2 <- dat4 %>% ungroup() %>%
  select(FamilyID_mom, num.ep.collected, clutch.size.collected,
                      prop.ep.collected, mates.collected.ep)

# add column for collected or replacement
pat.change1$clutch <- "replacement"
pat.change2$clutch <- "collected"

# make colnames the same
colnames(pat.change1)[5] <- "mates.ep"
colnames(pat.change2)[2:5] <- colnames(pat.change1)[2:5]

# combine
pat.change3 <- rbind(pat.change2, pat.change1)

# remove NA rows
pat.change3 <- subset(pat.change3, pat.change3$FamilyID_mom!="Cooks-27")

### Make plots for Figure S2

# clutch size (Fig S2 a)

# calculate line weights
clutch.lines <- data.frame(t1 = pat.change3$clutch.size[1:10], 
                    t2 = pat.change3$clutch.size[11:20])

clutch.lines2 <- data.frame(table(clutch.lines$t1, clutch.lines$t2)) %>%
  mutate(transition = paste(Var1, Var2))%>%
  pivot_longer(cols = c(Var1, Var2), names_to = "x", values_to = "y")%>%
  filter(Freq != 0 & !grepl("NA", transition ))

clutch.lines2$x[which(clutch.lines2$x=="Var1")] <- "collected"
clutch.lines2$x[which(clutch.lines2$x=="Var2")] <- "replacement"
clutch.lines2$y <- as.numeric(as.character(clutch.lines2$y))
colnames(clutch.lines2)[3:4] <- c("clutch", "clutch.size")

# make plot for change in clutch size
ggplot(clutch.lines2, aes(x = clutch, y = clutch.size)) +  
  geom_line(aes(linewidth=(Freq), group=transition), alpha=0.5) +
  scale_linewidth(range=c(1,3), breaks=c(1,2,3)) +
  geom_count(inherit.aes=F, data=pat.change3, 
             aes(x=clutch, y=clutch.size)) +
  xlab("Clutch type") +
  ylab("Number of offspring sampled")

ggsave("output-files/FigS2a_change in clutch size.png", w=3, h=2.5, scale=1.2)
  

# proportion EP offspring (Fig S2 b)

# table for prop EP line weights
propEP.lines <- data.frame(t1 = pat.change3$prop.ep[1:10], 
                           t2 = pat.change3$prop.ep[11:20])

propEP.lines2 <- data.frame(table(propEP.lines$t1, propEP.lines$t2)) %>%
  mutate(transition = paste(Var1, Var2))%>%
  pivot_longer(cols = c(Var1, Var2), names_to = "x", values_to = "y")%>%
  filter(Freq != 0 & !grepl("NA", transition ))

propEP.lines2$x[which(propEP.lines2$x=="Var1")] <- "collected"
propEP.lines2$x[which(propEP.lines2$x=="Var2")] <- "replacement"
propEP.lines2$y <- as.numeric(as.character(propEP.lines2$y))
colnames(propEP.lines2)[3:4] <- c("clutch", "prop.ep")

ggplot(propEP.lines2, aes(x = clutch, y = prop.ep)) +  
  geom_line(aes(linewidth=(Freq), group=transition), alpha=0.5) +
  scale_linewidth(range=c(1,1), breaks=c(1)) +
  geom_count(inherit.aes=F, data=pat.change3, 
             aes(x=clutch, y=prop.ep)) +
  xlab("Clutch type") +
  ylab("Proportion of EP offspring")

ggsave("output-files/FigS2b_change in prop EP.png", w=3, h=2.5, scale=1.2)


# Number EP  young (Fig S2 c)

# table for num EP young line weight
numEP.lines <- data.frame(t1 = pat.change3$num.ep[1:10], 
                           t2 = pat.change3$num.ep[11:20])

numEP.lines2 <- data.frame(table(numEP.lines$t1, numEP.lines$t2)) %>%
  mutate(transition = paste(Var1, Var2))%>%
  pivot_longer(cols = c(Var1, Var2), names_to = "x", values_to = "y")%>%
  filter(Freq != 0 & !grepl("NA", transition ))

numEP.lines2$x[which(numEP.lines2$x=="Var1")] <- "collected"
numEP.lines2$x[which(numEP.lines2$x=="Var2")] <- "replacement"
numEP.lines2$y <- as.numeric(as.character(numEP.lines2$y))
colnames(numEP.lines2)[3:4] <- c("clutch", "num.ep")

ggplot(numEP.lines2, aes(x = clutch, y = num.ep)) +  
  geom_line(aes(linewidth=(Freq), group=transition), alpha=0.5) +
  scale_linewidth(range=c(1,2), breaks=c(1,2)) +
  geom_count(inherit.aes=F, data=pat.change3, 
             aes(x=clutch, y=num.ep)) +
  xlab("Clutch type") +
  ylab("Number of EP offspring")

ggsave("output-files/FigS2c_change in num EP.png", w=3, h=2.5, scale=1.2)

# number of EP sires (Fig S3 d)

# table for calculating num EP sires line weights
sires.lines <- data.frame(t1 = pat.change3$mates.ep[1:10], 
                          t2 = pat.change3$mates.ep[11:20])

sires.lines2 <- data.frame(table(sires.lines$t1, sires.lines$t2)) %>%
  mutate(transition = paste(Var1, Var2))%>%
  pivot_longer(cols = c(Var1, Var2), names_to = "x", values_to = "y")%>%
  filter(Freq != 0 & !grepl("NA", transition ))

sires.lines2$x[which(sires.lines2$x=="Var1")] <- "collected"
sires.lines2$x[which(sires.lines2$x=="Var2")] <- "replacement"
sires.lines2$y <- as.numeric(as.character(sires.lines2$y))
colnames(sires.lines2)[3:4] <- c("clutch", "mates.ep")

ggplot(sires.lines2, aes(x = clutch, y = mates.ep)) +  
  geom_line(aes(linewidth=(Freq), group=transition), alpha=0.5) +
  scale_linewidth(range=c(1,4), breaks=c(1,2,3,4)) +
  geom_count(inherit.aes=F, data=pat.change3, 
             aes(x=clutch, y=mates.ep)) +
  xlab("Clutch type") +
  ylab("Number of EP sires")

ggsave("output-files/FigS2d_change in sires.png", w=3, h=2.5, scale=1.2)



# number of EP sires
ggplot(pat.change3, aes(x=clutch, y=mates.ep)) + 
  geom_count() + geom_line(aes(group=FamilyID_mom)) +
  xlab("Clutch type") +
  ylab("Number of EP sires") +
  ggtitle("Change in number of EP sires \nbetween time points")

ggsave("generated-files/change in num sires.png", h=3.5, w=4)

# calculate summary stats for mating -------------------------------------------

# collected
mean(dat4$prop.ep.collected) #0.363
sd(dat4$prop.ep.collected) # 0.465
range(dat4$prop.ep.collected, na.rm=T) # 0,1

mean(dat4$num.ep.collected) # 1.636364
sd(dat4$num.ep.collected) # 2.15744

mean(dat4$mates.collected.ep, na.rm=T) #0.45
range(dat4$mates.collected.ep, na.rm=T) #0, 1

# replacement
mean(dat4$prop.ep, na.rm=T) #0.4833
range(dat4$prop.ep, na.rm=T) # 0,1
sd(dat4$prop.ep, na.rm=T) # 0.392

mean(dat4$num.ep, na.rm=T) # 1.8
sd(dat4$num.ep, na.rm=T) # 1.47573

mean(dat4$mates.replacement.ep, na.rm=T) #0.9
range(dat4$mates.replacement.ep, na.rm=T) #0,2

# difference
mean(dat4$diff.ep.prop, na.rm=T) # 0.1833333
sd(dat4$diff.ep.prop, na.rm=T) # 0.692018

mean(dat4$diff.num.ep, na.rm=T) # 0.4
sd(dat4$diff.num.ep, na.rm=T) # 2.91357

mean(dat4$diff.ep.sires, na.rm=T) # 0.5
sd(dat4$diff.ep.sires, na.rm=T) # 0.9718253


# Histograms of variables-------------------------------------------------------

ggplot(dat3, aes(x=prop.ep)) + geom_histogram(color="black", fill="lightblue", binwidth = 0.1) +
  ggtitle("Proportion of EP young in replacement clutch") +
  ylab("Number of females") + xlab("Proportion of extra-pair young")

# ggsave("hist prop ep young.png", w=5, h=4)

ggplot(dat3, aes(x=num.ep)) + geom_histogram(color="black", fill="lightblue", binwidth=1) +
  ggtitle("Number of EP young in replacement clutch") +
  ylab("Number of females") + xlab("Number of extra-pair young")

ggplot(dat3, aes(x=num.ep.mates)) + geom_histogram(color="black", fill="lightblue", binwidth = 1) +
  ggtitle("Total number of EP mates during the season") +
  ylab("Number of females") + xlab("Number of extra-pair mates")

# ggsave("hist num ep mates.png", w=5, h=4)

# Scatterplots of relationships-------------------------------------------------

# Proportion of EP young in replacement clutch

ggplot(dat4, aes(mates.replacement.ep, prop.ep)) + geom_point() +
  xlab("Number of EP mates") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("Number EP mates vs. proportion of EP young")


ggplot(dat4, aes(log.kde90, prop.ep)) + geom_point() +
  xlab("Area of 90% KDE in km2") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("90% KDE area vs. \nproportion of EP young")

ggplot(dat4, aes(log.kde50, prop.ep)) + geom_point() +
  xlab("Area of 50% KDE in km2") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("50% KDE area vs. \nproportion of EP young")


ggplot(dat4, aes(proportion.close.50, prop.ep)) + geom_point() +
  xlab("Proportion of GPS points within 50m of barn") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("Points within 50m vs. proportion of EP young")

ggplot(dat4, aes(proportion.close.50, prop.ep)) + geom_point() +
  xlab("Proportion of GPS points within 50m of barn") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("Points within 50m vs. proportion of EP young") +
  geom_smooth(method=lm)

ggplot(dat4, aes(proportion.close.100, prop.ep)) + geom_point() +
  xlab("Proportion of GPS points within 100m of barn") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("Points within 100m vs. proportion of EP young")

ggplot(dat4, aes(group.size, prop.ep)) + geom_point(alpha=0.5) +
  xlab("Number of breeding pairs") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("Group size vs. proportion of EP young")+
  geom_smooth(method=lm, se=F)



# prop EP in replacement and group size

cor.test(dat4$prop.ep, dat4$group.size, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$prop.ep and dat4$group.size
# S = 177.67, p-value = 0.833
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.07680673 

# prop EP in collected vs. group size
cor.test(dat4$prop.ep.collected, dat4$group.size, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$prop.ep.collected and dat4$group.size
# S = 145.07, p-value = 0.3054
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3405829 



# EP mates in collected and group size
cor.test(dat4$group.size, dat4$mates.collected.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$mates.collected.ep
# S = 153.6, p-value = 0.367
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3018349 

# group size vs. num EP mates in replacement
ggplot(dat4, aes(group.size, mates.replacement.ep)) + 
  geom_point(alpha=0.3) +
  xlab("Number of breeding pairs") + 
  ylab("Number of EP mates in replacement") +
  ggtitle("Group size vs. Number of EP mates in replacement")+
  geom_smooth(method=lm, se=F)


cor.test(dat4$group.size, dat4$mates.replacement.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$mates.replacement.ep
# S = 169.91, p-value = 0.9349
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.02976467 

# Num EPY in replacement with possible confounders------------------------------

# group size
ggplot(dat4, aes(x=group.size, y=num.ep)) + 
  geom_point(position=position_jitter(w=0.1, h=0.05))+
  geom_smooth(method=lm, se=F)

cor.test(dat4$group.size, dat4$num.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$num.ep
# S = 183.03, p-value = 0.7638
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1092634 

# male color
ggplot(dat4, aes(x=B_avg.bright.male, y=num.ep)) + 
  geom_point() +
  geom_smooth(method=lm, se=F)

cor.test(dat4$B_avg.bright.male, dat4$num.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$num.ep
# S = 139.29, p-value = 0.6673
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1558245 

# Num EPY in collected with possible confounders -------------------------------

# group size
ggplot(dat4, aes(x=group.size, y=num.ep.collected)) + 
  geom_point(position=position_jitter(w=0.1, h=0.05))+
  geom_smooth(method=lm, se=F)

cor.test(dat4$group.size, dat4$num.ep.collected, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$num.ep.collected
# S = 139.97, p-value = 0.2714
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3637792

# male color
ggplot(dat4, aes(x=B_avg.bright.male, y=num.ep.collected)) + 
  geom_point() +
  geom_smooth(method=lm, se=F)

cor.test(dat4$B_avg.bright.male, dat4$num.ep.collected, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$num.ep.collected
# S = 111.75, p-value = 0.1242
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4920565

# group size vs. movement variables --------------------------------------------

ggplot(dat4, aes(x=group.size, y=log.kde90)) +
  geom_point() +
  ylab("Log transformed 90% KDE area") +
  xlab("Number of breeding pairs") +
  ggtitle("90% KDE vs. Group size") +
  geom_smooth(method=lm)

cor.test(dat4$group.size, dat4$log.kde90, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$log.kde90
# S = 156.1, p-value = 0.3863
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2904407 

ggplot(dat4, aes(x=group.size, y=log.kde50)) +
  geom_point() +
  ylab("Log transformed 50% KDE area") +
  xlab("Number of breeding pairs") +
  ggtitle("50% KDE vs. Group size")+
  geom_smooth(method=lm, se=F)



cor.test(dat4$group.size, dat4$log.kde50, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$log.kde50
# S = 144.49, p-value = 0.3014
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3432482 

ggplot(dat4, aes(x=group.size, y=log.max.dist)) +
  geom_point() +
  ylab("Estimated max distance (m)") +
  xlab("Number of breeding pairs") +
  ggtitle("Max distance vs. Group size")+
  geom_smooth(method=lm)

cor.test(dat4$group.size, dat4$est.max.dist, method="spearman")

ggplot(dat4, aes(x=group.size, y=proportion.close.50)) +
  geom_point() +
  ylab("Proportion GPS points within 50m of barn") +
  xlab("Number of breeding pairs") +
  ggtitle("Prop close 50m vs. Group size") +
  geom_smooth(method=lm, se=F)



cor.test(dat4$group.size, dat4$proportion.close.50, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$proportion.close.50
# S = 225.82, p-value = 0.9384
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.02646392 

ggplot(dat4, aes(x=group.size, y=proportion.close.100)) +
  geom_point() +
  ylab("Proportion GPS points within 100m of barn") +
  xlab("Number of breeding pairs") +
  ggtitle("Prop close 100m vs. Group size") +
  geom_smooth(method=lm, se=F)

cor.test(dat4$group.size, dat4$proportion.close.100, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$proportion.close.100
# S = 235.1, p-value = 0.841
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.06864963 

### Relationship between movement and plumage color----------------------------

# male plumage

ggplot(dat4, aes(x=B_avg.bright.male, y=log.kde90)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Log transformed 90% KDE area") +
  xlab("Social male belly average brightness") +
  ggtitle("90% KDE vs. Male plumage color")

cor.test(dat4$B_avg.bright.male, y=dat4$log.kde90, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$log.kde90
# S = 224, p-value = 0.9676
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.01818182 

ggplot(dat4, aes(x=B_avg.bright.male, y=log.kde50)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Log transformed 50% KDE area") +
  xlab("Social male belly average brightness") +
  ggtitle("50% KDE vs. Male plumage color")

cor.test(dat4$B_avg.bright.male, y=dat4$log.kde50, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$log.kde50
# S = 234, p-value = 0.8601
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.06363636 

ggplot(dat4, aes(x=B_avg.bright.male, y=log.max.dist)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Log transformed estimated max distance from barn") +
  xlab("Social male belly average brightness") +
  ggtitle("Max distance vs. Male plumage color")

cor.test(dat4$B_avg.bright.male, y=dat4$log.max.dist, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$log.max.dist
# S = 170, p-value = 0.5031
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2272727

ggplot(dat4, aes(x=B_avg.bright.male, y=proportion.close.50)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Proportion of GPS points within 50m of barn") +
  xlab("Social male belly average brightness") +
  ggtitle("Points within 50m vs. Male plumage color")

cor.test(dat4$B_avg.bright.male, y=dat4$proportion.close.50, 
         method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$proportion.close.50
# S = 205.97, p-value = 0.8522
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.06378149 

ggplot(dat4, aes(x=B_avg.bright.male, y=proportion.close.100)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Proportion of GPS points within 100m of barn") +
  xlab("Social male belly average brightness") +
  ggtitle("Points within 100m vs. Male plumage color")

cor.test(dat4$B_avg.bright.male, y=dat4$proportion.close.100, 
         method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$proportion.close.100
# S = 196, p-value = 0.7549
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1090909 


# female plumage

ggplot(dat4, aes(x=B_avg.bright.female, y=log.kde90)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Log transformed 90% KDE area") +
  xlab("Female belly average brightness") +
  ggtitle("90% KDE vs. Female plumage color")

cor.test(dat4$B_avg.bright.female, y=dat4$log.kde90, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.female and dat4$log.kde90
# S = 268, p-value = 0.5209
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2181818

ggplot(dat4, aes(x=B_avg.bright.female, y=log.kde50)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Log transformed 50% KDE area") +
  xlab("Female belly average brightness") +
  ggtitle("50% KDE vs. Female plumage color")

cor.test(dat4$B_avg.bright.female, y=dat4$log.kde50, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.female and dat4$log.kde50
# S = 246, p-value = 0.7343
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1181818 

ggplot(dat4, aes(x=B_avg.bright.female, y=log.max.dist)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Log transformed estimated max distance from barn") +
  xlab("Female belly average brightness") +
  ggtitle("Max distance vs. Female plumage color")



cor.test(dat4$B_avg.bright.female, y=dat4$log.max.dist, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.female and dat4$log.max.dist
# S = 292, p-value = 0.327
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.3272727 

ggplot(dat4, aes(x=B_avg.bright.female, y=proportion.close.50)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Proportion of GPS points within 50m of barn") +
  xlab("Female belly average brightness") +
  ggtitle("Points within 50m vs. Female plumage color")



cor.test(dat4$B_avg.bright.female, y=dat4$proportion.close.50, 
         method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.female and dat4$proportion.close.50
# S = 167.88, p-value = 0.4831
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2369027 

ggplot(dat4, aes(x=B_avg.bright.female, y=proportion.close.100)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Proportion of GPS points within 100m of barn") +
  xlab("Female belly average brightness") +
  ggtitle("Points within 100m vs. Female plumage color")

cor.test(dat4$B_avg.bright.female, y=dat4$proportion.close.100, 
         method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.female and dat4$proportion.close.100
# S = 242, p-value = 0.7757
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#  rho 
# -0.1 

### Relationship between mass and movement--------------------------------------

cor.test(dat4$mass, dat4$log.kde90, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$mass and dat4$log.kde90
# S = 144.83, p-value = 0.3037
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3416865 

cor.test(dat4$mass, dat4$log.kde50, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$mass and dat4$log.kde50
# S = 154.85, p-value = 0.3766
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2961283 

cor.test(dat4$mass, dat4$log.max.dist, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$mass and dat4$log.max.dist
# S = 185.92, p-value = 0.6493
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1548979 

cor.test(dat4$mass, dat4$proportion.close.50, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$mass and dat4$proportion.close.50
# S = 271.23, p-value = 0.4908
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2328767 

cor.test(dat4$mass, dat4$proportion.close.100, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$mass and dat4$proportion.close.100
# S = 228.02, p-value = 0.9153
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.03644656 

# mass and female plumage color
cor.test(dat4$mass, dat4$B_avg.bright.female, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$mass and dat4$B_avg.bright.female
# S = 354.31, p-value = 0.04606
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.6104799 

### Relationship between plumage and group size ###

ggplot(dat4, aes(x=group.size, y=B_avg.bright.female)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Female belly average brightness") +
  xlab("Number of breeding pairs at barn") +
  ggtitle("Female plumage vs. Group size")



cor.test(dat4$group.size, dat4$B_avg.bright.female, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$B_avg.bright.female
# S = 280.41, p-value = 0.4138
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2745985 



ggplot(dat4, aes(x=group.size, y=B_avg.bright.male)) +
  geom_point() +
  geom_smooth(method=lm, se=F)+
  ylab("Social male belly average brightness") +
  xlab("Number of breeding pairs at barn") +
  ggtitle("Social male plumage vs. Group size")



cor.test(dat4$group.size, dat4$B_avg.bright.male, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$B_avg.bright.male
# S = 204.9, p-value = 0.841
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.06864963 


# female and male plumage color
cor.test(dat4$B_avg.bright.female, dat4$B_avg.bright.male, method="spearman")


### Relationship between male plumage and EP mating----------------------------

cor.test(dat4$B_avg.bright.male, dat4$prop.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$prop.ep
# S = 128.33, p-value = 0.5371
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.2222603

cor.test(dat4$B_avg.bright.male, dat4$prop.ep.collected, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$prop.ep.collected
# S = 103.14, p-value = 0.09267
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5311965 

cor.test(dat4$B_avg.bright.male, dat4$mates.replacement.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$mates.replacement.ep
# S = 110.72, p-value = 0.3533
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3289758 

cor.test(dat4$B_avg.bright.male, dat4$mates.collected.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$mates.collected.ep
# S = 92.983, p-value = 0.0629
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5773503 


### Difference in EP mating ----------------------------------------------------

### movement
ggplot(dat4, aes(x=proportion.close.50, y=diff.ep.prop)) + geom_point() +
  geom_smooth(method=lm) + 
  ggtitle("Relationship between change in prop EP and \nmovement (Spearman -0.674, p=0.032)")



cor.test(dat4$proportion.close.50, dat4$diff.ep.prop,
         method="spearman", na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$proportion.close.50 and dat4$diff.ep.prop
# S = 276.34, p-value = 0.03231
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.6747752 

ggplot(dat4, aes(x=proportion.close.50, y=diff.num.ep)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP \nyoung and movement (Spearman -0.477, p=0.163)")



cor.test(dat4$proportion.close.50, dat4$diff.num.ep,
         method="spearman", na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$proportion.close.50 and dat4$diff.num.ep
# S = 243.72, p-value = 0.1633
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.4770665 

ggplot(dat4, aes(x=proportion.close.50, y=diff.ep.sires)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP sires \nand movement (Spearman -0.417, p=0.229)")



cor.test(dat4$proportion.close.50, dat4$diff.ep.sires,
         method="spearman", na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$proportion.close.50 and dat4$diff.ep.sires
# S = 233.96, p-value = 0.2294
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# -0.417917 


### female plumage
ggplot(dat4, aes(x=B_avg.bright.female, y=diff.ep.prop)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in prop EP and \nfemale plumage (Spearman 0.163, p=0.656)")



cor.test(dat4$diff.ep.prop, dat4$B_avg.bright.female, method="spearman",  
         na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$diff.ep.prop and dat4$B_avg.bright.female
# S = 138, p-value = 0.6567
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1636364 

ggplot(dat4, aes(x=B_avg.bright.female, y=diff.num.ep)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP young \nand female plumage (Spearman 0.396, p=0.257)")



cor.test(dat4$diff.num.ep, dat4$B_avg.bright.female, 
         method="spearman", na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$diff.num.ep and dat4$B_avg.bright.female
# S = 99.602, p-value = 0.2568
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3963488 

ggplot(dat4, aes(x=B_avg.bright.female, y=diff.ep.sires)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP sires \nand female plumage (Spearman 0.384, p=0.273)")



cor.test(dat4$diff.ep.sires, dat4$B_avg.bright.female, 
         method="spearman", na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$diff.ep.sires and dat4$B_avg.bright.female
# S = 101.62, p-value = 0.2732
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.384098 



# male plumage
ggplot(dat4, aes(x=B_avg.bright.male, y=diff.ep.prop)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in prop EP young \nand male plumage (Spearman -0.091, p=0.811)")



cor.test(dat4$B_avg.bright.male, dat4$diff.ep.prop, method="spearman", 
         na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$diff.ep.prop
# S = 180, p-value = 0.8114
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.09090909 


ggplot(dat4, aes(x=B_avg.bright.male, y=diff.num.ep)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP young \nand male plumage (Spearman -0.177, p=0.625)")



cor.test(dat4$B_avg.bright.male, dat4$diff.num.ep, method="spearman", 
         na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$diff.num.ep
# S = 194.18, p-value = 0.625
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1768326

ggplot(dat4, aes(x=B_avg.bright.male, y=diff.ep.sires)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP sires \nand male plumage (Spearman -0.019, p=0.957)")



cor.test(dat4$B_avg.bright.male, dat4$diff.ep.sires, method="spearman", 
         na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$B_avg.bright.male and dat4$diff.ep.sires
# S = 168.22, p-value = 0.9573
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.01953041


# Group size
ggplot(dat4, aes(x=group.size, y=diff.ep.prop)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in prop EP young \nand group size (Spearman -0.233, p=0.517)")



cor.test(dat4$group.size, dat4$diff.ep.prop, method="spearman", 
         na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$diff.ep.prop
# S = 203.45, p-value = 0.517
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2330462 

ggplot(dat4, aes(x=group.size, y=diff.num.ep)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP young \nand group size (Spearman -0.269, p=0.452)")



cor.test(dat4$group.size, dat4$diff.num.ep, method="spearman", 
         na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$diff.num.ep
# S = 209.38, p-value = 0.4524
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.2689527

ggplot(dat4, aes(x=group.size, y=diff.ep.sires)) + geom_point() +
  geom_smooth(method=lm) +
  ggtitle("Relationship between change in number of EP sires \nand group size (Spearman -0.143, p=0.692)")



cor.test(dat4$group.size, dat4$diff.ep.sires, method="spearman", 
         na.rm=T)
# Spearman's rank correlation rho
# 
# data:  dat4$group.size and dat4$diff.ep.sires
# S = 188.69, p-value = 0.6923
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1435727




