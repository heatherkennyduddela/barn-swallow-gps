
# GPS tagging 2021 - analysis for hypothesis testing
# Heather Kenny-Duddela

# libraries
library(tidyverse)
library(ggplot2)

### load data

# paternity data for each clutch
clutch <- read.csv("input-data/fert_2021_by_clutchID.csv")

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
areas <- read.csv("input-data/KDE probs_updated.csv")


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
# Manuscript Figure S3

ggplot(clutch.move, aes(x=FamilyID_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Female ID") + ylab("Number of\n fertilizations") + 
  facet_grid(Brood_ind1~.) +
  ggtitle("Fertilization types across broods for 2021 tagged females")

ggsave("generated-files/FigS3_fert by brood 2021 tagged females.png", w=5.5, h=4)
ggsave("generated-files/FigS3_fert by brood 2021 tagged females.pdf", w=5.5, h=4)


# calculate proportion wp in collected and replacement broods with SD

# just collected brood
clutch.wp.move.col <- subset(clutch.wp.move, clutch.wp.move$Brood_ind1=="collected")
mean(clutch.wp.move.col$prop.wp) # average wp in collected 0.636 +- 0.466
sd(clutch.wp.move.col$prop.wp)

# just replacement brood
clutch.wp.move.rep <- subset(clutch.wp.move, clutch.wp.move$Brood_ind1=="1")
mean(clutch.wp.move.rep$prop.wp) # average wp in collected 0.516 +- 0.392
sd(clutch.wp.move.rep$prop.wp)


###############################################################################

# Make the dataframe for the models

# take only replacement clutches (brood 1) from the tagged females
dat <- subset(clutch.wp.move, clutch.wp.move$Brood_ind1==1)

# add back in missing Cooks-27 whose replacement clutch failed
dat <- rbind(dat, clutch.wp.move[7,])
dat[11,c(1,3,5:9)] <- NA

# add areas data
dat2 <- left_join(dat, areas, by="FamilyID_mom")

# change column names 
colnames(dat2)[24:25] <- c("proportion.close.50","proportion.close.100")

# correct Cooks site for the Hepp bird
dat2[dat2$FamilyID_mom=="Cooks-31",]$Site_ind1 <- "Cooks"

# load data with number of overall mates
mates <- read.csv("input-data/num mates chicks clutches per female 2021.csv")
mates.replace <- read.csv("input-data/num mates replacement clutch.csv")

colnames(mates.replace)[c(5,7)] <- c("mates.replacement", "chicks.replacement")

# use belly avg brightness
color <- read.csv("input-data/plumage-color.csv")
# pull out only nestID and female and male color
color2 <- color[,c(1,5,7)]
colnames(color2)[1] <- "bird"

# add estimated max distance from barn
est.max <- read.csv("input-data/estimated max distance from resampling_original.csv")
colnames(est.max)[c(2,6)] <- c("bird","est.max.dist")

dat2.2 <- left_join(dat2, color2, by="bird")

dat2.3 <- left_join(dat2.2, mates, by="FamilyID_mom")

dat2.4 <- left_join(dat2.3, est.max[,c(2,6)], by="bird")

dat3 <- left_join(dat2.4, mates.replace[,c(2,5,7)], by="FamilyID_mom")

# add columns for EP mates, rather than total including WP male
dat3$num.ep.mates <- dat3$num.mates -1
dat3$mates.replacement.ep <- dat3$mates.replacement -1

colnames(dat3)[c(26,27)] <- c("B_avg.bright.female", "B_avg.bright.male")

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

# calculate proportion of EP young over whole season

# number of each fertilization type
clutch.female <- clutch %>%
  group_by(Band_mom, fert_type) %>%
  summarise(tot.fert = sum(fert))

# total number of fertilizations per female
clutch.female.total <- clutch %>%
  group_by(Band_mom) %>%
  summarise(tot.clutch = sum(fert))

# combine total number of each type with grant total
clutch.female2 <- left_join(clutch.female, clutch.female.total, by="Band_mom")

# pull out only wp fertilizations
clutch.female.wp <- subset(clutch.female2, clutch.female2$fert_type=="wp")
# calculate proportion ep across the season
# subtract number wp from total to get number ep, then divide by total
clutch.female.wp$prop.ep.season <- (clutch.female.wp$tot.clutch-clutch.female.wp$tot.fert)/clutch.female.wp$tot.clutch

# add to main data table
dat4 <- left_join(dat3, clutch.female.wp[,c(1,5)], by="Band_mom")

ggplot(dat4, aes(prop.ep, prop.ep.season)) + geom_point()

# save final table
# write.csv(dat4, "generated-files/table for movement and mating models_v5.csv", row.names = F) 


# Correlation plot among all variables------------------------------------------

library(psych)

# Manuscript Figure S1
pairs.panels(dat4[,c(39:40,24,25, 26, 27)],
             method="spearman",
             density=T,
             ellipses=F,
             smooth=F)


# correlation between replacement and season prop ep
ggplot(dat4, aes(x=prop.ep, y=prop.ep.season)) + geom_point(size=3) + 
  xlab("Proportion EP young in replacement") +
  ylab("Proportion EP young over season") +
  ggtitle("Mating patterns of females \nacross time scales")

# ggsave("female mating across timescales.png")


cor.test(dat4$prop.ep, dat4$prop.ep.season, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$prop.ep and dat4$prop.ep.season
# S = 114.07, p-value = 0.3855
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3086949 


# calculate summary stats for mating and movement------------------------------

mean(dat4$prop.ep, na.rm=T) #0.4833
range(dat4$prop.ep, na.rm=T) # 0,1
mean(dat4$num.ep.mates, na.rm=T) #1.636364
range(dat4$num.ep.mates, na.rm=T) #0, 3
mean(dat4$mates.replacement.ep, na.rm=T) #0.7
range(dat4$mates.replacement.ep, na.rm=T) #0, 2



# Histograms of variables-------------------------------------------------------

ggplot(dat3, aes(x=prop.ep)) + geom_histogram(color="black", fill="lightblue", binwidth = 0.1) +
  ggtitle("Proportion of EP young in replacement clutch") +
  ylab("Number of females") + xlab("Proportion of extra-pair young")

# ggsave("hist prop ep young.png", w=5, h=4)

ggplot(dat3, aes(x=num.ep.mates)) + geom_histogram(color="black", fill="lightblue", binwidth = 1) +
  ggtitle("Total number of EP mates during the season") +
  ylab("Number of females") + xlab("Number of extra-pair mates")

# ggsave("hist num ep mates.png", w=5, h=4)

# Scatterplots of relationships-------------------------------------------------

# Proportion of EP young in replacement clutch

ggplot(dat4, aes(num.ep.mates, prop.ep)) + geom_point() +
  xlab("Number of EP mates") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("Number EP mates vs. proportion of EP young")


ggplot(dat4, aes(prob_90/100, prop.ep)) + geom_point() +
  xlab("Area of 90% KDE in km2") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("90% KDE area vs. \nproportion of EP young")

# ggsave("movement plumage and paternity plots 2021/KDE 90% vs. prop ep young.png", w=4, h=3)

ggplot(dat4, aes(prob_50/100, prop.ep)) + geom_point() +
  xlab("Area of 50% KDE in km2") + 
  ylab("proportion of EP young in \nreplacement clutch") +
  ggtitle("50% KDE area vs. \nproportion of EP young")

# ggsave("movement plumage and paternity plots 2021/KDE 50% vs. prop ep young.png", w=4, h=3)

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


# Number of EP mates over the whole season

ggplot(dat4, aes(prob_90, num.ep.mates)) + geom_point() +
  xlab("Area of 90% KDE") + 
  ylab("Number of EP mates") +
  ggtitle("90% KDE area vs. Number of EP mates")

ggplot(dat4, aes(prob_50, num.ep.mates, color=Site_ind1)) + geom_point() +
  xlab("Area of 50% KDE") + 
  ylab("Number of EP mates") +
  ggtitle("50% KDE area vs. Number of EP mates")

ggplot(dat4, aes(prob_50, num.ep.mates)) + geom_point() +
  xlab("Area of 50% KDE") + 
  ylab("Number of EP mates") +
  ggtitle("50% KDE area vs. Number of EP mates") +
  geom_smooth(method=lm)

ggplot(dat4, aes(proportion.close.50, num.ep.mates)) + geom_point() +
  xlab("Proportion of GPS points within 50m of barn") + 
  ylab("Number of EP mates") +
  ggtitle("Points within 50m vs. Number of EP mates")

ggplot(dat4, aes(proportion.close.100, num.ep.mates)) + geom_point() +
  xlab("Proportion of GPS points within 100m of barn") + 
  ylab("Number of EP mates") +
  ggtitle("Points within 100m vs. Number of EP mates")

ggplot(dat4, aes(proportion.close.100, num.ep.mates)) + geom_point() +
  xlab("Proportion of GPS points within 100m of barn") + 
  ylab("Number of EP mates") +
  ggtitle("Points within 100m vs. Number of EP mates") +
  geom_smooth(method=lm)



################################################################################
--------------------------------------------------------------------------------




# Models!-----------------------------------------------------------------------

library(broom)# for calculating confidence intervals around coefficient estimates
library(AICcmodavg) # for calculating AICc for model comparison

# helpful threads about offsets and weights in Poisson and binomial models
# https://stats.stackexchange.com/questions/264071/how-is-a-poisson-rate-regression-equal-to-a-poisson-regression-with-correspondin

# https://stats.stackexchange.com/questions/297859/can-weights-and-offset-lead-to-similar-results-in-poisson-regression 

# Help about glms from UCLA stats page
# https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faq-how-do-i-interpret-odds-ratios-in-logistic-regression/

# https://stats.oarc.ucla.edu/r/dae/logit-regression/


# Proportion of EP in replacement clutch predicted by different movement components

# Model proportion EP as binomial with an offset for brood size
# for binomial family, response can be specified as numerical value 0 to 1 as proportion of 
# successful cases, with total number of cases given by weights



# 90% KDE
kde90 <- glm( prop.ep ~ scale(log.kde90) + B_avg.bright.female, 
              weights=clutch.size, data=dat4, family="binomial")

summary(kde90)
hist(resid(kde90, type="pearson"))
plot(fitted(kde90) ~ resid(kde90, type="pearson"))
confint(kde90)


# 50% KDE
kde50 <- glm( prop.ep ~ scale(log.kde50) + B_avg.bright.female, 
              weights=clutch.size, data=dat4, family="binomial")

summary(kde50)
confint(kde50)


# prop within 50m
prop50m <- glm(prop.ep ~ proportion.close.50 + B_avg.bright.female, 
               weights=clutch.size, data=dat4, 
               family="binomial")
summary(prop50m)
hist(resid(prop50m, type="pearson"))
plot(fitted(prop50m) ~ resid(prop50m, type="pearson"))
confint(prop50m)


# prop within 100m
prop100m <- glm(prop.ep ~ proportion.close.100 + B_avg.bright.female, weights=clutch.size, data=dat4,
                family="binomial")
summary(prop100m)
confint(prop100m)


# est max distance
estmax <- glm(prop.ep ~ scale(log.max.dist) + B_avg.bright.female, 
              weights=clutch.size, data=dat4,
                 family="binomial")
summary(estmax)
confint(estmax)



#### Model number of EP mates as Poisson --------------------------------------

# For the poisson model, compare models with and without the number of chicks
# covariate, to see which fits the data best



#-------------------------------------------------------------------------------
# Model number of EP mates in replacement clutch as Poisson

# first check for overdispersion
hist(dat4$mates.replacement.ep)
mean(dat4$mates.replacement.ep, na.rm=T) # 0.7
var(dat4$mates.replacement.ep, na.rm=T) # 0.455

null.rep.poisson <- glm(mates.replacement.ep ~ chicks.replacement , 
                        data=dat4, family="poisson")
AICc(null.rep.poisson) # 25.94659

# quasi-poisson can model
null.rep.quasi <- glm(mates.replacement.ep ~ chicks.replacement, data=dat4, family="quasipoisson")
summary(null.rep.quasi)  # (Dispersion parameter for quasipoisson family taken to be 0.6912519)

# negative binomial model
null.rep.nb <- glm.nb(mates.replacement.ep ~ chicks.replacement, data=dat4)
summary(null.rep.nb)
AICc(null.rep.nb) # 30.23244

### If anything the data is under dispersed, and the Poisson model performs better again. 
### proceed with fitting models using Poisson variance


# 90% KDE
mates.rep.kde90.cov <- glm(mates.replacement.ep ~ scale(log.kde90) + chicks.replacement +
                         B_avg.bright.female, 
                       data=dat4, family="poisson")
summary(mates.rep.kde90.cov)
confint(mates.rep.kde90.cov)
hist(resid(mates.rep.kde90.cov, type="pearson"))
plot(fitted(mates.rep.kde90.cov) ~ resid(mates.rep.kde90.cov, type="pearson"))

# 50% KDE
mates.rep.kde50.cov <- glm(mates.replacement.ep ~ scale(log.kde50) + chicks.replacement +
                             B_avg.bright.female, 
                           data=dat4, family="poisson")
summary(mates.rep.kde50.cov)
confint(mates.rep.kde50.cov)
hist(resid(mates.rep.kde50.cov, type="pearson"))
plot(fitted(mates.rep.kde50.cov) ~ resid(mates.rep.kde50.cov, type="pearson"))

# proportion points within 50
mates.rep.prop50m.cov <- glm(mates.replacement.ep ~ proportion.close.50 + chicks.replacement +
                             B_avg.bright.female, 
                           data=dat4, family="poisson")
summary(mates.rep.prop50m.cov)
confint(mates.rep.prop50m.cov)

# proportion points within 100m
mates.rep.prop100m.cov <- glm(mates.replacement.ep ~ proportion.close.100 + 
                                chicks.replacement +
                               B_avg.bright.female, 
                             data=dat4, family="poisson")
summary(mates.rep.prop100m.cov)
confint(mates.rep.prop100m.cov)

# est max distance
mates.rep.estmax.cov <- glm(mates.replacement.ep ~ scale(log.max.dist) + 
                              chicks.replacement +
                                B_avg.bright.female, 
                              data=dat4, family="poisson")
summary(mates.rep.estmax.cov)
confint(mates.rep.estmax.cov)

################################################################################

#-------------------------------------------------------------------------------
# AIC model selection - proportion EP in replacement clutch with female plumage color


propEP.null <- glm(prop.ep ~ scale(B_avg.bright.female), family = "binomial", data = dat4, weights = clutch.size)
summary(propEP.null)

# candidate set of models
propEP.models <- list(propEP.null, kde50, kde90, prop100m, prop50m, estmax)
# model names
propEP.modnames <- c("null","kde50","kde90","prop100m","prop50m","estmax")

propEP.table <- aictab(propEP.models, propEP.modnames)

# Manuscript Table S3
# top model is prop50m, but is the same as the null model!
write.csv(propEP.table, "generated-files/TableS3_AIC table for propEP_v5.csv", row.names=F)



#------------------------------------------------------------------------------
### AIC model selection for number ep mates in replacement clutch

mates.rep.null <- glm(mates.replacement.ep ~ 1, data=dat4, family="poisson")

mates.rep.cov <- glm(mates.replacement.ep ~ chicks.replacement, 
                     data=dat4, family="poisson")

mates.rep.plum.cov <- glm(mates.replacement.ep ~ chicks.replacement +
                            scale(B_avg.bright.female), 
                          data=dat4, family="poisson")
summary(mates.rep.plum.cov)



### tables for models with plumage color and number chicks

# candidate set of models
mates.rep.models.plum <- list(mates.rep.plum.cov, 
                                mates.rep.kde90.cov, mates.rep.kde50.cov, 
                                mates.rep.prop50m.cov, mates.rep.prop100m.cov,
                                mates.rep.estmax.cov)

mates.rep.modmanes.plum <- c("mates.rep.plum.cov", 
                               "mates.rep.kde90.cov", "mates.rep.kde50.cov", 
                               "mates.rep.prop50m.cov", "mates.rep.prop100m.cov",
                             "mates.rep.estmax.cov")

# Manuscript Table S4
mates.rep.plum.table <- aictab(mates.rep.models.plum, mates.rep.modmanes.plum)
write.csv(mates.rep.plum.table, "generated-files/TableS4_AIC table for mates in replacement_v5.csv", row.names=F)


#-------------------------------------------------------------------------------
# Compare the following models to test hypotheses for: 
# 1) space use
# 2) space use interacts with female plumage
# 3) space use interacts with male plumage


# LIKELIHOOD RATIO TESTS FOR HYPOTHESIS TESTING

# 1) movement only
# 2) movement + female plumage
# 3) movement * female plumage
# 4) movement + male plumage
# 5) movement * male plumage

# Explanation from Kayleigh about type II and type III tests
# Type III tests look at an individual term, keeping everything else in the model. In a model with an interaction term, this is generally inappropriate to use for testing main effects, since we would almost never want a model with an area-plumage interaction but no main effect term for area. That’s part of the general rule to always include lower-order terms in a model with interactions or polynomials.
# 
# Type II tests look at an individual term, keeping all other main effects (but not interactions) in the model. So it can be appropriate to look at the main effect results from a Type II model even if there is an interaction, although I would generally recommend just fitting the simpler model (without interaction) if the interaction term was not significant.
# 
# In my opinion, Type II is always sensible, whereas Type I rarely is and Type III only sometimes is. So I think car::Anova() gets it right by defaulting to Type II. The interaction term results are the same for both tests, since in that case the interaction term is tested keeping all the main effects in.
# 
# The general strategy I typically recommend for a model with an interaction term like yours is:
#   
#   Fit a model with A, B, and A:B.
# Test the A:B term (Type 2 and Type 3 will be identical)
# If A:B term is significant, this is your model. Calculate the summaries of interest, e.g. using emmeans()
# If A:B term is not significant, go to step 3.
# Fit a model with A and B (but no interaction)
# Test A and B in the model from (3); this can be with either Type 2 or 3 since there is no interaction they give same results. Calculate the summaries of interest.
# 
# The reason I typically go this route, as opposed to using a Type II test of A or B in the model from (1), is that it is much more difficult to interpret/summarize effects when there is an interaction. And if the interaction is not (statistically) significant, then you’re dealing with a lot of complexity that you don’t need, so its easier to go the simpler model (and you get a slight amount of statistical power back when you do that).


### Results from model summaries and confidence interval calculations below
### are reported in manuscript Table 2 and Tables S6, S7, and S8


#-----------------------------------------------------------------------------

# Models for prop EP young in replacement clutch
# Raw values of plumage used instead of scaled

# 1) Null
propEP.null <- glm(prop.ep ~ 1, weights=clutch.size, data=dat3, family="binomial")
summary(propEP.null)

# 2) just prop50m
propEP.prop50m <- glm(prop.ep ~ proportion.close.50, weights=clutch.size, 
                      data=dat3, family="binomial")
summary(propEP.prop50m)
confint(propEP.prop50m)


# 3) just female plumage
# Check if model will converge with unscaled plumage color - yes
propEP.plumF <- glm(prop.ep ~ B_avg.bright.female, weights=clutch.size, 
                    data=dat3, family="binomial")
summary(propEP.plumF)

# 4) just male plumage
propEP.plumM <- glm(prop.ep ~ B_avg.bright.male, weights=clutch.size, 
                    data=dat3, family="binomial")
summary(propEP.plumM)

# 5) prop50m + female plumage
propEP.prop50m.plumF <- glm(prop.ep ~ proportion.close.50 + B_avg.bright.female,
                            weights=clutch.size, 
                      data=dat3, family="binomial")
summary(propEP.prop50m.plumF)

# 6) prop50m * female plumage
propEP.prop50m.plumF.int <- glm(prop.ep ~ proportion.close.50 * B_avg.bright.female, 
                            weights=clutch.size, 
                            data=dat3, family="binomial")
summary(propEP.prop50m.plumF.int)
confint(propEP.prop50m.plumF.int)

# 7) prop50m + male plumage
propEP.prop50m.plumM <- glm(prop.ep ~ proportion.close.50 + B_avg.bright.male,
                            weights=clutch.size, 
                            data=dat3, family="binomial")
summary(propEP.prop50m.plumM)
confint(propEP.prop50m.plumM)

# 8) prop50m * male plumage
propEP.prop50m.plumM.int <- glm(prop.ep ~ proportion.close.50 * B_avg.bright.male, 
                                weights=clutch.size, 
                                data=dat3, family="binomial")
summary(propEP.prop50m.plumM.int)
confint(propEP.prop50m.plumM.int)

########### Likelihood ratio test for prop EP in replacement ------------------

library(car)

# 1) movement only
Anova(propEP.prop50m)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep
# LR Chisq Df Pr(>Chisq)   
# proportion.close.50   7.3378  1   0.006752 **
confint(propEP.prop50m)
# 2.5 %   97.5 %
#   (Intercept)          -0.1483768  1.94504
# proportion.close.50 -20.3613506 -2.63957
summary(propEP.prop50m)
# glm(formula = prop.ep ~ proportion.close.50, family = "binomial", 
#     data = dat3, weights = clutch.size)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3213  -0.9324   0.3345   1.1024   1.8888  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)           0.8487     0.5256   1.615   0.1064  
# proportion.close.50 -10.3358     4.3817  -2.359   0.0183 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 27.223  on 9  degrees of freedom
# Residual deviance: 19.885  on 8  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 33.503
# 
# Number of Fisher Scoring iterations: 5

# 2) movement + female plumage
Anova(propEP.prop50m.plumF)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep
# LR Chisq Df Pr(>Chisq)   
# proportion.close.50   7.3885  1   0.006564 **
#   B_avg.bright.female   0.0520  1   0.819588  
confint(propEP.prop50m.plumF)
# 2.5 %    97.5 %
#   (Intercept)          -4.7132233  5.391950
# proportion.close.50 -20.6726602 -2.686729
# B_avg.bright.female  -0.1688116  0.209322
summary(propEP.prop50m.plumF)
# glm(formula = prop.ep ~ proportion.close.50 + B_avg.bright.female, 
#     family = "binomial", data = dat3, weights = clutch.size)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.3933  -0.9169   0.2791   1.1364   1.8819  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)           0.28681    2.51509   0.114   0.9092  
# proportion.close.50 -10.45818    4.44307  -2.354   0.0186 *
#   B_avg.bright.female   0.02149    0.09412   0.228   0.8194  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 27.223  on 9  degrees of freedom
# Residual deviance: 19.833  on 7  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 35.451
# 
# Number of Fisher Scoring iterations: 5


# Here, the interaction term is significant at the 0.1 confidence level, so proceed to type III test next
# 3) movement * female plumage
Anova(propEP.prop50m.plumF.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: prop.ep
# LR Chisq Df Pr(>Chisq)  
# proportion.close.50                       6.2929  1    0.01212 *
#   B_avg.bright.female                       4.1567  1    0.04147 *
#   proportion.close.50:B_avg.bright.female   5.3427  1    0.02081 *
confint(propEP.prop50m.plumF.int)
# 2.5 %       97.5 %
#   (Intercept)                                1.3185340  35.82346847
# proportion.close.50                     -329.7985183 -28.45501651
# B_avg.bright.female                       -1.2951323  -0.01947026
# proportion.close.50:B_avg.bright.female    0.6946215  11.59544225
summary(propEP.prop50m.plumF.int)
# glm(formula = prop.ep ~ proportion.close.50 * B_avg.bright.female, 
#     family = "binomial", data = dat3, weights = clutch.size)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2232  -0.9894   0.5108   0.8418   1.3241  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                               15.9563     8.7111   1.832   0.0670 .
# proportion.close.50                     -152.4635    75.9991  -2.006   0.0448 *
#   B_avg.bright.female                       -0.5672     0.3216  -1.763   0.0778 .
# proportion.close.50:B_avg.bright.female    5.2449     2.7404   1.914   0.0556 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 27.223  on 9  degrees of freedom
# Residual deviance: 14.490  on 6  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 32.109
# 
# Number of Fisher Scoring iterations: 5

# 4) movement + male plumage
Anova(propEP.prop50m.plumM)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep
# LR Chisq Df Pr(>Chisq)   
# proportion.close.50   8.2022  1   0.004184 **
#   B_avg.bright.male     2.1881  1   0.139081 
confint(propEP.prop50m.plumM)
# 2.5 %     97.5 %
#   (Intercept)          -3.42039720  1.6065740
# proportion.close.50 -21.95856626 -3.2499351
# B_avg.bright.male    -0.02243126  0.1911544
summary(propEP.prop50m.plumM)
# glm(formula = prop.ep ~ proportion.close.50 + B_avg.bright.male, 
#     family = "binomial", data = dat3, weights = clutch.size)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.95050  -1.00030   0.06605   1.07523   2.07090  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)          -0.73532    1.23865  -0.594   0.5527  
# proportion.close.50 -11.25092    4.61272  -2.439   0.0147 *
#   B_avg.bright.male     0.07166    0.05152   1.391   0.1642  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 27.223  on 9  degrees of freedom
# Residual deviance: 17.697  on 7  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 33.315
# 
# Number of Fisher Scoring iterations: 5

# 5) movement * male plumage
# For male plumage color, the interaction is NOT significant from the type III test 
# so report the results from the type II test instead
Anova(propEP.prop50m.plumM.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: prop.ep
# LR Chisq Df Pr(>Chisq)
# proportion.close.50                   0.166110  1     0.6836
# B_avg.bright.male                     0.005104  1     0.9430
# proportion.close.50:B_avg.bright.male 0.074943  1     0.7843
confint(propEP.prop50m.plumM.int)
# 2.5 %     97.5 %
#   (Intercept)                            -13.892020  17.298436
# proportion.close.50                   -206.794006 131.103671
# B_avg.bright.male                       -0.756651   0.687209
# proportion.close.50:B_avg.bright.male   -6.429873   8.700483
summary(propEP.prop50m.plumM.int)
# glm(formula = prop.ep ~ proportion.close.50 * B_avg.bright.male, 
#     family = "binomial", data = dat3, weights = clutch.size)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.92824  -1.04830   0.01212   1.15982   2.09912  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                             1.35660    7.76069   0.175    0.861
# proportion.close.50                   -34.22398   84.24777  -0.406    0.685
# B_avg.bright.male                      -0.02569    0.35975  -0.071    0.943
# proportion.close.50:B_avg.bright.male   1.03157    3.77326   0.273    0.785
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 27.223  on 9  degrees of freedom
# Residual deviance: 17.622  on 6  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 35.24
# 
# Number of Fisher Scoring iterations: 5



################################################################################
#-------------------------------------------------------------------------------
# Models for prop EP offspring across whole season
# USE RAW PLUMAGE VALUES
# ALSO COMPARE RESULTS WITH 11 birds vs. 10 birds (birds from replacement clutch)

# Use proportion.close.50 as movement variable

# subset data to only use the 10 birds with replacement clutches
dat.rep <- subset(dat4, !is.na(dat4$Brood_ind1))


# 2) just prop50m
propEP.season.prop50m <- glm(prop.ep.season ~ proportion.close.50, weights=num.chicks, 
                      data=dat4, family="binomial")
summary(propEP.season.prop50m)
confint(propEP.season.prop50m)

# 2.2) prop50m with smaller data set - not significant any more!
propEP.season.prop50m.rep <- glm(prop.ep.season ~ proportion.close.50, weights=num.chicks, 
                             data=dat.rep, family="binomial")
summary(propEP.season.prop50m.rep)

# 3) just female plumage
propEP.season.plumF <- glm(prop.ep.season ~ B_avg.bright.female, weights=num.chicks, 
                    data=dat4, family="binomial")
summary(propEP.season.plumF)
confint(propEP.season.plumF)

# 3.2) just female plumage with smaller data set - still significant
propEP.season.plumF.rep <- glm(prop.ep.season ~ B_avg.bright.female, weights=num.chicks, 
                           data=dat.rep, family="binomial")
summary(propEP.season.plumF.rep)
confint(propEP.season.plumF.rep)

# 4) just male plumage
propEP.season.plumM <- glm(prop.ep.season ~ B_avg.bright.male, weights=num.chicks, 
                    data=dat4, family="binomial")
summary(propEP.season.plumM)
confint(propEP.season.plumM)

# 4.2) just male plumage with smaller data set - still significant
propEP.season.plumM.rep <- glm(prop.ep.season ~ B_avg.bright.male, weights=num.chicks, 
                           data=dat.rep, family="binomial")
summary(propEP.season.plumM.rep)
confint(propEP.season.plumM.rep)

# 5) prop50m + female plumage - female plumage sig but not movement
propEP.season.prop50m.plumF <- glm(prop.ep.season ~  proportion.close.50 + B_avg.bright.female,
                            weights=num.chicks, 
                            data=dat4, family="binomial")
summary(propEP.season.prop50m.plumF)
confint(propEP.season.prop50m.plumF)

# 5.5) prop50m + female plumage, smaller data - same results as above
propEP.season.prop50m.plumF.rep <- glm(prop.ep.season ~ proportion.close.50 + B_avg.bright.female,
                                   weights=num.chicks, 
                                   data=dat.rep, family="binomial")
summary(propEP.season.prop50m.plumF.rep)
confint(propEP.season.prop50m.plumF.rep)

# 6) prop50m * female plumage - nothing significant
propEP.season.prop50m.plumF.int <- glm(prop.ep.season ~ proportion.close.50 * B_avg.bright.female, 
                                weights=num.chicks, 
                                data=dat4, family="binomial")
summary(propEP.season.prop50m.plumF.int)
confint(propEP.season.prop50m.plumF.int)

# 6.2) prop50m * female plumage, smaller data - nothing significant
propEP.season.prop50m.plumF.int.rep <- glm(prop.ep.season ~ proportion.close.50 * B_avg.bright.female, 
                                       weights=num.chicks, 
                                       data=dat.rep, family="binomial")
summary(propEP.season.prop50m.plumF.int.rep)
confint(propEP.season.prop50m.plumF.int.rep)

# 7) prop50m + male plumage - male plumage sig but not movement
propEP.season.prop50m.plumM <- glm(prop.ep.season ~ proportion.close.50 + B_avg.bright.male,
                            weights=num.chicks, 
                            data=dat4, family="binomial")
summary(propEP.season.prop50m.plumM)
confint(propEP.season.prop50m.plumM)

# 7.2) prop50m + male plumage, smaller data - male plumage still sig
propEP.season.prop50m.plumM.rep <- glm(prop.ep.season ~ proportion.close.50 + B_avg.bright.male,
                                   weights=num.chicks, 
                                   data=dat.rep, family="binomial")
summary(propEP.season.prop50m.plumM.rep)
confint(propEP.season.prop50m.plumM.rep)

# 8) prop50m * male plumage - male plumage sig but not movement or interaction
propEP.season.prop50m.plumM.int <- glm(prop.ep.season ~ proportion.close.50 * B_avg.bright.male, 
                                weights=num.chicks, 
                                data=dat4, family="binomial")
summary(propEP.season.prop50m.plumM.int)
confint(propEP.season.prop50m.plumM.int)

# 8) prop50m * male plumage, smaller data - nothing significant
propEP.season.prop50m.plumM.int.rep <- glm(prop.ep.season ~ proportion.close.50 * B_avg.bright.male, 
                                       weights=num.chicks, 
                                       data=dat.rep, family="binomial")
summary(propEP.season.prop50m.plumM.int.rep)
confint(propEP.season.prop50m.plumM.int.rep)


# Anova tests for significance -------------------------------------------------

# 1) movement only
Anova(propEP.season.prop50m)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)  
# proportion.close.50   3.3201  1    0.06844 .
confint(propEP.season.prop50m)
# 2.5 %    97.5 %
#   (Intercept)         -0.7066788 0.4038400
# proportion.close.50 -8.3993009 0.2849376
summary(propEP.season.prop50m)
# glm(formula = prop.ep.season ~ proportion.close.50, family = "binomial", 
#     data = dat4, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.7752  -1.1844  -0.5029   1.0824   2.7799  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)          -0.1481     0.2820  -0.525   0.5993  
# proportion.close.50  -3.8718     2.1954  -1.764   0.0778 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 33.655  on 10  degrees of freedom
# Residual deviance: 30.335  on  9  degrees of freedom
# AIC: 59.557
# 
# Number of Fisher Scoring iterations: 4

# 1.2) smaller data set - not significant
Anova(propEP.season.prop50m.rep)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)
# proportion.close.50   1.1039  1     0.2934
confint(propEP.season.prop50m.rep)
# 2.5 %    97.5 %
#   (Intercept)         -1.067343 0.1527339
# proportion.close.50 -6.924251 1.9704169
summary(propEP.season.prop50m.rep)
# glm(formula = prop.ep.season ~ proportion.close.50, family = "binomial", 
#     data = dat.rep, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.5957  -1.0494  -0.3883   0.5374   3.0060  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          -0.4471     0.3096  -1.444    0.149
# proportion.close.50  -2.3258     2.2486  -1.034    0.301
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 24.227  on 9  degrees of freedom
# Residual deviance: 23.123  on 8  degrees of freedom
# AIC: 50.476
# 
# Number of Fisher Scoring iterations: 4



# 2) movement + female plumage
Anova(propEP.season.prop50m.plumF)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)   
# proportion.close.50   1.5382  1   0.214886   
# B_avg.bright.female   9.8943  1   0.001658 **
confint(propEP.season.prop50m.plumF)
# 2.5 %      97.5 %
#   (Intercept)          1.5280896  8.08296471
# proportion.close.50 -7.2237068  1.52412731
# B_avg.bright.female -0.3185398 -0.06693612
summary(propEP.season.prop50m.plumF)
# glm(formula = prop.ep.season ~ proportion.close.50 + B_avg.bright.female, 
#     family = "binomial", data = dat4, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.9870  -0.6818  -0.4967   0.2812   2.6413  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)          4.61810    1.65673   2.787  0.00531 **
#   proportion.close.50 -2.68794    2.21061  -1.216  0.22401   
# B_avg.bright.female -0.18514    0.06362  -2.910  0.00361 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 33.655  on 10  degrees of freedom
# Residual deviance: 20.441  on  8  degrees of freedom
# AIC: 51.663
# 
# Number of Fisher Scoring iterations: 4


# 2.2) movement + female plumage, smaller data set
# Anova(propEP.season.prop50m.plumF.rep)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)   
# proportion.close.50   0.2656  1   0.606322   
# B_avg.bright.female   9.7195  1   0.001823 **
confint(propEP.season.prop50m.plumF.rep)
# 2.5 %      97.5 %
#   (Intercept)          1.2028329  7.80083332
# proportion.close.50 -5.7844910  3.20487386
# B_avg.bright.female -0.3203844 -0.06622128
summary(propEP.season.prop50m.plumF.rep)
# glm(formula = prop.ep.season ~ proportion.close.50 + B_avg.bright.female, 
#     family = "binomial", data = dat.rep, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.8330  -0.5587  -0.2497  -0.1741   2.8857  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)          4.32628    1.66799   2.594  0.00949 **
#   proportion.close.50 -1.16305    2.27213  -0.512  0.60874   
# B_avg.bright.female -0.18579    0.06429  -2.890  0.00386 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 24.227  on 9  degrees of freedom
# Residual deviance: 13.403  on 7  degrees of freedom
# AIC: 42.756
# 
# Number of Fisher Scoring iterations: 4


# 3) movement * female plumage
Anova(propEP.season.prop50m.plumF.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)
# proportion.close.50                      0.10279  1     0.7485
# B_avg.bright.female                      2.62252  1     0.1054
# proportion.close.50:B_avg.bright.female  0.05006  1     0.8230
confint(propEP.season.prop50m.plumF.int)
# 2.5 %      97.5 %
#   (Intercept)                              -1.3134425 12.31271229
# proportion.close.50                     -64.4393381 45.20195618
# B_avg.bright.female                      -0.4819261  0.04380305
# proportion.close.50:B_avg.bright.female  -1.8286824  2.33147727
summary(propEP.season.prop50m.plumF.int)
# glm(formula = prop.ep.season ~ proportion.close.50 * B_avg.bright.female, 
#     family = "binomial", data = dat4, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.9560  -0.6942  -0.5872   0.2957   2.6553  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                               5.2902     3.4412   1.537    0.124
# proportion.close.50                      -8.8406    27.6295  -0.320    0.749
# B_avg.bright.female                      -0.2111     0.1328  -1.589    0.112
# proportion.close.50:B_avg.bright.female   0.2344     1.0485   0.224    0.823
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 33.655  on 10  degrees of freedom
# Residual deviance: 20.391  on  7  degrees of freedom
# AIC: 53.613
# 
# Number of Fisher Scoring iterations: 4


# 3.2) movement * female plumage, smaller sample size
Anova(propEP.season.prop50m.plumF.int.rep, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)
# proportion.close.50                      0.04649  1     0.8293
# B_avg.bright.female                      2.46673  1     0.1163
# proportion.close.50:B_avg.bright.female  0.03024  1     0.8619
confint(propEP.season.prop50m.plumF.int.rep)
# 2.5 %     97.5 %
#   (Intercept)                              -1.7969911 11.8144179
# proportion.close.50                     -60.9887225 48.2091970
# B_avg.bright.female                      -0.4758775  0.0506299
# proportion.close.50:B_avg.bright.female  -1.8837482  2.2629751
summary(propEP.season.prop50m.plumF.int.rep)
# glm(formula = prop.ep.season ~ proportion.close.50 * B_avg.bright.female, 
#     family = "binomial", data = dat.rep, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.8094  -0.5581  -0.2203  -0.1673   2.8978  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                               4.8477     3.4366   1.411    0.158
# proportion.close.50                      -5.9304    27.5253  -0.215    0.829
# B_avg.bright.female                      -0.2060     0.1331  -1.548    0.122
# proportion.close.50:B_avg.bright.female   0.1818     1.0456   0.174    0.862
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 24.227  on 9  degrees of freedom
# Residual deviance: 13.373  on 6  degrees of freedom
# AIC: 44.726
# 
# Number of Fisher Scoring iterations: 4



# 4) movement + male plumage
Anova(propEP.season.prop50m.plumM)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)    
# proportion.close.50   2.3823  1     0.1227    
# B_avg.bright.male    15.3266  1  9.043e-05 ***

confint(propEP.season.prop50m.plumM)
# 2.5 %     97.5 %
#   (Intercept)         -4.50029888 -1.3285205
# proportion.close.50 -8.08537721  0.8939015
# B_avg.bright.male    0.05257088  0.1775882
summary(propEP.season.prop50m.plumM)
# glm(formula = prop.ep.season ~ proportion.close.50 + B_avg.bright.male, 
#     family = "binomial", data = dat4, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.9867  -0.7175   0.1234   0.4142   2.7653  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -2.82954    0.80064  -3.534 0.000409 ***
#   proportion.close.50 -3.40452    2.26841  -1.501 0.133397    
# B_avg.bright.male    0.10982    0.03126   3.512 0.000444 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 33.655  on 10  degrees of freedom
# Residual deviance: 15.008  on  8  degrees of freedom
# AIC: 46.231
# 
# Number of Fisher Scoring iterations: 4


# 4.2) movement + male plumage, smaller sample size
Anova(propEP.season.prop50m.plumM.rep)
# Analysis of Deviance Table (Type II tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)   
# proportion.close.50   1.8082  1   0.178722   
# B_avg.bright.male     8.4824  1   0.003586 **
confint(propEP.season.prop50m.plumM.rep)
# 2.5 %     97.5 %
#   (Intercept)         -4.433098 -1.0083687
# proportion.close.50 -7.860059  1.3678872
# B_avg.bright.male    0.031016  0.1769565
summary(propEP.season.prop50m.plumM.rep)
# glm(formula = prop.ep.season ~ proportion.close.50 + B_avg.bright.male, 
#     family = "binomial", data = dat.rep, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.0239  -0.6775  -0.1761   0.4388   2.6949  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)         -2.61172    0.85841  -3.042  0.00235 **
#   proportion.close.50 -3.06403    2.33097  -1.314  0.18868   
# B_avg.bright.male    0.09763    0.03618   2.698  0.00697 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 24.227  on 9  degrees of freedom
# Residual deviance: 14.640  on 7  degrees of freedom
# AIC: 43.993
# 
# Number of Fisher Scoring iterations: 4


# 5) movement * male plumage
Anova(propEP.season.prop50m.plumM.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)   
# proportion.close.50                     0.3464  1    0.55615   
# B_avg.bright.male                       7.4567  1    0.00632 **
#   proportion.close.50:B_avg.bright.male   0.6202  1    0.43097  
confint(propEP.season.prop50m.plumM.int)
# 2.5 %     97.5 %
#   (Intercept)                            -7.78377211 -1.0748695
# proportion.close.50                   -23.68472793 55.0317350
# B_avg.bright.male                       0.03868799  0.3306405
# proportion.close.50:B_avg.bright.male  -2.56920549  0.8642876
summary(propEP.season.prop50m.plumM.int)
# glm(formula = prop.ep.season ~ proportion.close.50 * B_avg.bright.male, 
#     family = "binomial", data = dat4, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.0447  -0.6933  -0.1181   0.4557   2.5816  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                           -3.82217    1.60346  -2.384   0.0171 *
#   proportion.close.50                   10.77564   18.88694   0.571   0.5683  
# B_avg.bright.male                      0.15334    0.06923   2.215   0.0268 *
#   proportion.close.50:B_avg.bright.male -0.61765    0.82012  -0.753   0.4514  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 33.655  on 10  degrees of freedom
# Residual deviance: 14.388  on  7  degrees of freedom
# AIC: 47.611
# 
# Number of Fisher Scoring iterations: 4


# 5.2) movement * male plumage, smaller sample size
Anova(propEP.season.prop50m.plumM.int.rep, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: prop.ep.season
# LR Chisq Df Pr(>Chisq)
# proportion.close.50                    0.33211  1     0.5644
# B_avg.bright.male                      1.16523  1     0.2804
# proportion.close.50:B_avg.bright.male  0.40598  1     0.5240
confint(propEP.season.prop50m.plumM.int.rep)
# 2.5 %      97.5 %
#   (Intercept)                           -15.3913575   3.5739511
# proportion.close.50                   -69.8823805 134.3396566
# B_avg.bright.male                      -0.1889884   0.6852987
# proportion.close.50:B_avg.bright.male  -6.1507705   2.9941346
summary(propEP.season.prop50m.plumM.int.rep)
# glm(formula = prop.ep.season ~ proportion.close.50 * B_avg.bright.male, 
#     family = "binomial", data = dat.rep, weights = num.chicks)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.0535  -0.8230  -0.1090   0.5552   2.4627  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                            -5.5809     4.7864  -1.166    0.244
# proportion.close.50                    29.5228    51.5809   0.572    0.567
# B_avg.bright.male                       0.2352     0.2208   1.065    0.287
# proportion.close.50:B_avg.bright.male  -1.4615     2.3103  -0.633    0.527
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 24.227  on 9  degrees of freedom
# Residual deviance: 14.234  on 6  degrees of freedom
# AIC: 45.587
# 
# Number of Fisher Scoring iterations: 4


################################################################################
#-------------------------------------------------------------------------------


# Significance test for number of EP mates (with number of chicks accounted for) 
# across season
# movement component is proportion.close.50 (again)

# Complete null

numEP.null <- glm(num.ep.mates ~ 1, family="poisson", data=dat4)

# 1) Null
numEP.numchick <- glm(num.ep.mates ~ num.chicks, family = "poisson", data = dat4)
summary(numEP.numchick)

# 2) just close50
numEP.close50 <- glm(num.ep.mates ~ num.chicks + proportion.close.50, family = "poisson", data = dat4)
summary(numEP.close50)

# 2.2) just 50% KDE, smaller samples size
numEP.close50.rep <- glm(num.ep.mates ~ num.chicks + proportion.close.50, family = "poisson", data = dat.rep)

# 3) just female plumage
numEP.plumF <- glm(num.ep.mates ~ num.chicks + B_avg.bright.female,
                   family = "poisson", data = dat4)
summary(numEP.plumF)

# 3.2) just female plumage, smaller sample size
numEP.plumF.rep <- glm(num.ep.mates ~ num.chicks + B_avg.bright.female,
                   family = "poisson", data = dat.rep)

# 4) just male plumage
numEP.plumM <- glm(num.ep.mates ~ num.chicks + B_avg.bright.male,
                   family = "poisson", data = dat4)
summary(numEP.plumM)

# 4.2) just male plumage, smaller sample size
numEP.plumM.rep <- glm(num.ep.mates ~ num.chicks + B_avg.bright.male,
                   family = "poisson", data = dat.rep)

# 5) close50 + female plumage
numEP.close50.plumF <- glm(num.ep.mates ~ num.chicks + proportion.close.50 + 
                            B_avg.bright.female,
                          family = "poisson", data = dat4)
summary(numEP.close50.plumF)

# 5.2) close50 + female plumage, smaller sample size
numEP.close50.plumF.rep <- glm(num.ep.mates ~ num.chicks + proportion.close.50 + 
                           B_avg.bright.female,
                         family = "poisson", data = dat.rep)

# 6) close50 * female plumage
numEP.close50.plumF.int <- glm(num.ep.mates ~ num.chicks + proportion.close.50 * 
                                B_avg.bright.female,
                              family = "poisson", data = dat4)
summary(numEP.close50.plumF.int)

# 6.2) close50 * female plumage, smaller sample size
numEP.close50.plumF.int.rep <- glm(num.ep.mates ~ num.chicks + proportion.close.50 * 
                               B_avg.bright.female,
                             family = "poisson", data = dat.rep)

# 7) close50 + male plumage
numEP.close50.plumM <- glm(num.ep.mates ~ num.chicks + proportion.close.50 + 
                            B_avg.bright.male,
                          family = "poisson", data = dat4)
summary(numEP.close50.plumM)

# 7.2) close50 + male plumage, smaller sample size
numEP.close50.plumM.rep <- glm(num.ep.mates ~ num.chicks + proportion.close.50 + 
                           B_avg.bright.male,
                         family = "poisson", data = dat.rep)

# 8) close50 * male plumage
numEP.close50.plumM.int <- glm(num.ep.mates ~ num.chicks + proportion.close.50 *
                                B_avg.bright.male,
                              family = "poisson", data = dat4)
summary(numEP.close50.plumM.int)

# 8.2) close50 * male plumage, smaller sample size
numEP.close50.plumM.int.rep <- glm(num.ep.mates ~ num.chicks + 
                                     proportion.close.50 *
                               B_avg.bright.male,
                             family = "poisson", data = dat.rep)

################### Likelihood ratio tests for number of EP mates in whole season

# 1) movement only
Anova(numEP.close50)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks           0.01117  1     0.9158
# proportion.close.50  0.77540  1     0.3786

confint(numEP.close50)
# 2.5 %    97.5 %
#   (Intercept)         -1.375023 2.3893114
# num.chicks          -0.181923 0.1930371
# proportion.close.50 -8.812789 2.7996036
summary(numEP.close50)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50, 
#     family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.78782  -0.62129  -0.03686   0.58686   0.97981  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.600478   0.946908   0.634    0.526
# num.chicks           0.009982   0.094323   0.106    0.916
# proportion.close.50 -2.460913   2.907889  -0.846    0.397
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.5910  on 10  degrees of freedom
# Residual deviance: 6.7839  on  8  degrees of freedom
# AIC: 36.987
# 
# Number of Fisher Scoring iterations: 5



# 1.2) movement only, smaller sample size
Anova(numEP.close50.rep)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks           0.09310  1     0.7603
# proportion.close.50  0.44044  1     0.5069
confint(numEP.close50.rep)
# 2.5 %    97.5 %
#   (Intercept)         -2.0745848 2.2374374
# num.chicks          -0.1705865 0.2346101
# proportion.close.50 -8.3774215 3.5209150
summary(numEP.close50.rep)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50, 
#     family = "poisson", data = dat.rep)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.68185  -0.52222  -0.08305   0.37130   1.02344  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.26437    1.08077   0.245    0.807
# num.chicks           0.03093    0.10139   0.305    0.760
# proportion.close.50 -1.92209    2.97636  -0.646    0.518
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 6.5646  on 9  degrees of freedom
# Residual deviance: 6.1008  on 7  degrees of freedom
# AIC: 33.312
# 
# Number of Fisher Scoring iterations: 5



# 2) movement + female plumage
Anova(numEP.close50.plumF)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks           0.11289  1     0.7369
# proportion.close.50  0.51460  1     0.4732
# B_avg.bright.female  1.31278  1     0.2519
confint(numEP.close50.plumF)
# 2.5 %     97.5 %
#   (Intercept)         -1.2372989 6.14898477
# num.chicks          -0.1602915 0.21580770
# proportion.close.50 -8.3169561 3.18485766
# B_avg.bright.female -0.2233748 0.05665941
summary(numEP.close50.plumF)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 + 
#       B_avg.bright.female, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.44705  -0.51751  -0.03027   0.32588   1.17510  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          2.43499    1.87071   1.302    0.193
# num.chicks           0.03189    0.09451   0.337    0.736
# proportion.close.50 -1.99494    2.87953  -0.693    0.488
# B_avg.bright.female -0.08088    0.07129  -1.135    0.257
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.5910  on 10  degrees of freedom
# Residual deviance: 5.4712  on  7  degrees of freedom
# AIC: 37.674
# 
# Number of Fisher Scoring iterations: 5


# 2.2) movement + female plumage, smaller sample size
Anova(numEP.close50.plumF.rep)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks           0.30672  1     0.5797
# proportion.close.50  0.22758  1     0.6333
# B_avg.bright.female  1.49428  1     0.2216
confint(numEP.close50.plumF.rep)
# 2.5 %     97.5 %
#   (Intercept)         -1.6454073 6.12617975
# num.chicks          -0.1460100 0.26219810
# proportion.close.50 -7.8116364 4.00399251
# B_avg.bright.female -0.2376057 0.05265438
summary(numEP.close50.plumF.rep)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 + 
#       B_avg.bright.female, family = "poisson", data = dat.rep)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.3224  -0.4332  -0.1665   0.1470   1.2499  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          2.25326    1.96056   1.149    0.250
# num.chicks           0.05649    0.10204   0.554    0.580
# proportion.close.50 -1.38013    2.95412  -0.467    0.640
# B_avg.bright.female -0.08899    0.07383  -1.205    0.228
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 6.5646  on 9  degrees of freedom
# Residual deviance: 4.6065  on 6  degrees of freedom
# AIC: 33.818
# 
# Number of Fisher Scoring iterations: 5


# 3) movement * female plumage
Anova(numEP.close50.plumF.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)  
# num.chicks                                0.5217  1     0.4701  
# proportion.close.50                       1.7412  1     0.1870  
# B_avg.bright.female                       2.7362  1     0.0981 .
# proportion.close.50:B_avg.bright.female   1.5928  1     0.2069  

confint(numEP.close50.plumF.int)
# 2.5 %      97.5 %
#   (Intercept)                               -0.7525180 13.36627754
# num.chicks                                -0.1322801  0.29162454
# proportion.close.50                     -115.1660369 20.45497095
# B_avg.bright.female                       -0.5782366  0.04525841
# proportion.close.50:B_avg.bright.female   -0.8600075  4.21635346
summary(numEP.close50.plumF.int)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 * 
#       B_avg.bright.female, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.1800  -0.4073  -0.1418   0.3207   1.0335  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                               6.25494    3.53847   1.768   0.0771 .
# num.chicks                                0.07644    0.10619   0.720   0.4716  
# proportion.close.50                     -42.88701   33.60851  -1.276   0.2019  
# B_avg.bright.female                      -0.24940    0.15553  -1.604   0.1088  
# proportion.close.50:B_avg.bright.female   1.55372    1.25636   1.237   0.2162  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.5910  on 10  degrees of freedom
# Residual deviance: 3.8783  on  6  degrees of freedom
# AIC: 38.081
# 
# Number of Fisher Scoring iterations: 5


# 3.2) movement * female plumage, smaller sample size
Anova(numEP.close50.plumF.int.rep, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)  
# num.chicks                                1.0215  1    0.31216  
# proportion.close.50                       2.2837  1    0.13074  
# B_avg.bright.female                       3.5001  1    0.06136 .
# proportion.close.50:B_avg.bright.female   2.1777  1    0.14002  
confint(numEP.close50.plumF.int.rep)
# 2.5 %     97.5 %
#   (Intercept)                               -0.4236404 15.4296696
# num.chicks                                -0.1083780  0.3677656
# proportion.close.50                     -138.1065145 14.9941055
# B_avg.bright.female                       -0.7024257  0.0140380
# proportion.close.50:B_avg.bright.female   -0.6235464  5.0965104
summary(numEP.close50.plumF.int.rep)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 * 
#       B_avg.bright.female, family = "poisson", data = dat.rep)
# 
# Deviance Residuals: 
#   1         2         3         4         5         6         7         8         9  
# -0.13404   0.60993  -0.03659   0.12034  -0.31622   0.02848  -0.98087   0.68969   0.18253  
# 10  
# -0.67164  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                               7.1738     3.9227   1.829   0.0674 .
# num.chicks                                0.1172     0.1184   0.990   0.3223  
# proportion.close.50                     -52.9866    37.5054  -1.413   0.1577  
# B_avg.bright.female                      -0.3086     0.1765  -1.749   0.0803 .
# proportion.close.50:B_avg.bright.female   1.9588     1.3994   1.400   0.1616  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 6.5646  on 9  degrees of freedom
# Residual deviance: 2.4288  on 5  degrees of freedom
# AIC: 33.64
# 
# Number of Fisher Scoring iterations: 5


# 4) movement + male plumage
Anova(numEP.close50.plumM)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks           0.04115  1     0.8392
# proportion.close.50  0.74887  1     0.3868
# B_avg.bright.male    0.22268  1     0.6370
confint(numEP.close50.plumM)
# 2.5 %     97.5 %
#   (Intercept)         -2.75538670 2.67827989
# num.chicks          -0.17725742 0.21543764
# proportion.close.50 -8.83246523 2.86906137
# B_avg.bright.male   -0.04769029 0.06821057
summary(numEP.close50.plumM)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 + 
#       B_avg.bright.male, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.67284  -0.66835  -0.00847   0.45612   1.10515  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.15191    1.36476   0.111    0.911
# num.chicks           0.02000    0.09850   0.203    0.839
# proportion.close.50 -2.43862    2.92851  -0.833    0.405
# B_avg.bright.male    0.01388    0.02893   0.480    0.631
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.5910  on 10  degrees of freedom
# Residual deviance: 6.5613  on  7  degrees of freedom
# AIC: 38.764
# 
# Number of Fisher Scoring iterations: 5



# 4.2) movement + male plumage, smaller sample size
Anova(numEP.close50.plumM.rep)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks           0.08743  1     0.7675
# proportion.close.50  0.41332  1     0.5203
# B_avg.bright.male    0.01114  1     0.9159
confint(numEP.close50.plumM.rep)
# 2.5 %     97.5 %
#   (Intercept)         -2.6886398 3.26928524
# num.chicks          -0.1713643 -0.17474168
# proportion.close.50 -8.4126276 3.59623553
# B_avg.bright.male   -0.1053488 0.06711351
summary(numEP.close50.plumM.rep)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 + 
#       B_avg.bright.male, family = "poisson", data = dat.rep)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.7053  -0.5100  -0.0897   0.3737   1.0178  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.370206   1.474705   0.251    0.802
# num.chicks           0.030012   0.101567   0.295    0.768
# proportion.close.50 -1.876473   3.000182  -0.625    0.532
# B_avg.bright.male   -0.004435   0.042426  -0.105    0.917
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 6.5646  on 9  degrees of freedom
# Residual deviance: 6.0897  on 6  degrees of freedom
# AIC: 35.301
# 
# Number of Fisher Scoring iterations: 5



# 5) movement * male plumage
Anova(numEP.close50.plumM.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks                             0.11511  1     0.7344
# proportion.close.50                    0.73952  1     0.3898
# B_avg.bright.male                      1.20914  1     0.2715
# proportion.close.50:B_avg.bright.male  1.07608  1     0.2996
confint(numEP.close50.plumM.int)
# 2.5 %     97.5 %
#   (Intercept)                            -4.54799876  2.3266513
# num.chicks                             -0.16726847  0.2392844
# proportion.close.50                   -16.86991187 49.7681679
# B_avg.bright.male                      -0.04257601  0.1307053
# proportion.close.50:B_avg.bright.male  -2.21716390  0.5899022
summary(numEP.close50.plumM.int)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 * 
#       B_avg.bright.male, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.73943  -0.30141  -0.07708   0.25555   0.98009  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                           -0.87150    1.71646  -0.508    0.612
# num.chicks                             0.03450    0.10175   0.339    0.735
# proportion.close.50                   13.63378   16.36175   0.833    0.405
# B_avg.bright.male                      0.04907    0.04279   1.147    0.251
# proportion.close.50:B_avg.bright.male -0.68310    0.68579  -0.996    0.319
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.5910  on 10  degrees of freedom
# Residual deviance: 5.4852  on  6  degrees of freedom
# AIC: 39.688
# 
# Number of Fisher Scoring iterations: 5


# 5.2) movement * male plumage, smaller sample size
Anova(numEP.close50.plumM.int.rep, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: num.ep.mates
# LR Chisq Df Pr(>Chisq)
# num.chicks                             0.15094  1     0.6976
# proportion.close.50                    2.13788  1     0.1437
# B_avg.bright.male                      2.11842  1     0.1455
# proportion.close.50:B_avg.bright.male  2.23919  1     0.1346
confint(numEP.close50.plumM.int.rep)
# 2.5 %      97.5 %
#   (Intercept)                           -23.4076210   3.1592723
# num.chicks                             -0.1584412   0.2345329
# proportion.close.50                   -28.6470674 241.6599598
# B_avg.bright.male                      -0.1240059   1.0388904
# proportion.close.50:B_avg.bright.male -10.8216584   1.1951966
summary(numEP.close50.plumM.int.rep)
# glm(formula = num.ep.mates ~ num.chicks + proportion.close.50 * 
#       B_avg.bright.male, family = "poisson", data = dat.rep)
# 
# Deviance Residuals: 
#   1         2         3         4         5         6         7         8         9  
# -0.21884   0.42440  -0.24872   0.28479   0.32436  -0.03905  -1.69160   0.41450   0.49081  
# 10  
# -0.31393  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                           -8.25146    6.55982  -1.258    0.208
# num.chicks                             0.03826    0.09836   0.389    0.697
# proportion.close.50                   90.45014   66.67112   1.357    0.175
# B_avg.bright.male                      0.38581    0.28636   1.347    0.178
# proportion.close.50:B_avg.bright.male -4.12261    2.96238  -1.392    0.164
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 6.5646  on 9  degrees of freedom
# Residual deviance: 3.8505  on 5  degrees of freedom
# AIC: 35.062
# 
# Number of Fisher Scoring iterations: 5



################################################################################
#-------------------------------------------------------------------------------
# Models for number of EP mates in replacement clutch

# Best movement variable from AICc is proportion.close.50

# 0) complete null
mates.rep.null <- glm(mates.replacement.ep ~ 1, data=dat4, family="poisson")

# 1) Null with number of chicks
mates.rep.cov <- glm(mates.replacement.ep ~ chicks.replacement, 
                     data=dat4, family="poisson")

# 2) just close50
mates.rep.cov.close50 <- glm(mates.replacement.ep ~ chicks.replacement +
                       proportion.close.50, 
                     data=dat4, family="poisson")
confint(mates.rep.cov.close50)



# 3) just female plumage
mates.rep.cov.plumF <- glm(mates.replacement.ep ~ chicks.replacement +
                             B_avg.bright.female, 
                           data=dat4, family="poisson")

# 4) just male plumage
mates.rep.cov.plumM <- glm(mates.replacement.ep ~ chicks.replacement +
                             B_avg.bright.male, 
                           data=dat4, family="poisson")
summary(mates.rep.cov.plumM)
confint(mates.rep.cov.plumM)

# 5) close50 + female plumage
mates.rep.cov.plumF.close50 <- glm(mates.replacement.ep ~ chicks.replacement +
                             B_avg.bright.female +
                             proportion.close.50, 
                           data=dat4, family="poisson")

# 6) close50 * female plumage
mates.rep.cov.plumF.close50.int <- glm(mates.replacement.ep ~ chicks.replacement +
                                   B_avg.bright.female * proportion.close.50, 
                                 data=dat4, family="poisson")

# 7) close50 + male plumage
mates.rep.cov.plumM.close50 <- glm(mates.replacement.ep ~ chicks.replacement +
                                   B_avg.bright.male +
                                   proportion.close.50, 
                                 data=dat4, family="poisson")

# 8) close50 * male plumage
mates.rep.cov.plumM.close50.int <- glm(mates.replacement.ep ~ chicks.replacement +
                                       B_avg.bright.male * proportion.close.50, 
                                     data=dat4, family="poisson")


################### Likelihood ratio tests for number EP mates in replacement---

#1) movement only
Anova(mates.rep.cov.close50)
# Analysis of Deviance Table (Type II tests)
# 
# Response: mates.replacement.ep
# LR Chisq Df Pr(>Chisq)
# chicks.replacement   0.26127  1     0.6092
# proportion.close.50  0.30561  1     0.5804
confint(mates.rep.cov.close50)
# 2.5 %   97.5 %
#   (Intercept)          -4.8512422 1.826204
# chicks.replacement   -0.5523493 1.059064
# proportion.close.50 -12.6375908 5.534233
summary(mates.rep.cov.close50)
# glm(formula = mates.replacement.ep ~ chicks.replacement + proportion.close.50, 
#     family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.19220  -1.14247   0.09126   0.59566   0.90043  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          -0.9352     1.6376  -0.571    0.568
# chicks.replacement    0.1979     0.3930   0.504    0.615
# proportion.close.50  -2.3641     4.4287  -0.534    0.593
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.766  on 9  degrees of freedom
# Residual deviance: 7.313  on 7  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 25.927
# 
# Number of Fisher Scoring iterations: 5


# 2) movement + female plumage
Anova(mates.rep.cov.plumF.close50)
# Analysis of Deviance Table (Type II tests)
# 
# Response: mates.replacement.ep
# LR Chisq Df Pr(>Chisq)
# chicks.replacement   0.27349  1     0.6010
# B_avg.bright.female  0.39607  1     0.5291
# proportion.close.50  0.42309  1     0.5154
confint(mates.rep.cov.plumF.close50)
# 2.5 %    97.5 %
#   (Intercept)         -11.1035996 3.2230725
# chicks.replacement   -0.5478374 1.0975535
# B_avg.bright.female  -0.1459485 0.2886962
# proportion.close.50 -13.8488853 5.2445909
summary(mates.rep.cov.plumF.close50)
# glm(formula = mates.replacement.ep ~ chicks.replacement + B_avg.bright.female + 
#       proportion.close.50, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.34476  -0.93657   0.04559   0.59879   0.92580  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         -2.71350    3.42956  -0.791    0.429
# chicks.replacement   0.20463    0.39888   0.513    0.608
# B_avg.bright.female  0.06676    0.10615   0.629    0.529
# proportion.close.50 -2.87199    4.62757  -0.621    0.535
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.7660  on 9  degrees of freedom
# Residual deviance: 6.9169  on 6  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 27.531
# 
# Number of Fisher Scoring iterations: 5


# 3) movement * female plumage
Anova(mates.rep.cov.plumF.close50.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: mates.replacement.ep
# LR Chisq Df Pr(>Chisq)  
# chicks.replacement                        4.8351  1    0.02789 *
#   B_avg.bright.female                       3.9989  1    0.04553 *
#   proportion.close.50                       5.2331  1    0.02216 *
#   B_avg.bright.female:proportion.close.50   5.1004  1    0.02392 *
confint(mates.rep.cov.plumF.close50.int)
# 2.5 %       97.5 %
#   (Intercept)                               -2.2421953  48.04182150
# chicks.replacement                         0.1744703   3.99488422
# B_avg.bright.female                       -2.5110010  -0.01604308
# proportion.close.50                     -650.2479345 -29.96887932
# B_avg.bright.female:proportion.close.50    1.0229300  24.02411934
summary(mates.rep.cov.plumF.close50.int)
# glm(formula = mates.replacement.ep ~ chicks.replacement + B_avg.bright.female * 
#       proportion.close.50, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   1         2         3         4         5         6         7         8         9  
# -0.04947  -0.15049  -0.21805   0.31665  -0.64266   0.74746  -0.48169  -0.16465  -0.60289  
# 10  
# 0.22205  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                               17.4045    11.9406   1.458   0.1450  
# chicks.replacement                         1.7378     0.9107   1.908   0.0564 .
# B_avg.bright.female                       -0.9560     0.5932  -1.612   0.1070  
# proportion.close.50                     -260.5351   148.3449  -1.756   0.0790 .
# B_avg.bright.female:proportion.close.50    9.5568     5.4862   1.742   0.0815 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.7660  on 9  degrees of freedom
# Residual deviance: 1.8165  on 5  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 24.43
# 
# Number of Fisher Scoring iterations: 5


# 4) movement + male plumage
Anova(mates.rep.cov.plumM.close50)
# Analysis of Deviance Table (Type II tests)
# 
# Response: mates.replacement.ep
# LR Chisq Df Pr(>Chisq)
# chicks.replacement   0.08381  1     0.7722
# B_avg.bright.male    0.46955  1     0.4932
# proportion.close.50  0.32648  1     0.5677
confint(mates.rep.cov.plumM.close50)
# 2.5 %    97.5 %
#   (Intercept)          -5.26466601 1.8537517
# chicks.replacement   -0.66380859 0.9998989
# B_avg.bright.male    -0.07799313 0.1292006
# proportion.close.50 -13.48962018 5.6749567
summary(mates.rep.cov.plumM.close50)
# glm(formula = mates.replacement.ep ~ chicks.replacement + B_avg.bright.male + 
#       proportion.close.50, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.1696  -1.0313  -0.0741   0.5260   1.0250  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         -1.42758    1.68996  -0.845    0.398
# chicks.replacement   0.11631    0.40505   0.287    0.774
# B_avg.bright.male    0.03546    0.04966   0.714    0.475
# proportion.close.50 -2.54282    4.62653  -0.550    0.583
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.7660  on 9  degrees of freedom
# Residual deviance: 6.8434  on 6  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 27.457
# 
# Number of Fisher Scoring iterations: 5


# 5) movement * male plumage
Anova(mates.rep.cov.plumM.close50.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: mates.replacement.ep
# LR Chisq Df Pr(>Chisq)
# chicks.replacement                      0.4935  1     0.4824
# B_avg.bright.male                       1.3967  1     0.2373
# proportion.close.50                     1.1563  1     0.2822
# B_avg.bright.male:proportion.close.50   1.2159  1     0.2702
confint(mates.rep.cov.plumM.close50.int)
# 2.5 %     97.5 %
#   (Intercept)                           -44.8610149   6.306670
# chicks.replacement                     -0.5362225   1.246938
# B_avg.bright.male                      -0.3008268   1.873003
# proportion.close.50                   -82.3164433 430.151140
# B_avg.bright.male:proportion.close.50 -19.1874727   3.585031
summary(mates.rep.cov.plumM.close50.int)
# glm(formula = mates.replacement.ep ~ chicks.replacement + B_avg.bright.male * 
#       proportion.close.50, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   1         2         3         4         5         6         7         8         9  
# -0.33528  -0.80436   0.51326   0.07496  -0.62431   0.71112  -1.07358   0.69773  -1.29761  
# 10  
# 0.61681  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                           -12.9917    12.1325  -1.071    0.284
# chicks.replacement                      0.2973     0.4299   0.692    0.489
# B_avg.bright.male                       0.5310     0.5112   1.039    0.299
# proportion.close.50                   115.9156   120.9703   0.958    0.338
# B_avg.bright.male:proportion.close.50  -5.2877     5.3693  -0.985    0.325
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.7660  on 9  degrees of freedom
# Residual deviance: 5.6276  on 5  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 28.241
# 
# Number of Fisher Scoring iterations: 5


################################################################################
#-------------------------------------------------------------------------------

# Visualize predictions from top model for proportion EP young in whole season

library(jtools)

# https://jtools.jacob-long.com/articles/summ.html 

# Manuscript Figure S7A
# effect of male plumage on propEP.season
effect_plot(propEP.season.prop50m.plumM, pred=B_avg.bright.male, interval=T,
            y.label="Proportion of EP young in whole season",
            x.label="Male plumage color (average belly brightness)",
            plot.points=T, int.type="confidence",
            main.title="Predicted effect of male color on \nproportion of extra-pair young in the season \nafter accounting for female movement (N=11)")
ggsave("generated-files/FigS7A_modeled effect of male plumage on propEP in whole season N=11_v5.png", h=5, w=6)
ggsave("generated-files/FigS7A_modeled effect of male plumage on propEP in whole season N=11_v5.pdf", h=5, w=6)

# Manuscript Figure S7B
effect_plot(propEP.season.prop50m.plumM.rep, pred=B_avg.bright.male, interval=T,
            y.label="Proportion of EP young in whole season",
            x.label="Male plumage color (average belly brightness)",
            plot.points=T, int.type="confidence",
            main.title="Predicted effect of male color on \nproportion of extra-pair young in the season \nafter accounting for female movement (N=10)")
ggsave("generated-files/FigS7B_modeled effect of male plumage on propEP in whole season N=10_v5.png", h=5, w=6)
ggsave("generated-files/FigS7B_modeled effect of male plumage on propEP in whole season N=10_v5.pdf", h=5, w=6)


# Manuscript Figure S6B
# effect of female plumage on season EP with n=10
effect_plot(propEP.season.prop50m.plumF.rep, pred=B_avg.bright.female, interval=T,
            y.label="Proportion of EP young in whole season",
            x.label="Female plumage color (average belly brightness)",
            plot.points=T, int.type="confidence",
            main.title="Predicted effect of female color on \nproportion of extra-pair young in the season \nafter accounting for female movement (N=10)")
ggsave("generated-files/FigS6B_modeled effect of female plumage on propEP in whole season N=10_v5.png", h=5, w=6)
ggsave("generated-files/FigS6B_modeled effect of female plumage on propEP in whole season N=10_v5.pdf", h=5, w=6)

# Manuscript Figure S6A
# effect of female plumage on season EP with N=11
effect_plot(propEP.season.prop50m.plumF, pred=B_avg.bright.female, interval=T,
            y.label="Proportion of EP young in whole season",
            x.label="Female plumage color (average belly brightness)",
            plot.points=T, int.type="confidence",
            main.title="Predicted effect of female color on \nproportion of extra-pair young in the season \nafter accounting for female movement (N=11)")
ggsave("geneated-files/FigS6A_modeled effect of female plumage on propEP in whole season N=11_v5.png", h=5, w=6)
ggsave("generated-files/FigS6A_modeled effect of female plumage on propEP in whole season N=11_v5.pdf", h=5, w=6)




#-------------------------------------------------------------------------------
# Plots for proprotion EP young in replacement

# Manuscript Figure 2A
# top model for propEP, only has proportion.close.50
effect_plot(propEP.prop50m, pred=proportion.close.50, interval=T,
            y.label="Proportion of EP young in replacement clutch",
            x.label="Proportion of GPS points within 50m of barn",
            plot.points=T, int.type="confidence",
            main.title="Predicted effect of proportion of close points on \nproportion of extra-pair young")

ggsave("generated-files/Fig2A_modeled effect of points 50m on propEP_v5.png", h=5, w=6)
ggsave("generated-files/Fig2A_modeled effect of points 50m on propEP_v5.pdf", h=5, w=6)


### Visualize effect of interaction movement with female plumage

library(interactions)

# http://cran.nexr.com/web/packages/jtools/vignettes/interactions.html 

# interact_plot now part of the "interactions" library

# Try getting predictions to plot myself

library(effects)
?effect()

# make table for predicted values
pred_table <- effect(term="proportion.close.50:B_avg.bright.female", 
                     mod=propEP.prop50m.plumF.int, x.var="proportion.close.50",
                     xlevels=list(B_avg.bright.female=c(23, 26, 29, 30, 32),
                                  proportion.close.50=seq(0,1,0.01))) %>% as_tibble 
# fill in plumage values
pred_table$plumage_label <- NA
pred_table$plumage_label[which(pred_table$B_avg.bright.female==23)] <- "23"
pred_table$plumage_label[which(pred_table$B_avg.bright.female==26)] <- "26"
pred_table$plumage_label[which(pred_table$B_avg.bright.female==29)] <- "29"
pred_table$plumage_label[which(pred_table$B_avg.bright.female==30)] <- "30"
pred_table$plumage_label[which(pred_table$B_avg.bright.female==32)] <- "32"

# plot
ggplot(pred_table, aes(x=proportion.close.50, y=fit, 
                       color=B_avg.bright.female,
                       linetype=plumage_label)) +
  geom_path(size=1.5) + xlim(0,0.3) + 
  scale_color_continuous(low="#330000",high="tan")

# Manuscript Figure S4
# plot with raw data
ggplot() +
  geom_path(data=pred_table, aes(x=proportion.close.50, y=fit, 
                       color=B_avg.bright.female,
                       linetype=plumage_label),
            size=1.5) +
  xlim(0,0.3) + 
  scale_color_continuous(low="#330000",high="tan") +
  geom_point(data=dat4, 
             aes(x=proportion.close.50, y=prop.ep, 
                 color=B_avg.bright.female, size=clutch.size),
             alpha=0.8) +
  ylab("Proportion EP young in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  ggtitle("Model predictions for female plumage color \ninteraction effect with points close") +
  theme_light()

ggsave("generated-files/FigS4_female plumage interation plot values 23, 26, 29, 30, 32_v5.png", h=5, w=6.5)
ggsave("generated-files/FigS4_female plumage interation plot values 23, 26, 29, 30, 32_v5.pdf", h=5, w=6.5)

# Manuscript Figure 2B
# plot with just 3 plumage lines (23, 26, 29)
pred_table_small <- subset(pred_table, pred_table$B_avg.bright.female==23 |
                          pred_table$B_avg.bright.female==26 |
                          pred_table$B_avg.bright.female==29)

ggplot() +
  geom_path(data=pred_table_small, aes(x=proportion.close.50, y=fit, 
                                 color=B_avg.bright.female,
                                 linetype=plumage_label),size=1.5) +
  xlim(0,0.3) + 
  scale_color_continuous(low="#330000",high="tan") +
  geom_point(data=dat4, 
             aes(x=proportion.close.50, y=prop.ep, 
                 color=B_avg.bright.female, size=clutch.size),
             alpha=0.8) +
  ylab("Proportion EP young in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  ggtitle("Model predictions for female plumage color \ninteraction effect with points close") +
  theme_light()

ggsave("generated-files/Fig2B_female plumage interation plot values 23, 26, 29_v5.png", h=5, w=6.5)
ggsave("generated-files/Fig2B_female plumage interation plot values 23, 26, 29_v5.pdf", h=5, w=6.5)




#-------------------------------------------------------------------------------
# Plot of interaction effect, female plumage with close50 on number EP mates 
# in replacement clutch

# make table for predicted values
pred_table_mates.rep <- effect(term="B_avg.bright.female:proportion.close.50", 
                     mod=mates.rep.cov.plumF.close50.int, x.var="proportion.close.50",
                     xlevels=list(B_avg.bright.female=c(23, 26, 29, 30),
                                  proportion.close.50=seq(0,0.3,0.01))) %>% as_tibble 
# fill in plumage values
pred_table_mates.rep$plumage_label <- NA
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==23)] <- "23"
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==26)] <- "26"
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==29)] <- "29"
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==30)] <- "30"


# plot with raw data
# for some reason some of the fit values are WAY too high. Need to add ylim to plot
ggplot() +
  geom_path(data=pred_table_mates.rep, aes(x=proportion.close.50, y=fit, 
                                 color=B_avg.bright.female,
                                 linetype=plumage_label), size=1.5) +
  xlim(0,0.3) + scale_color_continuous(low="#330000",high="tan") +
  geom_point(data=dat4, 
             aes(x=proportion.close.50, y=mates.replacement.ep, 
                 color=B_avg.bright.female, size=3),
             alpha=0.8) +
  ylab("Number EP mates in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  ggtitle("Model predictions for female plumage color \ninteraction effect with points close") +
  theme_light() + ylim(0,15)  

# Manuscript Figure S5
# plot with ylim=3
ggplot() +
  geom_path(data=pred_table_mates.rep, aes(x=proportion.close.50, y=fit, 
                                           color=B_avg.bright.female,
                                           linetype=plumage_label), size=1.5) +
  xlim(0,0.3) + scale_color_continuous(low="#330000",high="tan") +
  geom_point(data=dat4, 
             aes(x=proportion.close.50, y=mates.replacement.ep, 
                 color=B_avg.bright.female), size=4,
             alpha=0.8) +
  ylab("Number EP mates in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  ggtitle("Model predictions for female plumage color \ninteraction effect with points close on number of EP mates") +
  theme_light() + ylim(0,3)

ggsave("generated-files/FigS5_female plumage interation plot num EP values 23, 26, 29, 30_v5.png", h=5, w=6.5)
ggsave("generated-files/FigS5_female plumage interation plot num EP values 23, 26, 29, 30_v5.pdf", h=5, w=6.5)


