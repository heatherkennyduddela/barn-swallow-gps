
################################################################################
# Script to fit models of EP mating in replacement clutch for 2021
# GPS tagged females
# Heather Kenny-Duddela
################################################################################

# set working directory
setwd("C:/Users/heath/Documents/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/replacement-clutch-models")


# load data
dat4 <- read.csv("input-files/table for movement and mating models_BES_R.csv")

# check correlations between prop EPO, num EPO, and num EP sires

cor.test(dat4$prop.ep, dat4$num.ep, method="spearman")
# Spearman's rank correlation rho
# 
# # data:  dat4$prop.ep and dat4$num.ep
# # S = 28.28, p-value = 0.003053
# # alternative hypothesis: true rho is not equal to 0
# # sample estimates:
# #      rho 
# # 0.828609

cor.test(dat4$prop.ep, dat4$mates.replacement.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$prop.ep and dat4$mates.replacement.ep
# S = 21.784, p-value = 0.00113
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8679749 

cor.test(dat4$num.ep, dat4$mates.replacement.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$num.ep and dat4$mates.replacement.ep
# S = 36.045, p-value = 0.007582
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.781548


# Models!-----------------------------------------------------------------------

library(MASS) # for the glm.nb function
library(broom)# for calculating confidence intervals around coefficient estimates
library(AICcmodavg) # for calculating AICc
library(car) # for Anova function
library(dplyr) # for data table manipulation
library(ggplot2) # for plotting
library(jtools) # for effect plots

# helpful threads about offsets and weights in Poisson and binomial models
# https://stats.stackexchange.com/questions/264071/how-is-a-poisson-rate-regression-equal-to-a-poisson-regression-with-correspondin

# https://stats.stackexchange.com/questions/297859/can-weights-and-offset-lead-to-similar-results-in-poisson-regression 

# Help about glms from UCLA stats page
# https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faq-how-do-i-interpret-odds-ratios-in-logistic-regression/

# https://stats.oarc.ucla.edu/r/dae/logit-regression/


#---------------------------------------------------------------------------------
# Proportion of EP in replacement clutch predicted by different movement components

# Model proportion EP as binomial with a weights term for brood size
# for binomial family, response can be specified as numerical value 0 to 1 as 
# proportion of successful cases, with total number of cases given by weights


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


# prop past 50m
prop50m <- glm(prop.ep ~ proportion.past.50 + B_avg.bright.female, 
               weights=clutch.size, data=dat4, 
               family="binomial")
summary(prop50m)
hist(resid(prop50m, type="pearson"))
plot(fitted(prop50m) ~ resid(prop50m, type="pearson"))
confint(prop50m)


# prop past 100m
prop100m <- glm(prop.ep ~ proportion.past.100 + B_avg.bright.female, 
                weights=clutch.size, data=dat4,
                family="binomial")
summary(prop100m)
confint(prop100m)


# est max distance
estmax <- glm(prop.ep ~ scale(log.max.dist) + B_avg.bright.female, 
              weights=clutch.size, data=dat4,
              family="binomial")
summary(estmax)
confint(estmax)


####----------------------------------------------------------------------------
# Model number of EP young in replacement clutch as Poisson

# first check for overdispersion
hist(dat4$num.ep)
mean(dat4$num.ep, na.rm=T) # 1.8
var(dat4$num.ep, na.rm=T) # 2.17, only slightly overdispersed

# fit null poisson model
null.numEP.poisson <- glm(num.ep ~ 1, data=dat4, family="poisson")
AICc(null.numEP.poisson) # 37.21893

# fit null negative binomial model
null.numEP.nb <- glm.nb(num.ep ~ 1, data=dat4)
summary(null.numEP.nb)
AICc(null.numEP.nb) # 40.3671

### Poisson model performs better

# 90% KDE
numEP.kde90 <- glm(num.ep ~ scale(log.kde90) + 
                             B_avg.bright.female, 
                           data=dat4, family="poisson")
summary(numEP.kde90)
confint(numEP.kde90)
hist(resid(numEP.kde90, type="pearson"))
plot(fitted(numEP.kde90) ~ resid(numEP.kde90, type="pearson"))

# 50% KDE
numEP.kde50 <- glm(num.ep ~ scale(log.kde50) + 
                             B_avg.bright.female, 
                           data=dat4, family="poisson")
summary(numEP.kde50)
confint(numEP.kde50)
hist(resid(numEP.kde50, type="pearson"))
plot(fitted(numEP.kde50) ~ resid(numEP.kde50, type="pearson"))

# proportion points past 50
numEP.prop50m <- glm(num.ep ~ proportion.past.50 + 
                               B_avg.bright.female, 
                             data=dat4, family="poisson")
summary(numEP.prop50m)
confint(numEP.prop50m)

# proportion points past 100m
numEP.prop100m <- glm(num.ep ~ proportion.past.100 + 
                                B_avg.bright.female, 
                              data=dat4, family="poisson")
summary(numEP.prop100m)
confint(numEP.prop100m)

# est max distance
numEP.estmax <- glm(num.ep ~ scale(log.max.dist) + 
                              B_avg.bright.female, 
                            data=dat4, family="poisson")
summary(numEP.estmax)
confint(numEP.estmax)


####----------------------------------------------------------------------------
# Model number of EP sires in replacement clutch as Poisson

# first check for overdispersion
hist(dat4$mates.replacement.ep)
mean(dat4$mates.replacement.ep, na.rm=T) # 0.9
var(dat4$mates.replacement.ep, na.rm=T) # 0.544

null.rep.poisson <- glm(mates.replacement.ep ~ chicks.replacement , 
                        data=dat4, family="poisson")
AICc(null.rep.poisson) # 27.955

# quasi-poisson model
null.rep.quasi <- glm(mates.replacement.ep ~ chicks.replacement, data=dat4, family="quasipoisson")
summary(null.rep.quasi)  # (Dispersion parameter for quasipoisson family taken to be 0.6730484)

# negative binomial model
null.rep.nb <- glm.nb(mates.replacement.ep ~ chicks.replacement, data=dat4)
summary(null.rep.nb)
AICc(null.rep.nb) # 32.24147

### If anything the data is under dispersed, and the Poisson model performs better again. 
### proceed with fitting models using Poisson variance


# 90% KDE
mates.rep.kde90.cov <- glm(mates.replacement.ep ~ scale(log.kde90) + 
                             chicks.replacement +
                             B_avg.bright.female, 
                           data=dat4, family="poisson")
summary(mates.rep.kde90.cov)
confint(mates.rep.kde90.cov)
hist(resid(mates.rep.kde90.cov, type="pearson"))
plot(fitted(mates.rep.kde90.cov) ~ resid(mates.rep.kde90.cov, type="pearson"))

# 50% KDE
mates.rep.kde50.cov <- glm(mates.replacement.ep ~ scale(log.kde50) + 
                             chicks.replacement +
                             B_avg.bright.female, 
                           data=dat4, family="poisson")
summary(mates.rep.kde50.cov)
confint(mates.rep.kde50.cov)
hist(resid(mates.rep.kde50.cov, type="pearson"))
plot(fitted(mates.rep.kde50.cov) ~ resid(mates.rep.kde50.cov, type="pearson"))

# proportion points past 50
mates.rep.prop50m.cov <- glm(mates.replacement.ep ~ proportion.past.50 + 
                               chicks.replacement +
                               B_avg.bright.female, 
                             data=dat4, family="poisson")
summary(mates.rep.prop50m.cov)
confint(mates.rep.prop50m.cov)

# proportion points past 100m
mates.rep.prop100m.cov <- glm(mates.replacement.ep ~ proportion.past.100 + 
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


propEP.null <- glm(prop.ep ~ B_avg.bright.female, family = "binomial", 
                   data = dat4, weights = clutch.size)
summary(propEP.null)

# candidate set of models
propEP.models <- list(propEP.null, kde50, kde90, prop100m, prop50m, estmax)
# model names
propEP.modnames <- c("null","kde50","kde90","prop100m","prop50m","estmax")

propEP.table <- aictab(propEP.models, propEP.modnames)

# Manuscript Table S6
# top model is prop50m
write.csv(propEP.table, "output-files/TableS6_AIC table for propEP_BES.csv", row.names=F)


#-------------------------------------------------------------------------------
# AIC model selection - number EP young in replacement clutch with female plumage color


numEP.null <- glm(num.ep ~ B_avg.bright.female, family = "poisson", data = dat4)
summary(numEP.null)

# candidate set of models
numEP.models <- list(numEP.null, numEP.kde50, numEP.kde90, numEP.prop100m, 
                     numEP.prop50m, numEP.estmax)
# model names
numEP.modnames <- c("numEP.null","numEP.kde50","numEP.kde90","numEP.prop100m",
                    "numEP.prop50m","numEP.estmax")

numEP.table <- aictab(numEP.models, numEP.modnames)

# Manuscript Table S7
# top model is prop50m, though not better than null
write.csv(numEP.table, "output-files/TableS7_AIC table for numEP_BES.csv", row.names=F)


#------------------------------------------------------------------------------
### AIC model selection for number ep mates in replacement clutch

mates.rep.plum.cov <- glm(mates.replacement.ep ~ chicks.replacement +
                            B_avg.bright.female, 
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

# Manuscript Table S8
mates.rep.plum.table <- aictab(mates.rep.models.plum, mates.rep.modmanes.plum)

write.csv(mates.rep.plum.table, "output-files/TableS8_AIC table for mates in replacement_BES.csv", row.names=F)




################################################################################
#-------------------------------------------------------------------------------
# Compare the following models to test hypothesis that space use influences
# EPP alone, and after accounting for plumage color



# LIKELIHOOD RATIO TESTS FOR HYPOTHESIS TESTING

# 1) movement only
# 2) movement + female plumage



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

# 1) prop50
propEP.prop50m <- glm(prop.ep ~ proportion.past.50,
                            weights=clutch.size, 
                            data=dat4, family="binomial")


# 2) prop50m + female plumage
propEP.prop50m.plumF <- glm(prop.ep ~ proportion.past.50 + B_avg.bright.female,
                            weights=clutch.size, 
                            data=dat4, family="binomial")
summary(propEP.prop50m.plumF)



########### Likelihood ratio test for prop EP in replacement ------------------

### save table of results

# get table of coefficients
table.propEP.prop50m <- summary(propEP.prop50m)
# Anova results
anova.propEP.prop50m <- Anova(propEP.prop50m)
anova.propEP.prop50m2 <- rbind(rep(NA, 3), anova.propEP.prop50m)
# combine coefficients, and CIs, and Anova
table2.propEP.prop50m <- cbind(table.propEP.prop50m$coefficients, 
                               confint(propEP.prop50m), 
                               anova.propEP.prop50m2)

# get table of coefficients
table.propEP.prop50m.plumF <- summary(propEP.prop50m.plumF)
# anova result
anova.propEP.prop50m.plumF <- Anova(propEP.prop50m.plumF)
anova.propEP.prop50m.plumF2 <- rbind(rep(NA, 3), anova.propEP.prop50m.plumF)
# combine coefficients, and CIs and anova
table2.propEP.prop50m.plumF <- cbind(table.propEP.prop50m.plumF$coefficients, 
                               confint(propEP.prop50m.plumF),
                               anova.propEP.prop50m.plumF2)

table2.propEP.all <- rbind(table2.propEP.prop50m, rep(NA,9),
                           table2.propEP.prop50m.plumF)

write.csv(table2.propEP.all, "output-files/mod results propEP both.csv")

## To interpret the odds ratio in words, we can do the following: 
# 1) Take beta and multiple by 0.1 (for 10 percentage points)
# 2) Then exponentiate to get the odds ratio

odds.ratio <- exp(10.46*0.1) # 2.846

# lower CI
exp(2.69*0.1) # 1.308655
# upper CI
exp(20.67*0.1) # 7.901084



# Goodness of fit test----------------------------------------------------------

# For binomial model, we use the Hosmer-Lemeshow goodness of fit test
# https://thestatsgeek.com/2014/02/16/the-hosmer-lemeshow-goodness-of-fit-test-for-logistic-regression/

library(ResourceSelection)

# suggested that number of group for this test should be at least 1 larger than
# the number of parameters in the model. 

## close50 only

# three groups
h.3 <- hoslem.test(propEP.prop50m$y, fitted(propEP.prop50m), g=3)
h.3
# Hosmer and Lemeshow goodness of fit (GOF) test
# 
# data:  propEP.prop50m$y, fitted(propEP.prop50m)
# X-squared = 0.57052, df = 1, p-value = 0.4501


## close50.plumF

# four groups
h.4 <- hoslem.test(propEP.prop50m.plumF$y, 
                   fitted(propEP.prop50m.plumF), g=4)
h.4
# Hosmer and Lemeshow goodness of fit (GOF) test
# 
# data:  propEP.prop50m.plumF$y, fitted(propEP.prop50m.plumF)
# X-squared = 0.99548, df = 2, p-value = 0.6079

# High p-vals in both cases indicates no evidence of poor model fit, although
# the small sample size means the test has low power to detect misspecifications


################################################################################
#-------------------------------------------------------------------------------
# Models for number of EP young in replacement clutch

# close 50
numEP.rep.close50 <- glm(num.ep ~ proportion.past.50, family="poisson",
                               data=dat4)

# close 50 + female plumage
numEP.rep.close50.plumF <- glm(num.ep ~ B_avg.bright.female +
                                 proportion.past.50, family="poisson",
                               data=dat4)




### Likelihood ratio tests for number EP young --------------------------------

### save table of results

# get table of coefficients
table.numEP.rep.close50 <- summary(numEP.rep.close50)
# Anova results
anova.numEP.rep.close50 <- Anova(numEP.rep.close50)
anova.numEP.rep.close502 <- rbind(rep(NA, 3), anova.numEP.rep.close50)
# combine coefficients, and CIs, and Anova
table2.numEP.rep.close50 <- cbind(table.numEP.rep.close50$coefficients, 
                               confint(numEP.rep.close50), 
                               anova.numEP.rep.close502)

# get table of coefficients
table.numEP.rep.close50.plumF <- summary(numEP.rep.close50.plumF)
# anova result
anova.numEP.rep.close50.plumF <- Anova(numEP.rep.close50.plumF)
anova.numEP.rep.close50.plumF2 <- rbind(rep(NA, 3), anova.numEP.rep.close50.plumF)
# combine coefficients, and CIs and anova
table2.numEP.rep.close50.plumF <- cbind(
  table.numEP.rep.close50.plumF$coefficients, 
  confint(numEP.rep.close50.plumF),
  anova.numEP.rep.close50.plumF2)

table2.numEP.rep.all <- rbind(table2.numEP.rep.close50, rep(NA,9),
                           table2.numEP.rep.close50.plumF)

write.csv(table2.numEP.rep.all, "output-files/mod results numEP rep both.csv")

# Goodness of fit---------------------------------------------------------------

# For Poisson model, use deviance goodness of fit
# https://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/ 

## numEP close50
pchisq(numEP.rep.close50$deviance, df=numEP.rep.close50$df.residual, 
       lower.tail=FALSE)
# 0.139485, small p-val indicates that the model fit is bad. This p-val
# suggests our model fit is okay but not great


## numEP close50 plumF
pchisq(numEP.rep.close50.plumF$deviance, 
       df=numEP.rep.close50.plumF$df.residual, 
       lower.tail=FALSE)
# 0.09847568, this model fits slightly worse than the simpler model


################################################################################
#-------------------------------------------------------------------------------
# Models for number of EP mates in replacement clutch

# Best movement variable from AICc is proportion.close.50

# 1) close50
mates.rep.cov.close50 <- glm(mates.replacement.ep ~ chicks.replacement +
                                     proportion.past.50, 
                                   data=dat4, family="poisson")


# 2) close50 + female plumage
mates.rep.cov.close50.plumF <- glm(mates.replacement.ep ~ chicks.replacement +
                                     B_avg.bright.female +
                                     proportion.past.50, 
                                   data=dat4, family="poisson")


################### Likelihood ratio tests for number EP mates in replacement---

### save table of results

# get table of coefficients
table.mates.rep.cov.close50 <- summary(mates.rep.cov.close50)
# Anova results
anova.mates.rep.cov.close50 <- Anova(mates.rep.cov.close50)
anova.mates.rep.cov.close502 <- rbind(rep(NA, 3), anova.mates.rep.cov.close50)
# combine coefficients, and CIs, and Anova
table2.mates.rep.cov.close50 <- cbind(table.mates.rep.cov.close50$coefficients, 
                                  confint(mates.rep.cov.close50), 
                                  anova.mates.rep.cov.close502)

# get table of coefficients
table.mates.rep.cov.close50.plumF <- summary(mates.rep.cov.close50.plumF)
# anova result
anova.mates.rep.cov.close50.plumF <- Anova(mates.rep.cov.close50.plumF)
anova.mates.rep.cov.close50.plumF2 <- rbind(rep(NA, 3), anova.mates.rep.cov.close50.plumF)
# combine coefficients, and CIs and anova
table2.mates.rep.cov.close50.plumF <- cbind(
  table.mates.rep.cov.close50.plumF$coefficients, 
  confint(mates.rep.cov.close50.plumF),
  anova.mates.rep.cov.close50.plumF2)

table2.mates.rep.cov.all <- rbind(table2.mates.rep.cov.close50, rep(NA,9),
                              table2.mates.rep.cov.close50.plumF)

write.csv(table2.mates.rep.cov.all, "output-files/mod results mates rep both.csv")

# Goodness of fit --------------------------------------------------------------

## mates close50
pchisq(mates.rep.cov.close50$deviance, 
       df=mates.rep.cov.close50$df.residual, 
       lower.tail=FALSE)
# 0.5399633, higher p-val supports our null that the current model fit is good

## mates close50 plumF
pchisq(mates.rep.cov.close50.plumF$deviance, 
       df=mates.rep.cov.close50.plumF$df.residual, 
       lower.tail=FALSE)
# 0.4335756, again here the model fit is pretty good

################################################################################
#-------------------------------------------------------------------------------

## Plot replacement prop EPO and replacement num EPO effects

# Proportion EPO

# make new data to predict from
pred.propEPO.y <- data.frame(
  proportion.past.50=seq(0,1,0.001),
  B_avg.bright.female=rep(mean(dat4$B_avg.bright.female), 1001))

# add predictions
pred.propEPO.y$fit <- predict.glm(propEP.prop50m.plumF, 
                        newdata=pred.propEPO.y[,1:2],
                        type="response")
# add the fitted values and SE on the scale of the model (not response)
link <- predict.glm(propEP.prop50m.plumF, 
              newdata=pred.propEPO.y[,1:2],
              type="link", se.fit=T)
pred.propEPO.y$fit_link <- link$fit
pred.propEPO.y$se_link <- link$se.fit
# get inverse-link function
fam <- family(propEP.prop50m.plumF)
ilink <- fam$linkinv
# calculate CIs as back-transformed fit_link +/- 2*se_link
pred.propEPO.y$upper <- ilink(pred.propEPO.y$fit_link + 
                                2*pred.propEPO.y$se_link)
pred.propEPO.y$lower <- ilink(pred.propEPO.y$fit_link - 
                                2*pred.propEPO.y$se_link)

# plot fit line and confidence intervals
ggplot() +
  geom_line(data=pred.propEPO.y, 
            aes(x=proportion.past.50, y=fit), size=2, 
            color="darkblue") +
  geom_ribbon(data=pred.propEPO.y, 
              aes(ymin=lower, ymax=upper, 
                  x=proportion.past.50), alpha=0.2, fill="darkblue") +
  xlim(0.7,1) + ylim(0,1) +
  geom_point(data=dat4, aes(x=proportion.past.50, y=prop.ep,
                            size=clutch.size),alpha=0.5) +
  xlab("Proportion of GPS points farther than 50m from the barn") +
  ylab("Proportion of EP offspring in replacement clutch") +
  ggtitle("The effect of movement on the proportion of\nEPO in the replacement clutch")

# Figure 4a
ggsave("output-files/effect of movement on prop EPO.png", h=4, w=5)

## plots using jtools
library(jtools)

effect_plot(propEP.prop50m.plumF, pred=proportion.past.50, interval=T) +
  ylim(0,1)


# Number of EPO ----------------------------------------------------------------

# make new data to predict from
pred.numEPO.y <- data.frame(
  proportion.past.50=seq(0,1,0.001),
  B_avg.bright.female=rep(mean(dat4$B_avg.bright.female), 1001))

# add predictions
pred.numEPO.y$fit <- predict.glm(numEP.rep.close50.plumF, 
                                  newdata=pred.numEPO.y[,1:2],
                                  type="response")
# add the fitted values and SE on the scale of the model (not response)
link.numEPO <- predict.glm(numEP.rep.close50.plumF, 
                    newdata=pred.numEPO.y[,1:2],
                    type="link", se.fit=T)
pred.numEPO.y$fit_link <- link.numEPO$fit
pred.numEPO.y$se_link <- link.numEPO$se.fit
# get inverse-link function
fam.numEPO <- family(numEP.rep.close50.plumF)
ilink.numEPO <- fam.numEPO$linkinv
# calculate CIs as back-transformed fit_link +/- 2*se_link
pred.numEPO.y$upper <- ilink.numEPO(pred.numEPO.y$fit_link + 
                                2*pred.numEPO.y$se_link)
pred.numEPO.y$lower <- ilink.numEPO(pred.numEPO.y$fit_link - 
                                2*pred.numEPO.y$se_link)

# plot fit line and confidence intervals
ggplot() +
  geom_line(data=pred.numEPO.y, 
            aes(x=proportion.past.50, y=fit), 
            size=2, color="darkblue", linetype="dashed") +
  geom_ribbon(data=pred.numEPO.y, 
              aes(ymin=lower, ymax=upper, 
                  x=proportion.past.50), alpha=0.2, fill="darkblue") +
  xlim(0.7,1) + ylim(0,6) +
  geom_point(data=dat4, aes(x=proportion.past.50, y=num.ep), 
             size=4, alpha=0.5) +
  xlab("Proportion of GPS points farther than 50m from the barn") +
  ylab("Number of EP offspring in replacement clutch") +
  ggtitle("The effect of movement on the number of\nEPO in the replacement clutch")

# Figure 4b  
ggsave("output-files/effect of movement on num EPO.png", h=4, w=5)


effect_plot(numEP.rep.close50.plumF, pred=proportion.past.50, interval=T)



