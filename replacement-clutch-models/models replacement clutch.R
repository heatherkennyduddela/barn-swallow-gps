
################################################################################
# Script to fit models of EP mating in replacement clutch for 2021
# GPS tagged females
# Heather Kenny-Duddela
################################################################################

# set working directory
setwd("C:/Users/heath/Documents/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/replacement-clutch-models")


# load data
dat4 <- read.csv("input-files/table for movement and mating models_BES.csv")



# Models!-----------------------------------------------------------------------

library(MASS) # for the glm.nb function
library(broom)# for calculating confidence intervals around coefficient estimates
library(AICcmodavg) # for calculating AICc
library(car) # for Anova function
library(dplyr) # for data table manipulation
library(ggplot2) # for plotting 

# helpful threads about offsets and weights in Poisson and binomial models
# https://stats.stackexchange.com/questions/264071/how-is-a-poisson-rate-regression-equal-to-a-poisson-regression-with-correspondin

# https://stats.stackexchange.com/questions/297859/can-weights-and-offset-lead-to-similar-results-in-poisson-regression 

# Help about glms from UCLA stats page
# https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faq-how-do-i-interpret-odds-ratios-in-logistic-regression/

# https://stats.oarc.ucla.edu/r/dae/logit-regression/


#---------------------------------------------------------------------------------
# Proportion of EP in replacement clutch predicted by different movement components

# Model proportion EP as binomial with a weights term for brood size
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

# proportion points within 50
numEP.prop50m <- glm(num.ep ~ proportion.close.50 + 
                               B_avg.bright.female, 
                             data=dat4, family="poisson")
summary(numEP.prop50m)
confint(numEP.prop50m)

# proportion points within 100m
numEP.prop100m <- glm(num.ep ~ proportion.close.100 + 
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

# quasi-poisson can model
null.rep.quasi <- glm(mates.replacement.ep ~ chicks.replacement, data=dat4, family="quasipoisson")
summary(null.rep.quasi)  # (Dispersion parameter for quasipoisson family taken to be 0.6730484)

# negative binomial model
null.rep.nb <- glm.nb(mates.replacement.ep ~ chicks.replacement, data=dat4)
summary(null.rep.nb)
AICc(null.rep.nb) # 32.24147

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
# Compare the following models to test hypothesis that space use interacts with 
# female plumage



# LIKELIHOOD RATIO TESTS FOR HYPOTHESIS TESTING

# 2) movement + female plumage
# 3) movement * female plumage


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


# 1) prop50m + female plumage
propEP.prop50m.plumF <- glm(prop.ep ~ proportion.close.50 + B_avg.bright.female,
                            weights=clutch.size, 
                            data=dat4, family="binomial")
summary(propEP.prop50m.plumF)

# 2) prop50m * female plumage
propEP.prop50m.plumF.int <- glm(prop.ep ~ proportion.close.50 * B_avg.bright.female, 
                                weights=clutch.size, 
                                data=dat4, family="binomial")
summary(propEP.prop50m.plumF.int)



########### Likelihood ratio test for prop EP in replacement ------------------


# 1) movement + female plumage
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


# 2) movement * female plumage
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



################################################################################
#-------------------------------------------------------------------------------
# Models for number of EP young in replacement clutch

# close 50 + female plumage
numEP.rep.plumF.close50 <- glm(num.ep ~ B_avg.bright.female +
                                 proportion.close.50, family="poisson",
                               data=dat4)

# close 50 * female plumage
numEP.rep.plumF.close50.int <- glm(num.ep ~ B_avg.bright.female *
                                 proportion.close.50, family="poisson",
                               data=dat4)



### Likelihood ratio tests for number EP young --------------------------------

# movement + female plumage
Anova(numEP.rep.plumF.close50)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep
# LR Chisq Df Pr(>Chisq)  
# B_avg.bright.female   0.2078  1    0.64849  
# proportion.close.50   3.1863  1    0.07426 .
confint(numEP.rep.plumF.close50)
# 2.5 %    97.5 %
#   (Intercept)          -3.2806612 3.5777293
# B_avg.bright.female  -0.1013904 0.1642695
# proportion.close.50 -12.8628397 0.4681680
summary(numEP.rep.plumF.close50)
# glm(formula = num.ep ~ B_avg.bright.female + proportion.close.50, 
#     family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.0052  -0.9728  -0.0644   0.7819   1.3117  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          0.22513    1.73198   0.130    0.897
# B_avg.bright.female  0.03045    0.06682   0.456    0.649
# proportion.close.50 -5.32598    3.32567  -1.601    0.109
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 15.250  on 9  degrees of freedom
# Residual deviance: 12.064  on 7  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 37.533
# 
# Number of Fisher Scoring iterations: 6


# movement * female plumage
Anova(numEP.rep.plumF.close50.int)
# Analysis of Deviance Table (Type II tests)
# 
# Response: num.ep
# LR Chisq Df Pr(>Chisq)  
# B_avg.bright.female                       0.2078  1    0.64849  
# proportion.close.50                       3.1863  1    0.07426 .
# B_avg.bright.female:proportion.close.50   0.1158  1    0.73364  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
confint(numEP.rep.plumF.close50.int)
# 2.5 %     97.5 %
#   (Intercept)                              -5.5358456  7.3590313
# B_avg.bright.female                      -0.2553697  0.2518696
# proportion.close.50                     -70.9357660 42.9747556
# B_avg.bright.female:proportion.close.50  -1.8445344  2.4704020
summary(numEP.rep.plumF.close50.int)
# glm(formula = num.ep ~ B_avg.bright.female * proportion.close.50, 
#     family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.9491  -0.9510  -0.1581   0.8555   1.2176  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                               1.174039   3.255869   0.361    0.718
# B_avg.bright.female                      -0.006966   0.128079  -0.054    0.957
# proportion.close.50                     -15.068732  28.720210  -0.525    0.600
# B_avg.bright.female:proportion.close.50   0.371260   1.085608   0.342    0.732
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 15.250  on 9  degrees of freedom
# Residual deviance: 11.948  on 6  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 39.417
# 
# Number of Fisher Scoring iterations: 6


################################################################################
#-------------------------------------------------------------------------------
# Models for number of EP mates in replacement clutch

# Best movement variable from AICc is proportion.close.50

# 1) close50 + female plumage
mates.rep.cov.plumF.close50 <- glm(mates.replacement.ep ~ chicks.replacement +
                                     B_avg.bright.female +
                                     proportion.close.50, 
                                   data=dat4, family="poisson")

# 2) close50 * female plumage
mates.rep.cov.plumF.close50.int <- glm(mates.replacement.ep ~ chicks.replacement +
                                         B_avg.bright.female * proportion.close.50, 
                                       data=dat4, family="poisson")




################### Likelihood ratio tests for number EP mates in replacement---


# 1) movement + female plumage
Anova(mates.rep.cov.plumF.close50)
# Analysis of Deviance Table (Type II tests)
# 
# Response: mates.replacement.ep
# LR Chisq Df Pr(>Chisq)
# chicks.replacement   0.10741  1     0.7431
# B_avg.bright.female  0.09018  1     0.7640
# proportion.close.50  1.10586  1     0.2930
confint(mates.rep.cov.plumF.close50)
# 2.5 %    97.5 %
#   (Intercept)          -5.4813916 4.6938295
# chicks.replacement   -0.7490114 0.5401911
# B_avg.bright.female  -0.1606925 0.2162301
# proportion.close.50 -15.8251888 3.5957166
summary(mates.rep.cov.plumF.close50)
# glm(formula = mates.replacement.ep ~ chicks.replacement + B_avg.bright.female + 
#       proportion.close.50, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.40717  -0.69851   0.00071   0.40443   0.86470  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         -0.04594    2.52415  -0.018    0.985
# chicks.replacement  -0.10483    0.31960  -0.328    0.743
# B_avg.bright.female  0.02787    0.09261   0.301    0.763
# proportion.close.50 -4.61644    4.73662  -0.975    0.330
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.4417  on 9  degrees of freedom
# Residual deviance: 5.9080  on 6  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 29.135
# 
# Number of Fisher Scoring iterations: 5


# 2) movement * female plumage
Anova(mates.rep.cov.plumF.close50.int, type="III")
# Analysis of Deviance Table (Type III tests)
# 
# Response: mates.replacement.ep
# LR Chisq Df Pr(>Chisq)  
# chicks.replacement                        2.5015  1    0.11374  
# B_avg.bright.female                       3.3595  1    0.06682 .
# proportion.close.50                       4.1426  1    0.04182 *
#   B_avg.bright.female:proportion.close.50   3.9334  1    0.04734 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
confint(mates.rep.cov.plumF.close50.int)
# 2.5 %      97.5 %
#   (Intercept)                               -0.86597571 35.48321299
# chicks.replacement                        -0.25242044  2.66429645
# B_avg.bright.female                       -1.78607036  0.04462712
# proportion.close.50                     -461.89904986 -6.34864949
# B_avg.bright.female:proportion.close.50    0.07398172 16.94080603
summary(mates.rep.cov.plumF.close50.int)
# glm(formula = mates.replacement.ep ~ chicks.replacement + B_avg.bright.female * 
#       proportion.close.50, family = "poisson", data = dat4)
# 
# Deviance Residuals: 
#   1         2         3         4         5         6         7         8         9        10  
# -0.02555  -0.21620  -0.14584   0.01191   0.44837   0.44604  -0.69392  -0.27251  -0.86167   0.45560  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)                               14.5125     8.8819   1.634   0.1023  
# chicks.replacement                         1.0762     0.7125   1.510   0.1309  
# B_avg.bright.female                       -0.7202     0.4470  -1.611   0.1072  
# proportion.close.50                     -195.2665   112.0351  -1.743   0.0814 .
# B_avg.bright.female:proportion.close.50    7.0629     4.1353   1.708   0.0876 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 7.4417  on 9  degrees of freedom
# Residual deviance: 1.9746  on 5  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 27.202
# 
# Number of Fisher Scoring iterations: 5



################################################################################
#-------------------------------------------------------------------------------



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


# Manuscript Figure S6
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
  scale_size(breaks = c(2,3,4,5),range = c(2,6)) +
  ylab("Proportion EP offspring in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  theme_light()

ggsave("output-files/FigS6_female plumage interation plot values 23, 26, 29, 30, 32_v5.png", h=2.7, w=3.3, scale=1.9, units="in")
ggsave("output-files/FigS6_female plumage interation plot values 23, 26, 29, 30, 32_v5.pdf", h=2.7, w=3.3, scale=1.9)

# Manuscript Figure 2A
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
  scale_size(breaks = c(2,3,4,5),range = c(2,6)) +
  ylab("Proportion EP offspring in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  theme_light()

ggsave("output-files/Fig2A_female plumage interation plot values 23, 26, 29_v5.png", h=2.5, w=3.3, scale=1.9)
ggsave("output-files/Fig2A_female plumage interation plot values 23, 26, 29_v5.pdf", h=2.5, w=3.3, scale=1.9)




#-------------------------------------------------------------------------------
# Plot of interaction effect, female plumage with close50 on number EP sires 
# in replacement clutch

# make table for predicted values
pred_table_mates.rep <- effect(term="B_avg.bright.female:proportion.close.50", 
                               mod=mates.rep.cov.plumF.close50.int, 
                               x.var="proportion.close.50",
                               xlevels=list(B_avg.bright.female=c(23, 26, 29, 30, 32),
                                            proportion.close.50=seq(0,0.3,0.01))) %>% as_tibble 
# fill in plumage values
pred_table_mates.rep$plumage_label <- NA
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==23)] <- "23"
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==26)] <- "26"
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==29)] <- "29"
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==30)] <- "30"
pred_table_mates.rep$plumage_label[which(pred_table_mates.rep$B_avg.bright.female==32)] <- "32"

# plot with raw data
# for some reason some of the fit values are WAY too high. Need to add ylim to plot
ggplot() +
  geom_path(data=pred_table_mates.rep, aes(x=proportion.close.50, y=fit, 
                                           color=B_avg.bright.female,
                                           linetype=plumage_label), size=1.5) +
  xlim(0,0.3) + scale_color_continuous(low="#330000",high="tan") +
  geom_point(data=dat4, 
             aes(x=proportion.close.50, y=mates.replacement.ep, 
                 color=B_avg.bright.female), size=3,
             alpha=0.8) +
  ylab("Number EP sires in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  ggtitle("Model predictions for number of EP sires\nin the replacement clutch") +
  theme_light() + ylim(0,15)  


# Manuscript Figure S7
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
  ylab("Number EP sires in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  theme_light() + ylim(0,3)

ggsave("output-files/FigS7_female plumage interation plot num EP values 23, 26, 29, 30, 32_v5.png", h=2.5, w=3.3, scale=1.9)
ggsave("output-files/FigS7_female plumage interation plot num EP values 23, 26, 29, 30, 32_v5.pdf", h=2.5, w=3.3, scale=1.9)


# Figure 2B in main text
# plot with just 3 plumage lines (23, 26, 29)
pred_table_mates.rep_small <- subset(pred_table_mates.rep, 
                             pred_table_mates.rep$B_avg.bright.female==23 |
                             pred_table_mates.rep$B_avg.bright.female==26 |
                             pred_table_mates.rep$B_avg.bright.female==29)

ggplot() +
  geom_path(data=pred_table_mates.rep_small, aes(x=proportion.close.50, y=fit, 
                                           color=B_avg.bright.female,
                                           linetype=plumage_label), size=1.5) +
  xlim(0,0.3) + scale_color_continuous(low="#330000",high="tan") +
  geom_point(data=dat4, 
             aes(x=proportion.close.50, y=mates.replacement.ep, 
                 color=B_avg.bright.female), size=4,
             alpha=0.8) +
  ylab("Number EP sires in replacement clutch") +
  xlab("Proportion GPS points within 50m of barn") +
  theme_light() + ylim(0,3)

ggsave("output-files/Fig2B_female plumage interation plot num EP values 23, 26, 29_v5.png", h=2.5, w=3.3, scale=1.9)
ggsave("output-files/Fig2B_female plumage interation plot num EP values 23, 26, 29_v5.pdf", h=2.5, w=3.3, scale=1.9)

