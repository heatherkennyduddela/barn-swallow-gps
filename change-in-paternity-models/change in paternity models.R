
################################################################################
# Script for fitting models for change in EP mating for 2021 GPS females
# Heather Kenny-Duddela
# June 14, 2024
################################################################################

# load libraries
library(broom)# for calculating confidence intervals around coefficient estimates
library(car) # for Anova function

# set working directory
setwd("C:/Users/heath/Documents/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/change-in-paternity-models")

# load data
dat4 <- read.csv("input-files/table for movement and mating models_BES_R.csv")

# check correlations between propEP, number EPO, and number EP sires

cor.test(dat4$diff.ep.prop, dat4$diff.num.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$diff.ep.prop and dat4$diff.num.ep
# S = 13.076, p-value = 0.0001567
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9207488 

cor.test(dat4$diff.ep.prop, dat4$diff.ep.sires , method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$diff.ep.prop and dat4$diff.ep.sires
# S = 17.838, p-value = 0.0005236
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.8918885 

cor.test(dat4$diff.num.ep, dat4$diff.ep.sires , method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$diff.num.ep and dat4$diff.ep.sires
# S = 10.454, p-value = 6.527e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9366433

#-------------------------------------------------------------------------------
# Models for change in paternity - proportion EP offspring

# USE RAW PLUMAGE VALUES
# Use proportion.close.50 as movement variable

# 1) prop50m
diff.propEP.close50 <- lm(diff.ep.prop ~  proportion.past.50, 
                                data=dat4)
hist(diff.propEP.close50$residuals)
summary(diff.propEP.close50)
# Multiple R-squared:  0.2234,	Adjusted R-squared:  0.1263


# 2) prop50m + female plumage 
diff.propEP.close50.plumF <- lm(diff.ep.prop ~  proportion.past.50 + 
                                        B_avg.bright.female, 
                                      data=dat4)
hist(diff.propEP.close50.plumF$residuals)
summary(diff.propEP.close50.plumF)
# Multiple R-squared:  0.3886,	Adjusted R-squared:  0.214

# Anova tests for significance -------------------------------------------------

### save table of results

# get table of coefficients
table.diff.propEP.close50 <- summary(diff.propEP.close50)
# Anova results
anova.diff.propEP.close50 <- Anova(diff.propEP.close50)
anova.diff.propEP.close502 <- rbind(rep(NA, 4), anova.diff.propEP.close50)
# combine coefficients, and CIs, and Anova
table2.diff.propEP.close50 <- cbind(table.diff.propEP.close50$coefficients, 
                                      confint(diff.propEP.close50), 
                                      anova.diff.propEP.close502[1:2, ])

# get table of coefficients
table.diff.propEP.close50.plumF <- summary(diff.propEP.close50.plumF)
# anova result
anova.diff.propEP.close50.plumF <- Anova(diff.propEP.close50.plumF)
anova.diff.propEP.close50.plumF2 <- rbind(rep(NA, 4), anova.diff.propEP.close50.plumF)
# combine coefficients, and CIs and anova
table2.diff.propEP.close50.plumF <- cbind(
  table.diff.propEP.close50.plumF$coefficients, 
  confint(diff.propEP.close50.plumF),
  anova.diff.propEP.close50.plumF2[1:3, ])

table2.diff.propEP.all <- rbind(table2.diff.propEP.close50, rep(NA,10),
                                  table2.diff.propEP.close50.plumF)

write.csv(table2.diff.propEP.all, 
          "output-files/mod results diff propEP both.csv")




################################################################################
#-------------------------------------------------------------------------------
# Models for change in number of EP sires
# movement component is proportion.close.50

# check correlations to decide if we should include clutch size as a covariate

# correlation between number of sires and clutch size for collected
ggplot(dat4, aes(x=chicks.collected, y=mates.collected.ep)) + 
  geom_point(position=position_jitter(h=0.01, w=0.01), 
             shape=1) +
  geom_smooth(method=lm, se=F)

cor.test(dat4$chicks.collected, dat4$mates.collected.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$chicks.collected and dat4$mates.collected.ep
# S = 318.68, p-value = 0.1664
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.4485426 

# correlation between number of sires and clutch size for replacement
ggplot(dat4, aes(x=chicks.replacement, y=mates.replacement.ep)) + 
  geom_point(position=position_jitter(h=0.01, w=0.01)) +
  geom_smooth(method=lm, se=F)

cor.test(dat4$chicks.replacement, dat4$mates.replacement.ep, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$chicks.replacement and dat4$mates.replacement.ep
# S = 195.17, p-value = 0.6131
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# -0.1828674 

# correlation between total number of chicks and change in number of sires
# make column for total number of chicks
dat4$tot.chicks <- dat4$chicks.collected + dat4$chicks.replacement

ggplot(dat4, aes(x=tot.chicks, y=diff.ep.sires)) + 
  geom_point(position=position_jitter(h=0, w=0.05)) +
  geom_smooth(method=lm, se=F) +
  xlab("Total chicks across both clutches") +
  ylab("Change in number of EP sires") +
  ggtitle("Relationship between total young sampled and change in EP sires")

cor.test(dat4$tot.chicks, dat4$diff.ep.sires, method="spearman")
# Spearman's rank correlation rho
# 
# data:  dat4$tot.chicks and dat4$diff.ep.sires
# S = 161.65, p-value = 0.9557
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.02028185 

# Conclusion: no need to include a covariate for number of chicks

#-------------------------------------------------------------------------------
# Fit the models

# 1) close50
diff.sires.close50 <- lm(diff.ep.sires ~ proportion.past.50, data = dat4)
hist(diff.sires.close50$residuals)
summary(diff.sires.close50)
# Multiple R-squared:   0.17,	Adjusted R-squared:  0.06621

# 2) close50 + female plumage
diff.sires.close50.plumF <- lm(diff.ep.sires ~ proportion.past.50 + 
                             B_avg.bright.female, data = dat4)
hist(diff.sires.close50.plumF$residuals)
summary(diff.sires.close50.plumF)
# Multiple R-squared:  0.3782,	Adjusted R-squared:  0.2006 

################### Likelihood ratio tests change in number of EP sires

### save table of results

# get table of coefficients
table.diff.sires.close50 <- summary(diff.sires.close50)
# Anova results
anova.diff.sires.close50 <- Anova(diff.sires.close50)
anova.diff.sires.close502 <- rbind(rep(NA, 4), anova.diff.sires.close50)
# combine coefficients, and CIs, and Anova
table2.diff.sires.close50 <- cbind(table.diff.sires.close50$coefficients, 
                                    confint(diff.sires.close50), 
                                    anova.diff.sires.close502[1:2, ])

# get table of coefficients
table.diff.sires.close50.plumF <- summary(diff.sires.close50.plumF)
# anova result
anova.diff.sires.close50.plumF <- Anova(diff.sires.close50.plumF)
anova.diff.sires.close50.plumF2 <- rbind(rep(NA, 4), anova.diff.sires.close50.plumF)
# combine coefficients, and CIs and anova
table2.diff.sires.close50.plumF <- cbind(
  table.diff.sires.close50.plumF$coefficients, 
  confint(diff.sires.close50.plumF),
  anova.diff.sires.close50.plumF2[1:3, ])

table2.diff.sires.all <- rbind(table2.diff.sires.close50, rep(NA,10),
                                table2.diff.sires.close50.plumF)

write.csv(table2.diff.sires.all, 
          "output-files/mod results diff sires both.csv")





################################################################################
#-------------------------------------------------------------------------------
# Models for change in number of EP young


# 1) close50
diff.numEP.close50 <- lm(diff.num.ep ~ proportion.past.50, data = dat4)
hist(diff.numEP.close50$residuals)
summary(diff.numEP.close50)
# Multiple R-squared:  0.1406,	Adjusted R-squared:  0.03315 


# 2) close50 + female plumage
diff.numEP.close50.plumF <- lm(diff.num.ep ~ proportion.past.50 + 
                                 B_avg.bright.female, data = dat4)
hist(diff.numEP.close50.plumF$residuals)
summary(diff.numEP.close50.plumF)
#  R-squared:  0.4024,	Adjusted R-squared:  0.2317 



#-------------------------------------------------------------------------------
# Likelihood ratio tests for change in number EP young

### save table of results

# get table of coefficients
table.diff.numEP.close50 <- summary(diff.numEP.close50)
# Anova results
anova.diff.numEP.close50 <- Anova(diff.numEP.close50)
anova.diff.numEP.close502 <- rbind(rep(NA, 4), anova.diff.numEP.close50)
# combine coefficients, and CIs, and Anova
table2.diff.numEP.close50 <- cbind(table.diff.numEP.close50$coefficients, 
                                   confint(diff.numEP.close50), 
                                   anova.diff.numEP.close502[1:2, ])

# get table of coefficients
table.diff.numEP.close50.plumF <- summary(diff.numEP.close50.plumF)
# anova result
anova.diff.numEP.close50.plumF <- Anova(diff.numEP.close50.plumF)
anova.diff.numEP.close50.plumF2 <- rbind(rep(NA, 4), anova.diff.numEP.close50.plumF)
# combine coefficients, and CIs and anova
table2.diff.numEP.close50.plumF <- cbind(
  table.diff.numEP.close50.plumF$coefficients, 
  confint(diff.numEP.close50.plumF),
  anova.diff.numEP.close50.plumF2[1:3, ])

table2.diff.numEP.all <- rbind(table2.diff.numEP.close50, rep(NA,10),
                               table2.diff.numEP.close50.plumF)

write.csv(table2.diff.numEP.all, 
          "output-files/mod results diff numEP both.csv")



