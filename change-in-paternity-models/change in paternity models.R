
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
dat4 <- read.csv("input-files/table for movement and mating models_BES.csv")



#-------------------------------------------------------------------------------
# Models for change in paternity - proportion EP offspring

# USE RAW PLUMAGE VALUES
# Use proportion.close.50 as movement variable


# 1) prop50m + female plumage - female plumage sig but not movement
diff.propEP.prop50m.plumF <- lm(diff.ep.prop ~  proportion.close.50 + 
                                        B_avg.bright.female, 
                                      data=dat4)
hist(diff.propEP.prop50m.plumF$residuals)


# 2) prop50m * female plumage - movement and interaction significant
diff.propEP.prop50m.plumF.int <- lm(diff.ep.prop ~ proportion.close.50 * 
                                            B_avg.bright.female, 
                                          data=dat4)
hist(diff.propEP.prop50m.plumF.int$residuals)




# Anova tests for significance -------------------------------------------------


# 1) movement + female plumage
Anova(diff.propEP.prop50m.plumF)
# Anova Table (Type II tests)
# 
# Response: diff.ep.prop
# Sum Sq Df F value Pr(>F)
# proportion.close.50 1.27539  1  3.3881 0.1082
# B_avg.bright.female 0.71223  1  1.8921 0.2114
# Residuals           2.63501  7 
confint(diff.propEP.prop50m.plumF)
# 2.5 %    97.5 %
#   (Intercept)         -4.48922318 1.8965115
# proportion.close.50 -8.66340929 1.0793773
# B_avg.bright.female -0.05069497 0.1916961
summary(diff.propEP.prop50m.plumF)
# lm(formula = diff.ep.prop ~ proportion.close.50 + B_avg.bright.female, 
#    data = dat4)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.75641 -0.45427 -0.09037  0.45609  0.87017 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)         -1.29636    1.35026  -0.960    0.369
# proportion.close.50 -3.79202    2.06011  -1.841    0.108
# B_avg.bright.female  0.07050    0.05125   1.376    0.211
# 
# Residual standard error: 0.6135 on 7 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.3886,	Adjusted R-squared:  0.214 
# F-statistic: 2.225 on 2 and 7 DF,  p-value: 0.1787



# 2) movement * female plumage
Anova(diff.propEP.prop50m.plumF.int, type="III")
# Anova Table (Type III tests)
# 
# Response: diff.ep.prop
# Sum Sq Df F value Pr(>F)
# (Intercept)                             0.14601  1  0.4224 0.5398
# proportion.close.50                     0.72078  1  2.0850 0.1989
# B_avg.bright.female                     0.07698  1  0.2227 0.6537
# proportion.close.50:B_avg.bright.female 0.56079  1  1.6222 0.2499
# Residuals                               2.07422  6   
confint(diff.propEP.prop50m.plumF.int)
# 2.5 %     97.5 %
#   (Intercept)                              -4.9259945  8.4889465
# proportion.close.50                     -84.0017740 21.6535687
# B_avg.bright.female                      -0.3120622  0.2111612
# proportion.close.50:B_avg.bright.female  -0.9688231  3.0722276
summary(diff.propEP.prop50m.plumF.int)
# lm(formula = diff.ep.prop ~ proportion.close.50 * B_avg.bright.female, 
#    data = dat4)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.7755 -0.4400  0.1246  0.4091  0.5814 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                               1.78148    2.74120   0.650    0.540
# proportion.close.50                     -31.17410   21.58953  -1.444    0.199
# B_avg.bright.female                      -0.05045    0.10692  -0.472    0.654
# proportion.close.50:B_avg.bright.female   1.05170    0.82575   1.274    0.250
# 
# Residual standard error: 0.588 on 6 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.5187,	Adjusted R-squared:  0.2781 
# F-statistic: 2.156 on 3 and 6 DF,  p-value: 0.1944



################################################################################
#-------------------------------------------------------------------------------
# Models for change in number of EP sires
# movement component is proportion.close.50

# check correlations to decide if we should include clutch size as a covariate

# correlation between number of sires and clutch size for collected
ggplot(dat4, aes(x=chicks.collected, y=mates.collected.ep)) + 
  geom_point(position=position_jitter(h=0.1, w=0.1)) +
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
  geom_point(position=position_jitter(h=0.1, w=0.1)) +
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



# 1) close50 + female plumage
diff.sires.close50.plumF <- lm(diff.ep.sires ~ proportion.close.50 + 
                             B_avg.bright.female, data = dat4)
hist(diff.sires.close50.plumF$residuals)




# 2) close50 * female plumage
diff.sires.close50.plumF.int <- lm(diff.ep.sires ~ proportion.close.50 * 
                                 B_avg.bright.female, data = dat4)
hist(diff.sires.close50.plumF.int$residuals)





################### Likelihood ratio tests change in number of EP sires


# 1) movement + female plumage
Anova(diff.sires.close50.plumF)
# Anova Table (Type II tests)
# 
# Response: diff.ep.sires
# Sum Sq Df F value Pr(>F)
# proportion.close.50 2.0804  1  2.7556 0.1409
# B_avg.bright.female 1.7705  1  2.3450 0.1695
# Residuals           5.2849  7  
confint(diff.sires.close50.plumF)
# 2.5 %    97.5 %
#   (Intercept)          -6.47269655 2.5708177
# proportion.close.50 -11.74201240 2.0557782
# B_avg.bright.female  -0.06048393 0.2827917
summary(diff.sires.close50.plumF)
# lm(formula = diff.ep.sires ~ proportion.close.50 + B_avg.bright.female, 
#    data = dat4)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.0313 -0.4693 -0.1124  0.3978  1.2073 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)         -1.95094    1.91225  -1.020    0.342
# proportion.close.50 -4.84312    2.91754  -1.660    0.141
# B_avg.bright.female  0.11115    0.07259   1.531    0.170
# 
# Residual standard error: 0.8689 on 7 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.3782,	Adjusted R-squared:  0.2006 
# F-statistic: 2.129 on 2 and 7 DF,  p-value: 0.1895


# 2) movement * female plumage
Anova(diff.sires.close50.plumF.int, type="III")
# Anova Table (Type III tests)
# 
# Response: diff.ep.sires
# Sum Sq Df F value Pr(>F)
# (Intercept)                             0.3081  1  0.4547 0.5252
# proportion.close.50                     1.5167  1  2.2385 0.1852
# B_avg.bright.female                     0.1366  1  0.2016 0.6692
# proportion.close.50:B_avg.bright.female 1.2195  1  1.7998 0.2283
# Residuals                               4.0654  6  
confint(diff.sires.close50.plumF.int)
# 2.5 %     97.5 %
#   (Intercept)                               -6.8026449 11.9781338
# proportion.close.50                     -119.1799405 28.7364452
# B_avg.bright.female                       -0.4334587  0.2990487
# proportion.close.50:B_avg.bright.female   -1.2778357  4.3795934
summary(diff.sires.close50.plumF.int)
# lm(formula = diff.ep.sires ~ proportion.close.50 * B_avg.bright.female, 
#    data = dat4)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.85060 -0.47734 -0.04397  0.29482  1.31168 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                               2.5877     3.8377   0.674    0.525
# proportion.close.50                     -45.2218    30.2251  -1.496    0.185
# B_avg.bright.female                      -0.0672     0.1497  -0.449    0.669
# proportion.close.50:B_avg.bright.female   1.5509     1.1560   1.342    0.228
# 
# Residual standard error: 0.8231 on 6 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.5217,	Adjusted R-squared:  0.2826 
# F-statistic: 2.182 on 3 and 6 DF,  p-value: 0.1912



################################################################################
#-------------------------------------------------------------------------------
# Models for change in number of EP young


# 1) close50 + female plumage
diff.numEP.close50.plumF <- lm(diff.num.ep ~ proportion.close.50 + 
                                 B_avg.bright.female, data = dat4)
hist(diff.numEP.close50.plumF$residuals)




# 2) close50 * female plumage
diff.numEP.close50.plumF.int <- lm(diff.num.ep ~ proportion.close.50 * 
                                     B_avg.bright.female, data = dat4)
hist(diff.numEP.close50.plumF.int$residuals)


#-------------------------------------------------------------------------------
# Likelihood ratio tests for change in number EP young


# 1) close50 + female plumage
Anova(diff.numEP.close50.plumF)
# Anova Table (Type II tests)
# 
# Response: diff.num.ep
# Sum Sq Df F value Pr(>F)
# proportion.close.50 16.815  1  2.5782 0.1524
# B_avg.bright.female 20.006  1  3.0674 0.1233
# Residuals           45.654  7 
confint(diff.numEP.close50.plumF)
# 2.5 %    97.5 %
#   (Intercept)         -21.3899611 5.1903367
# proportion.close.50 -34.0458553 6.5080018
# B_avg.bright.female  -0.1308261 0.8781146
summary(diff.numEP.close50.plumF)
# lm(formula = diff.num.ep ~ proportion.close.50 + B_avg.bright.female, 
#    data = dat4)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1126 -2.0195  0.2224  1.4352  3.4275 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)          -8.0998     5.6204  -1.441    0.193
# proportion.close.50 -13.7689     8.5751  -1.606    0.152
# B_avg.bright.female   0.3736     0.2133   1.751    0.123
# 
# Residual standard error: 2.554 on 7 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.4024,	Adjusted R-squared:  0.2317 
# F-statistic: 2.357 on 2 and 7 DF,  p-value: 0.165


# 2) close50 * female plumage
Anova(diff.numEP.close50.plumF.int, type="III")
# Anova Table (Type III tests)
# 
# Response: diff.num.ep
# Sum Sq Df F value Pr(>F)
# (Intercept)                              0.174  1  0.0240 0.8819
# proportion.close.50                      3.484  1  0.4816 0.5137
# B_avg.bright.female                      0.525  1  0.0725 0.7967
# proportion.close.50:B_avg.bright.female  2.244  1  0.3101 0.5978
# Residuals                               43.410  6 
confint(diff.numEP.close50.plumF.int)
# 2.5 %    97.5 %
#   (Intercept)                              -32.628621  28.74170
# proportion.close.50                     -310.213887 173.13542
# B_avg.bright.female                       -1.065099   1.32853
# proportion.close.50:B_avg.bright.female   -7.139807  11.34709
summary(diff.numEP.close50.plumF.int)
# lm(formula = diff.num.ep ~ proportion.close.50 * B_avg.bright.female, 
#    data = dat4)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -3.151 -1.989  0.157  1.561  2.941 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                              -1.9435    12.5404  -0.155    0.882
# proportion.close.50                     -68.5392    98.7672  -0.694    0.514
# B_avg.bright.female                       0.1317     0.4891   0.269    0.797
# proportion.close.50:B_avg.bright.female   2.1036     3.7776   0.557    0.598
# 
# Residual standard error: 2.69 on 6 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.4318,	Adjusted R-squared:  0.1477 
# F-statistic:  1.52 on 3 and 6 DF,  p-value: 0.3027


