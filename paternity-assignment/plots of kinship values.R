
################################################################################
# Plots of kinship values
# Heather Kenny-Duddela
# Nov 19, 2024
################################################################################

# set working directory
setwd("~/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/paternity-assignment")

# load data
kin <- read.csv("output-files/kin_2021_parent_offspring_assigned_wide.csv")

# libraries
library(ggplot2)
library(ggbreak)

# Histograms of kinship values of moms
ggplot(kin, aes(x=pi_HAT_mom)) + 
  geom_histogram(fill="maroon", color="black", binwidth=0.005) +
  ggtitle("Relatedness values of mothers with offspring") +
  xlab("Estimated kinship coefficients between mothers and offspring")

# check the kinship value higher than 0.5, an egg from BlueCloud-09
kin.check <- subset(kin, kin$pi_HAT_mom > 0.5)

kin.check2 <- subset(kin, kin$FamilyID_ind1=="BlueCloud-09")


# make long table with mom and dad kinship values
kin.mom <- kin[,c(2,15)]
kin.mom$parent <- "Mothers"
colnames(kin.mom)[2] <- "pi_HAT"

kin.dad <- kin[,c(2,25)]
kin.dad$parent <- "Sires"
colnames(kin.dad)[2] <- "pi_HAT"

kin.long <- rbind(kin.mom, kin.dad)

# set colors to create a legend
cols <- c("Mothers"="maroon", "Sires"="orange")

# plot with facet
ggplot(kin.long, aes(x=pi_HAT)) + 
  geom_histogram(aes(fill=parent), color="black", binwidth=0.005) +
  ggtitle("Relatedness values of parents with offspring") +
  xlab("Estimated kinship coefficients between parents and offspring") +
  xlim(0.36, 0.56) +
  scale_x_break(c(0.44,0.54), scales = 0.1, ticklabels=c(0.55)) +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank()) +
  scale_fill_manual(name="Legend", values=cols) +
  facet_grid(parent~.)

# Figure S1
ggsave("output-files/kinship values parents offspring.png", h=3, w=5)

# t-test
t.test(pi_HAT ~ parent, data=kin.long, alternative="two.sided")
# Welch Two Sample t-test
# 
# data:  pi_HAT by parent
# t = -0.33242, df = 263.41, p-value = 0.7398
# alternative hypothesis: true difference in means between group Mothers and group Sires is not equal to 0
# 95 percent confidence interval:
#   -0.004308313  0.003063723
# sample estimates:
#   mean in group Mothers   mean in group Sires 
# 0.3971400             0.3977623







