

# Heather Kenny-Duddela
# Feb 16, 2023

# first pass at paternity assignment from lcMLkin output

# set working directory
setwd("~/CU Boulder/BARS fieldwork/2022 Field and Lab/Paternity assignment methods/barn-swallow-gps/paternity-assignment")

# read in lcMLkin output file for 2021 data
kin <- read.table("input-data/CO-2021.thin10k.GL-recalc.relate", header=T)

# read in families table
fam.21 <- read.csv("input-data/fam_clutch21_updated broods.csv")


# load libraries
library(tidyverse)
library(ggplot2)


# k0_hat = prob no alleles shared that are IBD (unrelated)
# k1_hat = prob 1 allele shared that is IBD
# k2_hat = prob 2 alleles shared that are IBD
# pi_hat = coefficient of relatedness (r)

# coefficient of relatedness (r) = 2*phi = (k1)/2 + k2

# pull out one focal individual
# Blue Cloud 09
kin.83404 <- subset(kin, kin$Ind1=="CO_83404")

# try calculating r for the test individual
# see that the calculation for r matches pi_hat
kin.83404$r <- kin.83404$k1_hat/2 + kin.83404$k2_hat

# plot k0 vs. pi_hat
plot(kin.83404$k0_hat, kin.83404$pi_HAT)


# make plot for all individuals
plot(kin$k0_hat, kin$pi_HAT)

# more test birds
# Cathy's 02
kin.82984 <- subset(kin, kin$Ind1=="CO_82984")
plot(kin.82984$k0_hat, kin.82984$pi_HAT)

# cooks 23 - no father assigned, social dad unrelated
kin.82969 <- subset(kin, kin$Ind1=="CO_82969")
plot(kin.82969$k0_hat, kin.82969$pi_HAT )

kin.83018 <- subset(kin, kin$Ind1 == "CO_83018")
plot(kin.83018$k0_hat, kin.83018$pi_HAT)


# seems like good relatedness cut-off values are 
# <0.9 for k0 and 
# >0.1 for pi_hat


# Add band, site, family, sex, type to kinship table for ind1 and ind2
# this way we can filter to only parent-offspring relationships
# and can identify social vs. extra-pair dams and sires

# replace "-" in seqID with "_"
fam.21$Seq.ID.test <- gsub("-","_",fam.21$Seq.ID)


# relevant info from fam.21
# band, site, nest, type, familyID, brood, clutch_id, seq.ID.test
seq.cols <- fam.21[,c(2,3,5,6,9,10,11,12)]

colnames(seq.cols) <- paste(colnames(seq.cols),"ind1",sep="_")
colnames(seq.cols)[8] <- "Ind1"

kin1 <- left_join(kin, seq.cols, by="Ind1")

seq.cols2 <- fam.21[,c(2,3,5,6,9,10,11,12)]

colnames(seq.cols2) <- paste(colnames(seq.cols2),"ind2",sep="_")
colnames(seq.cols2)[8] <- "Ind2"

kin2 <- left_join(kin1, seq.cols2, by="Ind2")



# filter so ind1 is all kids
kin.kids <- subset(kin2, kin2$Type_ind1=="kid" | kin2$Type_ind1=="egg")

# next filter so ind2 is only parents (due to the way the pairwise combinations
# were made, it didn't work to filter ind1 as parents and then ind2 as kids)
kin.kids.par <- subset(kin.kids, kin.kids$Type_ind2=="dad" | kin.kids$Type_ind2=="mom")

# check to make sure I got all possible kid-parent combinations from the table
# this object is the same size as kin.kids.par
kin.po <- subset(kin2, kin2$Type_ind1=="kid" & kin2$Type_ind2=="mom" |
                     kin2$Type_ind1=="kid" & kin2$Type_ind2=="dad" |
                     kin2$Type_ind1=="egg" & kin2$Type_ind2=="mom" |
                     kin2$Type_ind1=="egg" & kin2$Type_ind2=="dad" |
                     kin2$Type_ind1=="mom" & kin2$Type_ind2=="egg" |
                     kin2$Type_ind1=="mom" & kin2$Type_ind2=="kid" |
                     kin2$Type_ind1=="dad" & kin2$Type_ind2=="egg" |
                     kin2$Type_ind1=="dad" & kin2$Type_ind2=="kid")

# number of unique kids = 150
length(unique(kin.po$Ind1))

# pull out some test dads

dad.06912 <- subset(kin.kids.par, kin.kids.par$Ind2=="CO_06912")
plot(dad.06912$k0_hat, dad.06912$pi_HAT) # this guy has only 2 offspring

dad.57505 <- subset(kin.kids.par, kin.kids.par$Ind2=="CO_57505")
plot(dad.57505$k0_hat, dad.57505$pi_HAT) # this guy has 12 offspring (including collected eggs)


# filter so k0<0.9 and pi_hat > 0.1
kin.close <- subset(kin.kids.par, kin.kids.par$k0_hat<0.9 & kin.kids.par$pi_HAT>0.1)

# check that test dads retain the expected offspring
dad.assign.06912 <- subset(kin.close, kin.close$Ind2=="CO_06912")

dad.assign.57505 <- subset(kin.close, kin.close$Ind2=="CO_57505")


# Identify cases where parent familyID does not match kid familyID (extra-pair sire)
index <- which(kin.close$FamilyID_ind1 != kin.close$FamilyID_ind2)

# look at mismatched rows
mismatch <- kin.close[index,]

# some mismatched parents have pi_hat values lower than 0.3
# I think this is because they are grandparents or aunts/uncles instead of parents
plot(mismatch$k0_hat, mismatch$pi_HAT)


# try subsetting to test specific sib/PO relationships between different adults
# first get table of just moms and dads
adults <- subset(kin2, kin2$Type_ind1!="egg" & kin2$Type_ind1!="kid")
# check another way (same size table)
adults2 <- subset(kin2, kin2$Type_ind1=="mom" & kin2$Type_ind2=="dad" |
                    kin2$Type_ind1=="mom" & kin2$Type_ind2=="mom" |
                    kin2$Type_ind1=="dad" & kin2$Type_ind2=="mom" |
                    kin2$Type_ind1=="dad" & kin2$Type_ind2=="dad")


# now filter for close relatives
# some of the low-value relationships in mismatch are explained, but not all...
adults.close <- subset(adults, adults$k0_hat<0.9 & adults$pi_HAT>0.1)


mismatch.site <- kin.close[which(kin.close$Site_ind1 != kin.close$Site_ind2), ]

# look at matched rows (moms and within-pair sires)
match <- kin.close[which(kin.close$FamilyID_ind1 == kin.close$FamilyID_ind2), ]

# all matched parents have pi_hat values greater than 0.3
plot(match$k0_hat, match$pi_HAT)
# save minimum relatedness value for matched parents
match.par.min <- min(match$pi_HAT)

# range of matched mom pi_hat values
# 0.367 to 0.557
range(subset(match$pi_HAT, match$Type_ind2=="mom"))

match.mom <- subset(match, match$Type_ind2=="mom")

mom.kid.relate.min <- min(match.mom$pi_HAT) #0.367
# subset to only include relatedness values greater than or equal to the mom min
kin.close.cutoff <- subset(kin.close, kin.close$pi_HAT >= 0.367)

# check how many kids have both, one, and no parents assigned
assign.par <- kin.close.cutoff %>% group_by(Ind1) %>%
  summarise(mom.assign = length(which(Type_ind2=="mom")),
            dad.assign = length(which(Type_ind2 == "dad")))

# check mom and dad assignments for the larger close relations table
# not using the minimum relatedness value cutoff
assign.par.check <- kin.close %>% group_by(Ind1) %>%
  summarise(mom.assign = length(which(Type_ind2=="mom")),
            dad.assign = length(which(Type_ind2 == "dad")))

# discrepancy between the cutoff close relationships and the regular close
# 25 unassigned with close, 28 unassigned with cutoff
# figure out which offspring lost the dad assignment

# 25 offspring with unassigned dads
length(which(assign.par.check$dad.assign==0))

# 28 unassigned from the cutoff table
length(which(assign.par$dad.assign==0))

# combine cutoff and non-cutoff tables to identify kids that are assigned in 
# one but unassigned in the other table
assign.par.check2 <- left_join(assign.par, assign.par.check[c(1,3)], by="Ind1")
# subset cases where assignments are mismatched between the two tables
assign.par.check3 <- subset(assign.par.check2, assign.par.check2$dad.assign.x !=
                              assign.par.check2$dad.assign.y)

# check relatedness values for newly unassigned kids
kin.close.83214 <- subset(kin.close, kin.close$Ind1=="CO_83214")


# pull out offspring with unassigned dads
dad.unk.index <- which(kin.close.cutoff$Ind1 %in% assign.par$Ind1[which(assign.par$dad.assign==0)])

dad.unk <- kin.close.cutoff[dad.unk.index,]

# the 28 unassigned offspring are from 7 different families
length(unique(dad.unk$FamilyID_ind1))



kin.cutoff.21 <- kin.close %>%
  group_by(Ind1, Type_ind2) %>%
  summarise(par.cutoff = Ind2[which(pi_HAT >= 0.367)]) # cutoff from matched parent-offspring relatedness values from 2021

colnames(kin.cutoff.21)[3] <- "Ind2"

kin.cutoff.21.info <- left_join(kin.cutoff.21[,c(1,3)], kin.close, by=c("Ind1","Ind2"))



# check number of kids with assigned moms and dads
# now there are 28 with no dad assigned,  and all have one mom assigned
assign.par.cutoff <- kin.cutoff.21.info %>% group_by(Ind1) %>%
  summarise(mom.assign = length(which(Type_ind2=="mom")),
            dad.assign = length(which(Type_ind2 == "dad")))

# make new unassigned dads table
dad.unk.index.cutoff <- which(kin.cutoff.21.info$Ind1 %in% 
                                assign.par.cutoff$Ind1[which(assign.par.cutoff$dad.assign==0)])

dad.unk.cutoff <- kin.cutoff.21.info[dad.unk.index.cutoff,]

# check number of clutches for unassigned kids (7)
length(unique(dad.unk.cutoff$FamilyID_ind1))


# check again for extra-pair sires after filtering extra dads
kin.ep <- kin.cutoff.21.info[which(kin.cutoff.21.info$FamilyID_ind1 != kin.cutoff.21.info$FamilyID_ind2), ]

# The only true EP dad between sites was Blue Cloud 18 (band 2640-97548) who 
# sired offspring at Struthers and Make Believe
# Interesting because his social mate was GPS tagged and then she disappeared and was never recaptured
# He had no social offspring to take care of so was probably pursuing lots of EP matings
# Hepp and Cooks are so close together, it is not clear that they should be treated as different sites
kin.ep.site <- kin.cutoff.21.info[which(kin.cutoff.21.info$FamilyID_ind1 != 
                                          kin.cutoff.21.info$FamilyID_ind2 &
                                   kin.cutoff.21.info$Site_ind1 != 
                                     kin.cutoff.21.info$Site_ind2), ]

kin.hi.info <- kin.cutoff.21.info

# check number of assigned moms and dads
length(unique(subset(kin.cutoff.21.info$Band_ind2, kin.cutoff.21.info$Type_ind2=="dad")))

length(unique(subset(kin.cutoff.21.info$Band_ind2, kin.cutoff.21.info$Type_ind2=="mom")))

# number of unique clutches
length(unique(kin.cutoff.21.info$clutch_id_ind1))

# number of males sampled
length(unique(c(subset(kin2$Ind1, kin2$Type_ind1=="dad"),
                subset(kin2$Ind2, kin2$Type_ind2=="dad"))))

# number of females sampled
length(unique(c(subset(kin2$Ind1, kin2$Type_ind1=="mom"),
                subset(kin2$Ind2, kin2$Type_ind2=="mom"))))



# check number of genetic offspring per mom, and number of genetic offspring per known dad
offspring.by.par <- kin.hi.info %>%
  group_by(Ind2, Type_ind2, FamilyID_ind2, Site_ind2) %>%
  summarise(kids.genetic = n( ),
            within.pair = length(which(FamilyID_ind2 == FamilyID_ind1)),
            extra.pair = length(which(FamilyID_ind2 != FamilyID_ind1)),
            same.site.ep = length(which(FamilyID_ind2 != FamilyID_ind1 &
                                          Site_ind1 == Site_ind2)),
            diff.site.ep = length(which(FamilyID_ind2 != FamilyID_ind1 &
                                          Site_ind1 != Site_ind2)))



####################################################################################
# For kids with unassigned dads, check whether they are full or half sibs

# take indiv from full kin table that are in the unassigned dads table
dad.unk.sibs <- kin2[which(kin2$Ind1 %in% dad.unk.cutoff$Ind1), ]

# filter to keep unassigned kids compared to unassigned kids
dad.unk.sibs2 <- dad.unk.sibs[which(dad.unk.sibs$Ind2 %in% dad.unk.cutoff$Ind1), ]

# filter to keep only pairs that are within the same family
dad.unk.sibs.fam <- subset(dad.unk.sibs2, dad.unk.sibs2$FamilyID_ind1 == dad.unk.sibs2$FamilyID_ind2)



# solution proposed by Adi Schuerg (Toni's brother)
# Convert pairwise list into distance matrix for each clutch
# Set element for full sibs = 0 and element for half sibs = 1
# Then find the product of each row, and add them up to get the number of fathers.
# This works because the presence of a full sib sets the row to zero and reduces the
# total number of possible fathers by one. 


#-------------------------------------------------------------------------------
# try doing this for Cooks-23 and Cooks-31 to see if share the same unk dad

# change Hepp site to "Cooks" because the two sites are so close together
dad.unk.sibs3 <- dad.unk.sibs2
dad.unk.sibs3$Site_ind1[which(dad.unk.sibs3$Site_ind1=="Hepp (near Cook)")] <- "Cooks"
dad.unk.sibs3$Site_ind2[which(dad.unk.sibs3$Site_ind2=="Hepp (near Cook)")] <- "Cooks"

# from the unassigned kids compared with unassigned kids table, keep matches from the same site
dad.unk.sibs.site <- subset(dad.unk.sibs3, dad.unk.sibs2$Site_ind1 == dad.unk.sibs3$Site_ind2)

# relatedness values for dad.unk.sibs.site has a clear cutoff, where it jumps 
# from 0.155 to 0.020
# k_0 also jumps from 0.8 to 0.964
plot(dad.unk.sibs.site$k0_hat, dad.unk.sibs.site$pi_HAT)

cooks.test <- subset(
  dad.unk.sibs.site, dad.unk.sibs.site$FamilyID_ind1 == "Cooks-23" |
                       dad.unk.sibs.site$FamilyID_ind1 == "Cooks-31" |
                       dad.unk.sibs.site$FamilyID_ind2 == "Cooks-23" |
                       dad.unk.sibs.site$FamilyID_ind2== "Cooks-31")

#-------------------------------------------------------------------------------
# Calculate number of unknown sires per clutch for collected clutch only

collect <- subset(dad.unk.sibs.fam, dad.unk.sibs.fam$Brood_ind1 == "collected" &
                    dad.unk.sibs.fam$Brood_ind2 == "collected")

unk.fam.collect <- as.data.frame(unique(c(collect$FamilyID_ind1, collect$FamilyID_ind2)))
colnames(unk.fam.collect)[1] <- "FamilyID"
# add column to store number of dads
unk.fam.collect$num.dads <- NA

for (i in 1:length(unk.fam.collect$FamilyID)) {
  # subset by family id
  fam <- subset(collect, collect$FamilyID_ind1 == unk.fam.collect$FamilyID[i] |
                  collect$FamilyID_ind2 == unk.fam.collect$FamilyID[i])
  
  # get relatedness data
  rel <- fam[,c(1,2,6)]
  
  # use reshape to go from a pairwise list to a distance matrix
  rel.r <- reshape(rel, direction="wide", idvar="Ind2", timevar="Ind1")
  
  # remove first column which is just labels
  rel.r.mat <- as.matrix(rel.r[,-1]) 
  
  # give full sibs as 0
  rel.r.bin <- rel.r.mat
  rel.r.bin[rel.r.mat >= 0.239] <- 0 # cutoff value from range of relatedness from genetic full sibs in 2022
  
  # give half sibs as 1
  rel.r.bin[rel.r.mat < 0.239] <- 1
  
  # calculate row products
  row.prod <- apply(rel.r.bin, 1, prod, na.rm=T)
  
  # sum up and add 1 to get the number of fathers
  num.dad <- sum(row.prod) + 1
  
  # save in storage
  unk.fam.collect$num.dads[i] <- num.dad
  
}



#-------------------------------------------------------------------------------
# Calculate number of unknown sires per clutch for replacement clutch only

replace <- subset(dad.unk.sibs.fam, dad.unk.sibs.fam$Brood_ind1 == 1 &
                    dad.unk.sibs.fam$Brood_ind2 ==1)

unk.fam.replace <- as.data.frame(unique(c(replace$FamilyID_ind1, replace$FamilyID_ind2)))
colnames(unk.fam.replace)[1] <- "FamilyID"
# add column to store number of dads
unk.fam.replace$num.dads <- NA

for (i in 1:length(unk.fam.replace$FamilyID)) {
  # subset by family id
  fam <- subset(replace, replace$FamilyID_ind1 == unk.fam.replace$FamilyID[i] |
                  replace$FamilyID_ind2 == unk.fam.replace$FamilyID[i])
  
  # get relatedness data
  rel <- fam[,c(1,2,6)]
  
  # use reshape to go from a pairwise list to a distance matrix
  rel.r <- reshape(rel, direction="wide", idvar="Ind2", timevar="Ind1")
  
  # remove first column which is just labels
  rel.r.mat <- as.matrix(rel.r[,-1]) 
  
  # give full sibs as 0
  rel.r.bin <- rel.r.mat
  rel.r.bin[rel.r.mat >= 0.239] <- 0 # cutoff value from range of relatedness from genetic full sibs in 2022
  
  # give half sibs as 1
  rel.r.bin[rel.r.mat < 0.239] <- 1
  
  # calculate row products
  row.prod <- apply(rel.r.bin, 1, prod, na.rm=T)
  
  # sum up an add 1 to get the number of fathers
  num.dad <- sum(row.prod) + 1
  
  # save in storage
  unk.fam.replace$num.dads[i] <- num.dad
  
}

#-------------------------------------------------------------------------------
# count number of unk sires in 2nd broods

brood2 <- subset(dad.unk.sibs.fam, dad.unk.sibs.fam$Brood_ind1 == 2 &
                    dad.unk.sibs.fam$Brood_ind2 ==2)

unk.fam.brood2 <- as.data.frame(unique(c(brood2$FamilyID_ind1, brood2$FamilyID_ind2)))
colnames(unk.fam.brood2)[1] <- "FamilyID"
# add column to store number of dads
unk.fam.brood2$num.dads <- NA

for (i in 1:length(unk.fam.brood2$FamilyID)) {
  # subset by family id
  fam <- subset(brood2, brood2$FamilyID_ind1 == unk.fam.brood2$FamilyID[i] |
                  brood2$FamilyID_ind2 == unk.fam.brood2$FamilyID[i])
  
  # get relatedness data
  rel <- fam[,c(1,2,6)]
  
  # use reshape to go from a pairwise list to a distance matrix
  rel.r <- reshape(rel, direction="wide", idvar="Ind2", timevar="Ind1")
  
  # remove first column which is just labels
  rel.r.mat <- as.matrix(rel.r[,-1]) 
  
  # give full sibs as 0
  rel.r.bin <- rel.r.mat
  rel.r.bin[rel.r.mat >= 0.239] <- 0 # cutoff value from range of relatedness from genetic full sibs in 2022
  
  # give half sibs as 1
  rel.r.bin[rel.r.mat < 0.239] <- 1
  
  # calculate row products
  row.prod <- apply(rel.r.bin, 1, prod, na.rm=T)
  
  # sum up an add 1 to get the number of fathers
  num.dad <- sum(row.prod) + 1
  
  # save in storage
  unk.fam.brood2$num.dads[i] <- num.dad
  
}

#-------------------------------------------------------------------------------
# calculate number of dads for each female across the whole season

### Loop to calculate number of unknown dads within each family (within each mom)

# object for number of families (female) with unassigned kids
unk.fam <- as.data.frame(unique(dad.unk.sibs.fam$FamilyID_ind1))
test <-  as.data.frame(unique(c(dad.unk.sibs.fam$FamilyID_ind1, 
                                dad.unk.sibs.fam$FamilyID_ind2)))
colnames(unk.fam)[1] <- "FamilyID"
# add column to store number of dads
unk.fam$num.dads <- NA

for (i in 1:length(unk.fam$FamilyID)) {
  # subset by family id
  fam <- subset(dad.unk.sibs.fam, dad.unk.sibs.fam$FamilyID_ind1 == unk.fam$FamilyID[i])
  
  # get relatedness data
  rel <- fam[,c(1,2,6)]
  
  # use reshape to go from a pairwise list to a distance matrix
  rel.r <- reshape(rel, direction="wide", idvar="Ind2", timevar="Ind1")
  
  # remove first column which is just labels
  rel.r.mat <- as.matrix(rel.r[,-1]) 
  
  # give full sibs as 0
  rel.r.bin <- rel.r.mat
  rel.r.bin[rel.r.mat >= 0.239] <- 0 # cutoff value from range of relatedness from genetic full sibs in 2022
  
  # give half sibs as 1
  rel.r.bin[rel.r.mat < 0.239] <- 1
  
  # calculate row products
  row.prod <- apply(rel.r.bin, 1, prod, na.rm=T)
  
  # sum up an add 1 to get the number of fathers
  num.dad <- sum(row.prod) + 1
  
  # save in storage
  unk.fam$num.dads[i] <- num.dad
  
}



#-------------------------------------------------------------------------------
# Find number of shared ferts for each genetic fam
#-------------------------------------------------------------------------------

# switch to wide format where each row is a kid, and there are cols for dad and mom
# first just do dads
po.dad.wide21 <- subset(kin.hi.info, kin.hi.info$Type_ind2=="dad")
# pull out only cols that need to get "dad" labels added
po.dad.wide21.2 <- po.dad.wide21[,c(1,2,3:7,15:18)]
colnames(po.dad.wide21.2)[2:11] <- c("SeqID_dad", "k0_hat_dad", "k1_hat_dad", "k2_hat_dad", "pi_HAT_dad",
                                  "nbSNP_dad", "Band_dad","Site_dad", "Type_dad" , "FamilyID_dad")

# do same for moms
po.mom.wide21 <- subset(kin.hi.info, kin.hi.info$Type_ind2=="mom")
po.mom.wide21.2 <- po.mom.wide21[,c(1,2,3:7,15:18)]
colnames(po.mom.wide21.2)[2:11] <- c("SeqID_mom","k0_hat_mom", "k1_hat_mom", "k2_hat_mom", "pi_HAT_mom",
                                  "nbSNP_mom", "Band_mom","Site_mom", "Type_mom" , "FamilyID_mom")

# pull out columns with relevant kid info for kids with assigned parents, and no duplicated kids
# Ind1 are kids, Ind2 are parents
po.kids.21 <- kin.hi.info[-(which(duplicated(kin.hi.info$Ind1))),c(1,8:14)]

po.wide.21.m <- left_join(po.kids.21, po.mom.wide21.2, by=c("Ind1"))

po.wide.21 <- left_join(po.wide.21.m, po.dad.wide21.2, by="Ind1")

po.wide.21$genetic_fam <- paste(po.wide.21$Band_dad, po.wide.21$Band_mom, sep="_")

# move genetic family label near the front of the table
po.wide.21 <- po.wide.21[,c(1, 29, 2:28)]


# save wide file
write.csv(po.wide.21, file="output-files/kin_2021_parent_offspring_assigned_wide.csv")


# count up number of offspring produced by each genetic family
shared.fert.list.21 <- po.wide.21 %>% 
  group_by(genetic_fam, Band_dad, Band_mom, 
           Site_dad, Site_ind1, FamilyID_dad, FamilyID_mom) %>%
  summarise(fert = n())

hist(shared.fert.list.21$fert)

# label fertilization types
shared.fert.list.21$fert_type <- NA
shared.fert.list.21$fert_type[which(shared.fert.list.21$FamilyID_dad == 
                                   shared.fert.list.21$FamilyID_mom)] <- "wp"
shared.fert.list.21$fert_type[which(shared.fert.list.21$FamilyID_dad != 
                                     shared.fert.list.21$FamilyID_mom &
                                   shared.fert.list.21$Site_dad == 
                                     shared.fert.list.21$Site_ind1)] <- "ep_same"
shared.fert.list.21$fert_type[which(shared.fert.list.21$FamilyID_dad != 
                                     shared.fert.list.21$FamilyID_mom &
                                   shared.fert.list.21$Site_dad != 
                                     shared.fert.list.21$Site_ind1)] <- "ep_diff"
shared.fert.list.21$fert_type[which(is.na(shared.fert.list.21$FamilyID_mom))] <- "mom_unk"
shared.fert.list.21$fert_type[which(is.na(shared.fert.list.21$FamilyID_dad))] <- "dad_unk"


# save shared fertilizations list
write.csv(shared.fert.list.21, file="output-files/fert_2021_by_genetic_family.csv", row.names=F)


# summarise number and proportion of each offspring type
total.fert.21 <- sum(shared.fert.list.21$fert)

fert.type.summary.21 <- shared.fert.list.21 %>%
  group_by(fert_type) %>%
  summarise(num.fert = sum(fert),
            prop.fert = num.fert/total.fert.21)

# summarize fertilization types within nest sites
fert.site.sum.21<- shared.fert.list.21 %>%
  group_by(Site_ind1, fert_type) %>%
  summarise(num.fert = sum(fert))

fert.site.sum.21$fert_type <- factor(fert.site.sum.21$fert_type, levels=c("wp","ep_same","mom_unk",
                                                                    "dad_unk","ep_diff"))




#-------------------------------------------------------------------------------
# at the nest/clutch level, look at proportion of outside site dads
# get a sense for whether patterns are usually lots of EP at few nests, or 
# some EP at most nests
#-------------------------------------------------------------------------------


# summarise at clutch level
shared.fert.clutch.21 <- po.wide.21 %>% 
  group_by(clutch_id_ind1, Band_dad, Band_mom, 
           Site_dad, Site_ind1, FamilyID_dad, FamilyID_mom, Brood_ind1) %>%
  summarise(fert = n())

# number of clutches
length(unique(shared.fert.clutch.21$clutch_id_ind1))

# assign fert categories
shared.fert.clutch.21$fert_type <- NA
shared.fert.clutch.21$fert_type[which(shared.fert.clutch.21$FamilyID_dad == 
                                     shared.fert.clutch.21$FamilyID_mom)] <- "wp"

shared.fert.clutch.21$fert_type[which(shared.fert.clutch.21$FamilyID_dad != 
                                     shared.fert.clutch.21$FamilyID_mom &
                                     shared.fert.clutch.21$Site_dad == 
                                     shared.fert.clutch.21$Site_ind1)] <- "ep_same"

shared.fert.clutch.21$fert_type[which(shared.fert.clutch.21$FamilyID_dad != 
                                        shared.fert.clutch.21$FamilyID_mom &
                                     shared.fert.clutch.21$Site_dad != 
                                       shared.fert.clutch.21$Site_ind1)] <- "ep_diff"

shared.fert.clutch.21$fert_type[which(is.na(shared.fert.clutch.21$FamilyID_mom))] <- "mom_unk"
shared.fert.clutch.21$fert_type[which(is.na(shared.fert.clutch.21$FamilyID_dad))] <- "dad_unk"


# summarise number and proportion of each offspring type
total.fert.clutch.21 <- sum(shared.fert.clutch.21$fert)

# number of females where we have both 1st and 2nd broods
fem.both.broods.21 <- shared.fert.clutch.21 %>%
  group_by(FamilyID_mom, Brood_ind1) %>%
  summarise(chicks = n())


fem.both.brood2.21 <- fem.both.broods.21 %>%
  group_by(FamilyID_mom) %>%
  summarise(brood_count = n())

# 14 females where we have more than one brood
length(subset(fem.both.brood2.21$FamilyID_mom, fem.both.brood2.21$brood_count>=2))  

# proportion of within-pair offspring per clutch
shared.fert.clutch2.21 <- shared.fert.clutch.21 %>% 
  group_by(clutch_id_ind1) %>%
  mutate(clutch_size = sum(fert))

# specify "collected" as brood number for collected eggs
shared.fert.clutch2.21$Brood_ind1[which(is.na(shared.fert.clutch2.21$Brood_ind1))] <- "collected"

# save fertilization types by clutch ID
write.csv(shared.fert.clutch2.21, file="output-files/fert_2021_by_clutchID.csv", row.names=F)
# calculate proportion of within-pair offspring per clutch
shared.fert.wp.21 <- subset(shared.fert.clutch2.21, shared.fert.clutch2.21$fert_type=="wp")
shared.fert.wp.21 <- shared.fert.wp.21 %>%
  mutate(prop.wp = fert/clutch_size)
# add proportion wp to the main fertilization table with all fert types
shared.fert.wp2.21 <- left_join(shared.fert.clutch2.21, shared.fert.wp.21[,c(1,12)], by="clutch_id_ind1")
# for clutches where there were no wp offspring, set prop_wp to zero
shared.fert.wp2.21$prop.wp[which(is.na(shared.fert.wp2.21$prop.wp))] <- 0

# reduce to one row per clutch, for plotting proportion WP
clutch.dup <- duplicated(shared.fert.wp2.21$clutch_id_ind1)
shared.fert.prop <- shared.fert.wp2.21[which(clutch.dup==FALSE),]

ggplot(shared.fert.prop, aes(x=prop.wp)) + 
  geom_histogram(fill="#440154FF", color="black", bins = 16) +
  xlab("Proportion of within-pair offspring") + ylab("Count") + 
  ggtitle("Distirbution of within-pair proportions per clutch 2021")



# plot across all nests

# change factor order to match previous plots
shared.fert.clutch2.21$fert_type <- factor(shared.fert.clutch2.21$fert_type, 
                                          levels=c("wp","ep_same","mom_unk",
                                                   "dad_unk","ep_diff"))

# reorder brood factors
shared.fert.clutch2.21$Brood_ind1 <- factor(shared.fert.clutch2.21$Brood_ind1, 
                                           levels=c("collected","1","2"))

# facet by brood 1 and 2 and collected
ggplot(shared.fert.clutch2.21, aes(x=FamilyID_mom, y=fert, fill=fert_type)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_viridis_d() +
  xlab("Family ID") + ylab("Number of\n fertilizations") + 
  facet_grid(Brood_ind1~.) +
  ggtitle("Fertilization types across broods 2021")


#-------------------------------------------------------------------------------
# Determine total number of mates for each female
#-------------------------------------------------------------------------------

# count number of mates per female
mates.by.fem <- shared.fert.clutch2.21 %>%
  group_by(Band_mom, FamilyID_mom) %>%
  summarise(num.mates = length(unique(Band_dad)),
            num.clutches = length(unique(clutch_id_ind1)),
            num.chicks = sum(fert))


# add in number of unknown mates
mates.by.fem.unk <- mates.by.fem

# match by family ID from unknown mates table
# subtract one from the number in the unk.fam table because the mates.by.fem
# table already includes one unknown sire from the NA present in the Band_dad
# column of the shared.fert.clutch2.21 table. So a female with offspring from 
# two unknown sires just needs one additional sire added to her total number
# of mates. 
mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[1])] <-
  mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[1])] +
  unk.fam$num.dads[1] -1

mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[2])] <-
  mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[2])] +
  unk.fam$num.dads[2] -1

mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[3])] <-
  mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[3])] +
  unk.fam$num.dads[3] -1

mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[4])] <-
  mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[4])] +
  unk.fam$num.dads[4] -1

mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[5])] <-
  mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[5])] +
  unk.fam$num.dads[5] -1

mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[6])] <-
  mates.by.fem.unk$num.mates[which(mates.by.fem.unk$FamilyID_mom == unk.fam$FamilyID[6])] +
  unk.fam$num.dads[6] -1


# save mates per female table
write.csv(mates.by.fem.unk, file="output-files/num mates chicks clutches per female 2021.csv", row.names=F)

#-------------------------------------------------------------------------------
# Now calculate total number of mates per clutch (rather than across the whole season)

# count number of mates per female per clutch
mates.by.fem.clutch <- shared.fert.clutch2.21 %>%
  group_by(Band_mom, FamilyID_mom, clutch_id_ind1, Brood_ind1) %>%
  summarise(num.mates = length(unique(Band_dad)),
            num.clutches = length(unique(clutch_id_ind1)),
            num.chicks = sum(fert))

# for collected clutches, no females had more than one UNK sire, so don't need to
# change the numbers counted in mates.by.fem.clutch
unk.fam.collect

mates.by.fem.clutch.collect <- subset(mates.by.fem.clutch, mates.by.fem.clutch$Brood_ind1=="collected")
mean(mates.by.fem.clutch.collect$num.mates)

# calculate mean EP mates in collected clutch
collect.ep <- mates.by.fem.clutch.collect$num.mates -1

# save collected clutch table
write.csv(mates.by.fem.clutch.collect, "output-files/num mates collected clutch.csv", row.names=F)

# for replacement clutch, do need to update
# subset to only replacement clutch
mates.by.fem.clutch.replace <- subset(mates.by.fem.clutch, mates.by.fem.clutch$Brood_ind1==1)

# only 2 replacement clutches had additional unknown dads, from unk.fam.replace
unk.fam.replace
# Cooks-23 had 2 unk dads total, so add 1 to the number of mates

mates.by.fem.clutch.replace.unk <- mates.by.fem.clutch.replace
mates.by.fem.clutch.replace.unk$num.mates[2] <- 3

# save mates in replacement clutch table
write.csv(mates.by.fem.clutch.replace.unk, "output-files/num mates replacement clutch.csv", row.names=F)

# for 2nd broods, Schaaps 80 had 2 UNK sires
unk.fam.brood2

# subset 2nd broods
mates.by.fem.clutch.brood2 <- subset(mates.by.fem.clutch, mates.by.fem.clutch$Brood_ind1==2)

mates.by.fem.clutch.brood2.unk <- mates.by.fem.clutch.brood2
mates.by.fem.clutch.brood2.unk$num.mates[6] <- 3

# updated table for number of mates per brood
mates.by.fem.clutch.update <- rbind(mates.by.fem.clutch.collect, 
                                    mates.by.fem.clutch.replace.unk,
                                    mates.by.fem.clutch.brood2.unk)

#-------------------------------------------------------------------------------
# plot number of mates per clutch for each female
#-------------------------------------------------------------------------------


# change factor order for brood
mates.by.fem.clutch.update$Brood_ind1 <- factor(mates.by.fem.clutch.update$Brood_ind1, 
                                     levels=c("collected","1","2"))

# plot number of mates over time by female
ggplot(mates.by.fem.clutch.update, aes(x=Brood_ind1, y=num.mates)) + ylim(0.5,3.5) +
  geom_point(aes(size=num.chicks/2)) +
  geom_line(aes(group=FamilyID_mom)) + facet_wrap("FamilyID_mom") +
  ggtitle("Mates per brood for 2021 females") +
  ylab("Number of mates") + xlab("Brood") 





