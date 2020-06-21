###########################################################################
###########################################################################
############### SEX DIFFERENCES IN THE TIMSS 2015 MATH DATA ###############
###########################################################################
###########################################################################


#### Table of Contents ####

# [1] Read and Format Data
# dataframes and variables needed throughout the script

# [2] Weighted Functions
# functions needed throughout the script, with descriptions

# [3] Means and Medians, Tail Proportion Ratios (TPRs), U3 Ratios (U3Rs)
# Means and medians of total group, females, males, as well as mean and median differences
# TPRs and LTPRs above mean and percentiles
# LTPR tail-center differences
# U3Rs and LU3Rs above percentiles
# LU3R tail-center differences

# [4] Other Effect Sizes
# Cohen's d
# U3
# Probability of superiority (PS)
# Variance ratio (VR), Log-transformed VR (LVR)
# Left and right VR and LVR (VR_L, VR_R, LVR_L, LVR_R)
# Mean absolute deviation (from the median) ratio (MADR), Log-transformed MADR (LMADR)
# Left and right MADR and LMADR (MADR_L, MADR_R, LMADR_L, LMADR_R)
# Gini's mean difference ratio (GMDR), Log-transformed GMDR (LGMDR)

# [5] Other Effect Sizes Adjusted for Age
# dataframe with scores linearly corrected for age
# age-corrected effect sizes

# [6] Standard Errors
# jackknife resampling variance
# imputation variance
# total variance and standard error
# for all of the following
# Means of total group, females, males and mean difference
# LTPRs and LTPR tail-center differences**
# LU3Rs and LU3R tail-center differences*
# d, U3, LVR, LVR_L, LVR_R, LMADR, LMADR_L, LMADR_R, LGMDR
# Probability of superiority***

# [7] 95% Confidence Intervals
# 95% CIs for Means, TPRs, U3Rs, other effect sizes

# [8] Output
# write dataframe with all needed variables

# ratios are log-transformed to place them on a linear scale

# *nontrivial runtime (more than a few seconds on a personal computer)
# a higher number of asterisks indicates longer runtime
# updates are printed while long loops are running




################################
##### Read and Format Data #####
################################

# in the selected file ASGAUSM6.sav, AUS is the country code for Australia
# see pp. 47-49 of the "TIMSS 2015 User Guide for the International Database" for country codes
# the initial 'A' in ASGAUSM6.sav is instead 'B' for 8th grade files
# the initial 'A' in columns ASDAGE, ASDMLOWP, ASMMAT01, etc. is instead 'B' for 8th grade files
# Armenia is in the database but is not in the report due to late testing

library(haven) # read SPSS

CNT <- 'AUS' # enter country code manually

Grade <- 4 # enter grade manually

TIMSS15 <- read_spss('TIMSS/2015/T15_4_1/ASGAUSM6.sav') # read data

T15 <- TIMSS15[c('IDSTUD', 'ITSEX', 'ASDAGE', 'HOUWGT', 'JKZONE', 'JKREP', 'ASDMLOWP',
                 'ASMMAT01', 'ASMMAT02', 'ASMMAT03', 'ASMMAT04', 'ASMMAT05')] # subset columns

remove(TIMSS15) # remove unneeded file

names(T15) <- c('Student', 'Sex', 'Age', 'HWt', 'JKZ', 'JKR', 'Low',
                'PV1', 'PV2', 'PV3', 'PV4', 'PV5') # rename columns

for (i in names(T15)) {attributes(T15[[i]])$label <- NULL} # delete column labels
T15 <- data.frame(zap_labels(T15)) # remove all labels and convert to dataframe

L <- length(unique(T15$JKZ)) # number of JK zones

Size <- nrow(T15) # sample size
T15$HWt <- Size/sum(T15$HWt)*T15$HWt # equalize sample size and sum of HWt

T15_F <- T15[which(T15$Sex == 1),] # female subset
T15_M <- T15[which(T15$Sex == 2),] # male subset

FSize <- nrow(T15_F) # female sample size
MSize <- nrow(T15_M) # male sample size




##############################
##### Weighted Functions #####
##############################

# arithmetic mean
wt.mn <- function(x, w) {mean(x*(length(x)/sum(w)*w))}


# variance
wt.var <- function(x, w) {
  mn <- wt.mn(x, w)
  n <- length(x)
  wt <- n/sum(w)*w
  return(sum(wt*(x-mn)^2)/(n-sum(wt^2)/n))
}
# mean squared deviation from the mean, with the weighted variant of Bessel's correction


# quantile
wt.qnt <- function(x, w, q) {
  k <- sum(w)*q
  if (q <= .5) {
    ascend <- order(x)
    a <- x[ascend]
    b <- w[ascend]
    wt <- 0
    for (i in 1:length(a)) {
      if (wt < k) {
        wt <- wt + b[i]
      } else {
        R <- b[i-1]
        S <- a[i-2]
        V <- a[i-1]
        break
      }
    }
    return(V-(wt-k)/R*(V-S))
  }
  else {
    descend <- order(x, decreasing = T)
    a <- x[descend]
    b <- w[descend]
    wt <- sum(w)
    for (i in 1:length(a)) {
      if (wt > k) {
        wt <- wt - b[i]
      } else {
        R <- b[i-1]
        S <- a[i-1]
        V <- a[i]
        break
      }
    }
    return(V+(k-wt)/R*(S-V))
  }
}
# sort the scores and weights in ascending order of scores
# add weights until the fraction of total weight exceeds q
# perform linear interpolation between the closest scores (S and V)
# for efficiency, reverse the process if q exceeds .5


# log-transformed tail proportion ratio (LTPR)
LTPRfn <- function(d1, d2, q, v) {
  tango1 <- d1[which(d1[,v] >= max(d1[,v][which(d1[,v] < q)])),]
  tango2 <- d2[which(d2[,v] >= max(d2[,v][which(d2[,v] < q)])),]
  mango1 <- tango1[order(tango1[,v]),]
  mango2 <- tango2[order(tango2[,v]),]
  
  if (mango1[1,v] == mango1[2,v]) {
    excess <- nrow(mango1[which(mango1[,v] == mango1[1,v]),])-1
    mango1 <- mango1[-c(1:excess),]
  }
  if (mango2[1,v] == mango2[2,v]) {
    excess <- nrow(mango2[which(mango2[,v] == mango2[1,v]),])-1
    mango2 <- mango2[-c(1:excess),]
  }
  
  slice1 <- (mango1[2,v]-q)/(mango1[2,v]-mango1[1,v])*mango1[2,'HWt']
  slice2 <- (mango2[2,v]-q)/(mango2[2,v]-mango2[1,v])*mango2[2,'HWt']
  
  return(log((sum(mango1[-c(1:2),'HWt'])+slice1)/sum(d1$HWt)/
               ((sum(mango2[-c(1:2),'HWt'])+slice2)/sum(d2$HWt))))
}
# M/F ratio of weight above a threshold q divided by the subgroup weight (see U3 description)
# log-transform the ratio to produce unbiased means


# log-transformed U3 ratio (LU3R)
LU3Rfn <- function(d1, d2, q, v) {
  k <- wt.qnt(x = d2[,v], w = d2$HWt, q = q)
  tango <- d1[which(d1[,v] >= max(d1[,v][which(d1[,v] < k)])),]
  mango <- tango[order(tango[,v]),]
  if (mango[1,v] == mango[2,v]) {
    excess <- nrow(mango[which(mango[,v] == mango[1,v]),])-1
    mango <- mango[-c(1:excess),]
  }
  slice <- (mango[2,v]-k)/(mango[2,v]-mango[1,v])*mango[2,'HWt']
  return(log((sum(mango[-c(1:2),'HWt'])+slice)/sum(d1$HWt)/(1-q)))
}
# share of male weight above a female subgroup quantile (see U3 description)
# divided by the natural share of female weight above that quantile (1-q)
# log-transform the ratio to produce unbiased means


# Cohen's d
dfn <- function(d1, d2, v) {
  (wt.mn(x = d1[,v], w = d1$HWt)-wt.mn(x = d2[,v], w = d2$HWt))/
    sqrt((wt.var(x = d1[,v], w = d1$HWt)+wt.var(x = d2[,v], w = d2$HWt))/2)
}
# raw mean difference divided by quadratic mean of standard deviations


U3fn <- function(d1, d2, v) {
  q <- wt.qnt(x = d2[,v], w = d2$HWt, q = .5)
  tango <- d1[which(d1[,v] >= max(d1[,v][which(d1[,v] < q)])),]
  mango <- tango[order(tango[,v]),]
  if (mango[1,v] == mango[2,v]) {
    excess <- nrow(mango[which(mango[,v] == mango[1,v]),])-1
    mango <- mango[-c(1:excess),]
  }
  slice <- (mango[2,v]-q)/(mango[2,v]-mango[1,v])*mango[2,'HWt']
  return((sum(mango[-c(1:2),'HWt'])+slice)/sum(d1$HWt))
}
# compute female median (q)
# store the male subset with scores higher than q,
# plus the row with the score immediately below q (tango)
# sort tango in ascending order (mango)
# discard excess in the rare case that more than one student has the score immediately below q
# divide the distance between the score immediately above q and q
# by the distance between the score immediately above q and the score immediately below q
# then multiply this proportion by the weight associated with the score above q (slice)
# add slice to all weight above the score immediately above q, then divide by total male weight
# this gives the precise share of male weight above the female median


PSfn <- function(d1, d2, v) {
  D1 <- d1[order(d1[,v]),]
  
  SuperM <- numeric(FSize)
  for (i in 1:FSize) {SuperM[i] <- sum(D1[which(D1[,v] > d2[i,v]),'HWt'])}
  
  equal <- intersect(D1[,v], d2[,v])
  EWtsM <- EWtsF <- numeric(length(equal))
  for (i in 1:length(equal)) {
    EWtsM[i] <- sum(D1[which(D1[,v] %in% equal[i]),'HWt'])
    EWtsF[i] <- sum(d2[which(d2[,v] %in% equal[i]),'HWt'])
  }
  
  return(sum(c(SuperM*d2$HWt, .5*EWtsM*EWtsF))/(sum(D1$HWt)*sum(d2$HWt)))
}
# for efficiency, sort male scores and weights in order of scores (D1)
# compute the male weight above each female score (SuperM)
# store the small subset of scores that are equal to a score in the other subgroup (equal)
# find the associated weights in males and females (EWtsM, EWtsF)
# sum weights of superior males multiplied by associated female weights
# and half of the sum of equal male and female weights multiplied
# then divide by the total weight of male-female pairs
# this gives the probability that a random male has a higher score than a random female


# log-transformed variance ratio (LVR)
LVRfn <- function(d1, d2, v) {
  log(wt.var(x = d1[,v], w = d1$HWt)/wt.var(x = d2[,v], w = d2$HWt))
}
# M/F ratio of variance
# log-transform the ratio to produce unbiased means


# log-transformed tail VR (LVR_T)
LVR_Tfn <- function(d1, d2, v, t) {
  q1 <- wt.mn(x = d1[,v], w = d1$HWt)
  q2 <- wt.mn(x = d2[,v], w = d2$HWt)
  if (t == 'L') {
    MBMn <- d1[which(d1[,v] < q1),]
    FBMn <- d2[which(d2[,v] < q2),]
  } else if (t == 'R') {
    MBMn <- d1[which(d1[,v] > q1),]
    FBMn <- d2[which(d2[,v] > q2),]
  }
  return(log(sum(MBMn$HWt*(MBMn[,v]-q1)^2)/sum(MBMn$HWt)/
               (sum(FBMn$HWt*(FBMn[,v]-q2)^2)/sum(FBMn$HWt))))
}
# subset the males and females below (t = 'L') or above (t = 'R') the subgroup mean
# compute the M/F ratio of mean squared deviation from the mean in the left or right tail
# Bessel's correction not applicable because in this case the sample mean's deviation
# from the population mean causes random not systematic error
# log-transform the ratio to produce unbiased means


# log-transformed mean absolute deviation (from the median) ratio (LMADR)
LMADRfn <- function(d1, d2, v) {
  log(wt.mn(x = abs(d1[,v]-wt.qnt(x = d1[,v], w = d1$HWt, q = .5)), w = d1$HWt)/
        wt.mn(x = abs(d2[,v]-wt.qnt(x = d2[,v], w = d2$HWt, q = .5)), w = d2$HWt))
}
# M/F ratio of mean absolute deviation from the median
# log-transform for unbiased means


# log-transformed tail MADR (LMADR_T)
LMADR_Tfn <- function(d1, d2, v, t) {
  q1 <- wt.qnt(x = d1[,v], w = d1$HWt, q = .5)
  q2 <- wt.qnt(x = d2[,v], w = d2$HWt, q = .5)
  if (t == 'L') {
    MBMd <- d1[which(d1[,v] < q1),]
    FBMd <- d2[which(d2[,v] < q2),]
    return(log(wt.mn(x = q1-MBMd[,v], w = MBMd$HWt)/
                 wt.mn(x = q2-FBMd[,v], w = FBMd$HWt)))
  } else if (t == 'R') {
    MBMd <- d1[which(d1[,v] > q1),]
    FBMd <- d2[which(d2[,v] > q2),]
    return(log(wt.mn(x = MBMd[,v]-q1, w = MBMd$HWt)/
                 wt.mn(x = FBMd[,v]-q2, w = FBMd$HWt)))
  }
}
# subset the males and females below (t = 'L') or above (t = 'R') the subgroup median
# compute the M/F ratio of mean absolute deviation from the median in the left or right tail
# log-transform the ratio to produce unbiased means


# log-transformed Gini's mean difference ratio (LGMDR)
LGMDRfn <- function(d1, d2, v) {
  a1 <- sort(d1[,v])
  a2 <- sort(d2[,v])
  b1 <- d1$HWt[order(d1[,v])]/sum(d1$HWt)
  b2 <- d2$HWt[order(d2[,v])]/sum(d2$HWt)
  s1 <- cumsum(b1)
  s2 <- cumsum(b2)
  y1 <- cumsum(a1*b1)
  y2 <- cumsum(a2*b2)
  z1 <- y1/y1[MSize]
  z2 <- y2/y2[FSize]
  Gini1 <- as.numeric(t(z1[-1]) %*% s1[-MSize] - t(z1[-MSize]) %*% s1[-1])
  Gini2 <- as.numeric(t(z2[-1]) %*% s2[-FSize] - t(z2[-FSize]) %*% s2[-1])
  return(log(Gini1*wt.mn(x = d1[,v], w = d1$HWt)/(Gini2*wt.mn(x = d2[,v], w = d2$HWt))))
}
# adapted from acid::weighted.gini by Alexander Sohn
# order scores in ascending order (a1/a2)
# order weights in ascending order of scores and divide by total subgroup weight (b1/b2)
# store cumul sum of b1/b2, expressing cumul wt as a fraction of total wt (s1/s2; Lorenz curve)
# store cumul sum of ordered scores multiplied by ordered weights (y1/y2)
# and divide by the last index of y1/y2, which is the sum of a1*b1/a2*b2 (z1/z2)
# compute the Gini index (Gini1/Gini2)
# multiply by two times the mean for the mean absolute difference (but the twos cancel)
# this produces the ratio of mean abs diff between students selected randomly with replacement
# log-transform the ratio to produce unbiased means




##############################################################################
##### Means and Medians, Tail Proportion Ratios (TPRs), U3 Ratios (U3Rs) #####
##############################################################################

### Means and Medians

# total group
Mn1 <- wt.mn(x = T15$PV1, w = T15$HWt) # PV1
Mn2 <- wt.mn(x = T15$PV2, w = T15$HWt) # PV2
Mn3 <- wt.mn(x = T15$PV3, w = T15$HWt) # PV3
Mn4 <- wt.mn(x = T15$PV4, w = T15$HWt) # PV4
Mn5 <- wt.mn(x = T15$PV5, w = T15$HWt) # PV5

Mns <- c(Mn1, Mn2, Mn3, Mn4, Mn5)
Mn <- mean(Mns)

Md1 <- wt.qnt(x = T15$PV1, w = T15$HWt, q = .5) # PV1
Md2 <- wt.qnt(x = T15$PV2, w = T15$HWt, q = .5) # PV2
Md3 <- wt.qnt(x = T15$PV3, w = T15$HWt, q = .5) # PV3
Md4 <- wt.qnt(x = T15$PV4, w = T15$HWt, q = .5) # PV4
Md5 <- wt.qnt(x = T15$PV5, w = T15$HWt, q = .5) # PV5

Mds <- c(Md1, Md2, Md3, Md4, Md5)
Md <- mean(Mds)


# females
Mn1_F <- wt.mn(x = T15_F$PV1, w = T15_F$HWt) # PV1
Mn2_F <- wt.mn(x = T15_F$PV2, w = T15_F$HWt) # PV2
Mn3_F <- wt.mn(x = T15_F$PV3, w = T15_F$HWt) # PV3
Mn4_F <- wt.mn(x = T15_F$PV4, w = T15_F$HWt) # PV4
Mn5_F <- wt.mn(x = T15_F$PV5, w = T15_F$HWt) # PV5

Mns_F <- c(Mn1_F, Mn2_F, Mn3_F, Mn4_F, Mn5_F)
Mn_F <- mean(Mns_F)

Md1_F <- wt.qnt(x = T15_F$PV1, w = T15_F$HWt, q = .5) # PV1
Md2_F <- wt.qnt(x = T15_F$PV2, w = T15_F$HWt, q = .5) # PV2
Md3_F <- wt.qnt(x = T15_F$PV3, w = T15_F$HWt, q = .5) # PV3
Md4_F <- wt.qnt(x = T15_F$PV4, w = T15_F$HWt, q = .5) # PV4
Md5_F <- wt.qnt(x = T15_F$PV5, w = T15_F$HWt, q = .5) # PV5

Mds_F <- c(Md1_F, Md2_F, Md3_F, Md4_F, Md5_F)
Md_F <- mean(Mds_F)


# males
Mn1_M <- wt.mn(x = T15_M$PV1, w = T15_M$HWt) # PV1
Mn2_M <- wt.mn(x = T15_M$PV2, w = T15_M$HWt) # PV2
Mn3_M <- wt.mn(x = T15_M$PV3, w = T15_M$HWt) # PV3
Mn4_M <- wt.mn(x = T15_M$PV4, w = T15_M$HWt) # PV4
Mn5_M <- wt.mn(x = T15_M$PV5, w = T15_M$HWt) # PV5

Mns_M <- c(Mn1_M, Mn2_M, Mn3_M, Mn4_M, Mn5_M)
Mn_M <- mean(Mns_M)

Md1_M <- wt.qnt(x = T15_M$PV1, w = T15_M$HWt, q = .5) # PV1
Md2_M <- wt.qnt(x = T15_M$PV2, w = T15_M$HWt, q = .5) # PV2
Md3_M <- wt.qnt(x = T15_M$PV3, w = T15_M$HWt, q = .5) # PV3
Md4_M <- wt.qnt(x = T15_M$PV4, w = T15_M$HWt, q = .5) # PV4
Md5_M <- wt.qnt(x = T15_M$PV5, w = T15_M$HWt, q = .5) # PV5

Mds_M <- c(Md1_M, Md2_M, Md3_M, Md4_M, Md5_M)
Md_M <- mean(Mds_M)


# difference
MnDf1 <- Mn1_M-Mn1_F # PV1
MnDf2 <- Mn2_M-Mn2_F # PV2
MnDf3 <- Mn3_M-Mn3_F # PV3
MnDf4 <- Mn4_M-Mn4_F # PV4
MnDf5 <- Mn5_M-Mn5_F # PV5

MnDfs <- c(MnDf1, MnDf2, MnDf3, MnDf4, MnDf5)
MnDf <- mean(MnDfs)

MdDf1 <- Md1_M-Md1_F # PV1
MdDf2 <- Md2_M-Md2_F # PV2
MdDf3 <- Md3_M-Md3_F # PV3
MdDf4 <- Md4_M-Md4_F # PV4
MdDf5 <- Md5_M-Md5_F # PV5

MdDfs <- c(MdDf1, MdDf2, MdDf3, MdDf4, MdDf5)
MdDf <- mean(MdDfs)



### TPRs and LTPRs: mean and every 5th percentile from 5 to 95
LTPRMn_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = Mn1, v = 'PV1') # PV1 mean
LTPRMn_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = Mn2, v = 'PV2') # PV2
LTPRMn_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = Mn3, v = 'PV3') # PV3
LTPRMn_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = Mn4, v = 'PV4') # PV4
LTPRMn_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = Mn5, v = 'PV5') # PV5
LTPR05_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .05), v = 'PV1') # PV1 5th
LTPR05_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .05), v = 'PV2') # PV2
LTPR05_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .05), v = 'PV3') # PV3
LTPR05_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .05), v = 'PV4') # PV4
LTPR05_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .05), v = 'PV5') # PV5
LTPR10_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .10), v = 'PV1') # PV1 10th
LTPR10_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .10), v = 'PV2') # PV2
LTPR10_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .10), v = 'PV3') # PV3
LTPR10_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .10), v = 'PV4') # PV4
LTPR10_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .10), v = 'PV5') # PV5
LTPR15_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .15), v = 'PV1') # PV1 15th
LTPR15_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .15), v = 'PV2') # PV2
LTPR15_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .15), v = 'PV3') # PV3
LTPR15_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .15), v = 'PV4') # PV4
LTPR15_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .15), v = 'PV5') # PV5
LTPR20_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .20), v = 'PV1') # PV1 20th
LTPR20_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .20), v = 'PV2') # PV2
LTPR20_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .20), v = 'PV3') # PV3
LTPR20_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .20), v = 'PV4') # PV4
LTPR20_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .20), v = 'PV5') # PV5
LTPR25_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .25), v = 'PV1') # PV1 25th
LTPR25_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .25), v = 'PV2') # PV2
LTPR25_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .25), v = 'PV3') # PV3
LTPR25_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .25), v = 'PV4') # PV4
LTPR25_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .25), v = 'PV5') # PV5
LTPR30_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .30), v = 'PV1') # PV1 30th
LTPR30_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .30), v = 'PV2') # PV2
LTPR30_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .30), v = 'PV3') # PV3
LTPR30_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .30), v = 'PV4') # PV4
LTPR30_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .30), v = 'PV5') # PV5
LTPR35_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .35), v = 'PV1') # PV1 35th
LTPR35_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .35), v = 'PV2') # PV2
LTPR35_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .35), v = 'PV3') # PV3
LTPR35_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .35), v = 'PV4') # PV4
LTPR35_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .35), v = 'PV5') # PV5
LTPR40_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .40), v = 'PV1') # PV1 40th
LTPR40_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .40), v = 'PV2') # PV2
LTPR40_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .40), v = 'PV3') # PV3
LTPR40_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .40), v = 'PV4') # PV4
LTPR40_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .40), v = 'PV5') # PV5
LTPR45_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .45), v = 'PV1') # PV1 45th
LTPR45_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .45), v = 'PV2') # PV2
LTPR45_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .45), v = 'PV3') # PV3
LTPR45_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .45), v = 'PV4') # PV4
LTPR45_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .45), v = 'PV5') # PV5
LTPR50_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .50), v = 'PV1') # PV1 50th
LTPR50_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .50), v = 'PV2') # PV2
LTPR50_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .50), v = 'PV3') # PV3
LTPR50_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .50), v = 'PV4') # PV4
LTPR50_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .50), v = 'PV5') # PV5
LTPR55_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .55), v = 'PV1') # PV1 55th
LTPR55_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .55), v = 'PV2') # PV2
LTPR55_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .55), v = 'PV3') # PV3
LTPR55_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .55), v = 'PV4') # PV4
LTPR55_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .55), v = 'PV5') # PV5
LTPR60_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .60), v = 'PV1') # PV1 60th
LTPR60_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .60), v = 'PV2') # PV2
LTPR60_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .60), v = 'PV3') # PV3
LTPR60_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .60), v = 'PV4') # PV4
LTPR60_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .60), v = 'PV5') # PV5
LTPR65_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .65), v = 'PV1') # PV1 65th
LTPR65_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .65), v = 'PV2') # PV2
LTPR65_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .65), v = 'PV3') # PV3
LTPR65_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .65), v = 'PV4') # PV4
LTPR65_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .65), v = 'PV5') # PV5
LTPR70_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .70), v = 'PV1') # PV1 70th
LTPR70_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .70), v = 'PV2') # PV2
LTPR70_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .70), v = 'PV3') # PV3
LTPR70_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .70), v = 'PV4') # PV4
LTPR70_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .70), v = 'PV5') # PV5
LTPR75_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .75), v = 'PV1') # PV1 75th
LTPR75_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .75), v = 'PV2') # PV2
LTPR75_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .75), v = 'PV3') # PV3
LTPR75_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .75), v = 'PV4') # PV4
LTPR75_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .75), v = 'PV5') # PV5
LTPR80_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .80), v = 'PV1') # PV1 80th
LTPR80_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .80), v = 'PV2') # PV2
LTPR80_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .80), v = 'PV3') # PV3
LTPR80_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .80), v = 'PV4') # PV4
LTPR80_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .80), v = 'PV5') # PV5
LTPR85_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .85), v = 'PV1') # PV1 85th
LTPR85_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .85), v = 'PV2') # PV2
LTPR85_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .85), v = 'PV3') # PV3
LTPR85_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .85), v = 'PV4') # PV4
LTPR85_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .85), v = 'PV5') # PV5
LTPR90_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .90), v = 'PV1') # PV1 90th
LTPR90_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .90), v = 'PV2') # PV2
LTPR90_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .90), v = 'PV3') # PV3
LTPR90_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .90), v = 'PV4') # PV4
LTPR90_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .90), v = 'PV5') # PV5
LTPR95_1 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV1, T15$HWt, .95), v = 'PV1') # PV1 95th
LTPR95_2 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV2, T15$HWt, .95), v = 'PV2') # PV2
LTPR95_3 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV3, T15$HWt, .95), v = 'PV3') # PV3
LTPR95_4 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV4, T15$HWt, .95), v = 'PV4') # PV4
LTPR95_5 <- LTPRfn(d1 = T15_M, d2 = T15_F, q = wt.qnt(T15$PV5, T15$HWt, .95), v = 'PV5') # PV5

LTPRMns <- c(LTPRMn_1, LTPRMn_2, LTPRMn_3, LTPRMn_4, LTPRMn_5) # mean
LTPR05s <- c(LTPR05_1, LTPR05_2, LTPR05_3, LTPR05_4, LTPR05_5) # 5th
LTPR10s <- c(LTPR10_1, LTPR10_2, LTPR10_3, LTPR10_4, LTPR10_5) # 10th
LTPR15s <- c(LTPR15_1, LTPR15_2, LTPR15_3, LTPR15_4, LTPR15_5) # 15th
LTPR20s <- c(LTPR20_1, LTPR20_2, LTPR20_3, LTPR20_4, LTPR20_5) # 20th
LTPR25s <- c(LTPR25_1, LTPR25_2, LTPR25_3, LTPR25_4, LTPR25_5) # 25th
LTPR30s <- c(LTPR30_1, LTPR30_2, LTPR30_3, LTPR30_4, LTPR30_5) # 30th
LTPR35s <- c(LTPR35_1, LTPR35_2, LTPR35_3, LTPR35_4, LTPR35_5) # 35th
LTPR40s <- c(LTPR40_1, LTPR40_2, LTPR40_3, LTPR40_4, LTPR40_5) # 40th
LTPR45s <- c(LTPR45_1, LTPR45_2, LTPR45_3, LTPR45_4, LTPR45_5) # 45th
LTPR50s <- c(LTPR50_1, LTPR50_2, LTPR50_3, LTPR50_4, LTPR50_5) # 50th
LTPR55s <- c(LTPR55_1, LTPR55_2, LTPR55_3, LTPR55_4, LTPR55_5) # 55th
LTPR60s <- c(LTPR60_1, LTPR60_2, LTPR60_3, LTPR60_4, LTPR60_5) # 60th
LTPR65s <- c(LTPR65_1, LTPR65_2, LTPR65_3, LTPR65_4, LTPR65_5) # 65th
LTPR70s <- c(LTPR70_1, LTPR70_2, LTPR70_3, LTPR70_4, LTPR70_5) # 70th
LTPR75s <- c(LTPR75_1, LTPR75_2, LTPR75_3, LTPR75_4, LTPR75_5) # 75th
LTPR80s <- c(LTPR80_1, LTPR80_2, LTPR80_3, LTPR80_4, LTPR80_5) # 80th
LTPR85s <- c(LTPR85_1, LTPR85_2, LTPR85_3, LTPR85_4, LTPR85_5) # 85th
LTPR90s <- c(LTPR90_1, LTPR90_2, LTPR90_3, LTPR90_4, LTPR90_5) # 90th
LTPR95s <- c(LTPR95_1, LTPR95_2, LTPR95_3, LTPR95_4, LTPR95_5) # 95th

LTPRMn <- mean(LTPRMns) # mean
LTPR05 <- mean(LTPR05s) # 5th
LTPR10 <- mean(LTPR10s) # 10th
LTPR15 <- mean(LTPR15s) # 15th
LTPR20 <- mean(LTPR20s) # 20th
LTPR25 <- mean(LTPR25s) # 25th
LTPR30 <- mean(LTPR30s) # 30th
LTPR35 <- mean(LTPR35s) # 35th
LTPR40 <- mean(LTPR40s) # 40th
LTPR45 <- mean(LTPR45s) # 45th
LTPR50 <- mean(LTPR50s) # 50th
LTPR55 <- mean(LTPR55s) # 55th
LTPR60 <- mean(LTPR60s) # 60th
LTPR65 <- mean(LTPR65s) # 65th
LTPR70 <- mean(LTPR70s) # 70th
LTPR75 <- mean(LTPR75s) # 75th
LTPR80 <- mean(LTPR80s) # 80th
LTPR85 <- mean(LTPR85s) # 85th
LTPR90 <- mean(LTPR90s) # 90th
LTPR95 <- mean(LTPR95s) # 95th

TPRMn <- exp(LTPRMn) # mean
TPR05 <- exp(LTPR05) # 5th
TPR10 <- exp(LTPR10) # 10th
TPR15 <- exp(LTPR15) # 15th
TPR20 <- exp(LTPR20) # 20th
TPR25 <- exp(LTPR25) # 25th
TPR30 <- exp(LTPR30) # 30th
TPR35 <- exp(LTPR35) # 35th
TPR40 <- exp(LTPR40) # 40th
TPR45 <- exp(LTPR45) # 45th
TPR50 <- exp(LTPR50) # 50th
TPR55 <- exp(LTPR55) # 55th
TPR60 <- exp(LTPR60) # 60th
TPR65 <- exp(LTPR65) # 65th
TPR70 <- exp(LTPR70) # 70th
TPR75 <- exp(LTPR75) # 75th
TPR80 <- exp(LTPR80) # 80th
TPR85 <- exp(LTPR85) # 85th
TPR90 <- exp(LTPR90) # 90th
TPR95 <- exp(LTPR95) # 95th



### LTPR tail-center differences

# LTPR difference between 95th percentile and median
Md95T_1 <- LTPR95_1-LTPR50_1 # PV1
Md95T_2 <- LTPR95_2-LTPR50_2 # PV2
Md95T_3 <- LTPR95_3-LTPR50_3 # PV3
Md95T_4 <- LTPR95_4-LTPR50_4 # PV4
Md95T_5 <- LTPR95_5-LTPR50_5 # PV5

Md95Ts <- c(Md95T_1, Md95T_2, Md95T_3, Md95T_4, Md95T_5)
Md95T <- mean(Md95Ts)

# LTPR difference between 90th percentile and median
Md90T_1 <- LTPR90_1-LTPR50_1 # PV1
Md90T_2 <- LTPR90_2-LTPR50_2 # PV2
Md90T_3 <- LTPR90_3-LTPR50_3 # PV3
Md90T_4 <- LTPR90_4-LTPR50_4 # PV4
Md90T_5 <- LTPR90_5-LTPR50_5 # PV5

Md90Ts <- c(Md90T_1, Md90T_2, Md90T_3, Md90T_4, Md90T_5)
Md90T <- mean(Md90Ts)

# LTPR difference between median and 10th percentile
Md10T_1 <- LTPR50_1-LTPR10_1 # PV1
Md10T_2 <- LTPR50_2-LTPR10_2 # PV2
Md10T_3 <- LTPR50_3-LTPR10_3 # PV3
Md10T_4 <- LTPR50_4-LTPR10_4 # PV4
Md10T_5 <- LTPR50_5-LTPR10_5 # PV5

Md10Ts <- c(Md10T_1, Md10T_2, Md10T_3, Md10T_4, Md10T_5)
Md10T <- mean(Md10Ts)

# LTPR difference between median and 5th percentile
Md05T_1 <- LTPR50_1-LTPR05_1 # PV1
Md05T_2 <- LTPR50_2-LTPR05_2 # PV2
Md05T_3 <- LTPR50_3-LTPR05_3 # PV3
Md05T_4 <- LTPR50_4-LTPR05_4 # PV4
Md05T_5 <- LTPR50_5-LTPR05_5 # PV5

Md05Ts <- c(Md05T_1, Md05T_2, Md05T_3, Md05T_4, Md05T_5)
Md05T <- mean(Md05Ts)


# LTPR difference between 95th percentile and mean
Mn95T_1 <- LTPR95_1-LTPRMn_1 # PV1
Mn95T_2 <- LTPR95_2-LTPRMn_2 # PV2
Mn95T_3 <- LTPR95_3-LTPRMn_3 # PV3
Mn95T_4 <- LTPR95_4-LTPRMn_4 # PV4
Mn95T_5 <- LTPR95_5-LTPRMn_5 # PV5

Mn95Ts <- c(Mn95T_1, Mn95T_2, Mn95T_3, Mn95T_4, Mn95T_5)
Mn95T <- mean(Mn95Ts)

# LTPR difference between 90th percentile and mean
Mn90T_1 <- LTPR90_1-LTPRMn_1 # PV1
Mn90T_2 <- LTPR90_2-LTPRMn_2 # PV2
Mn90T_3 <- LTPR90_3-LTPRMn_3 # PV3
Mn90T_4 <- LTPR90_4-LTPRMn_4 # PV4
Mn90T_5 <- LTPR90_5-LTPRMn_5 # PV5

Mn90Ts <- c(Mn90T_1, Mn90T_2, Mn90T_3, Mn90T_4, Mn90T_5)
Mn90T <- mean(Mn90Ts)

# LTPR difference between mean and 10th percentile
Mn10T_1 <- LTPRMn_1-LTPR10_1 # PV1
Mn10T_2 <- LTPRMn_2-LTPR10_2 # PV2
Mn10T_3 <- LTPRMn_3-LTPR10_3 # PV3
Mn10T_4 <- LTPRMn_4-LTPR10_4 # PV4
Mn10T_5 <- LTPRMn_5-LTPR10_5 # PV5

Mn10Ts <- c(Mn10T_1, Mn10T_2, Mn10T_3, Mn10T_4, Mn10T_5)
Mn10T <- mean(Mn10Ts)

# LTPR difference between mean and 5th percentile
Mn05T_1 <- LTPRMn_1-LTPR05_1 # PV1
Mn05T_2 <- LTPRMn_2-LTPR05_2 # PV2
Mn05T_3 <- LTPRMn_3-LTPR05_3 # PV3
Mn05T_4 <- LTPRMn_4-LTPR05_4 # PV4
Mn05T_5 <- LTPRMn_5-LTPR05_5 # PV5

Mn05Ts <- c(Mn05T_1, Mn05T_2, Mn05T_3, Mn05T_4, Mn05T_5)
Mn05T <- mean(Mn05Ts)



### U3Rs and LU3Rs: every 5th percentile from 5 to 95
LU3R05_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .05, v = 'PV1') # PV1 5th
LU3R05_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .05, v = 'PV2') # PV2
LU3R05_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .05, v = 'PV3') # PV3
LU3R05_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .05, v = 'PV4') # PV4
LU3R05_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .05, v = 'PV5') # PV5
LU3R10_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .10, v = 'PV1') # PV1 10th
LU3R10_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .10, v = 'PV2') # PV2
LU3R10_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .10, v = 'PV3') # PV3
LU3R10_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .10, v = 'PV4') # PV4
LU3R10_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .10, v = 'PV5') # PV5
LU3R15_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .15, v = 'PV1') # PV1 15th
LU3R15_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .15, v = 'PV2') # PV2
LU3R15_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .15, v = 'PV3') # PV3
LU3R15_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .15, v = 'PV4') # PV4
LU3R15_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .15, v = 'PV5') # PV5
LU3R20_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .20, v = 'PV1') # PV1 20th
LU3R20_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .20, v = 'PV2') # PV2
LU3R20_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .20, v = 'PV3') # PV3
LU3R20_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .20, v = 'PV4') # PV4
LU3R20_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .20, v = 'PV5') # PV5
LU3R25_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .25, v = 'PV1') # PV1 25th
LU3R25_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .25, v = 'PV2') # PV2
LU3R25_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .25, v = 'PV3') # PV3
LU3R25_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .25, v = 'PV4') # PV4
LU3R25_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .25, v = 'PV5') # PV5
LU3R30_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .30, v = 'PV1') # PV1 30th
LU3R30_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .30, v = 'PV2') # PV2
LU3R30_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .30, v = 'PV3') # PV3
LU3R30_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .30, v = 'PV4') # PV4
LU3R30_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .30, v = 'PV5') # PV5
LU3R35_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .35, v = 'PV1') # PV1 35th
LU3R35_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .35, v = 'PV2') # PV2
LU3R35_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .35, v = 'PV3') # PV3
LU3R35_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .35, v = 'PV4') # PV4
LU3R35_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .35, v = 'PV5') # PV5
LU3R40_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .40, v = 'PV1') # PV1 40th
LU3R40_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .40, v = 'PV2') # PV2
LU3R40_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .40, v = 'PV3') # PV3
LU3R40_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .40, v = 'PV4') # PV4
LU3R40_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .40, v = 'PV5') # PV5
LU3R45_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .45, v = 'PV1') # PV1 45th
LU3R45_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .45, v = 'PV2') # PV2
LU3R45_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .45, v = 'PV3') # PV3
LU3R45_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .45, v = 'PV4') # PV4
LU3R45_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .45, v = 'PV5') # PV5
LU3R50_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .50, v = 'PV1') # PV1 50th
LU3R50_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .50, v = 'PV2') # PV2
LU3R50_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .50, v = 'PV3') # PV3
LU3R50_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .50, v = 'PV4') # PV4
LU3R50_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .50, v = 'PV5') # PV5
LU3R55_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .55, v = 'PV1') # PV1 55th
LU3R55_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .55, v = 'PV2') # PV2
LU3R55_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .55, v = 'PV3') # PV3
LU3R55_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .55, v = 'PV4') # PV4
LU3R55_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .55, v = 'PV5') # PV5
LU3R60_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .60, v = 'PV1') # PV1 60th
LU3R60_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .60, v = 'PV2') # PV2
LU3R60_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .60, v = 'PV3') # PV3
LU3R60_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .60, v = 'PV4') # PV4
LU3R60_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .60, v = 'PV5') # PV5
LU3R65_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .65, v = 'PV1') # PV1 65th
LU3R65_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .65, v = 'PV2') # PV2
LU3R65_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .65, v = 'PV3') # PV3
LU3R65_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .65, v = 'PV4') # PV4
LU3R65_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .65, v = 'PV5') # PV5
LU3R70_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .70, v = 'PV1') # PV1 70th
LU3R70_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .70, v = 'PV2') # PV2
LU3R70_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .70, v = 'PV3') # PV3
LU3R70_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .70, v = 'PV4') # PV4
LU3R70_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .70, v = 'PV5') # PV5
LU3R75_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .75, v = 'PV1') # PV1 75th
LU3R75_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .75, v = 'PV2') # PV2
LU3R75_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .75, v = 'PV3') # PV3
LU3R75_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .75, v = 'PV4') # PV4
LU3R75_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .75, v = 'PV5') # PV5
LU3R80_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .80, v = 'PV1') # PV1 80th
LU3R80_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .80, v = 'PV2') # PV2
LU3R80_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .80, v = 'PV3') # PV3
LU3R80_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .80, v = 'PV4') # PV4
LU3R80_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .80, v = 'PV5') # PV5
LU3R85_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .85, v = 'PV1') # PV1 85th
LU3R85_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .85, v = 'PV2') # PV2
LU3R85_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .85, v = 'PV3') # PV3
LU3R85_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .85, v = 'PV4') # PV4
LU3R85_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .85, v = 'PV5') # PV5
LU3R90_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .90, v = 'PV1') # PV1 90th
LU3R90_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .90, v = 'PV2') # PV2
LU3R90_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .90, v = 'PV3') # PV3
LU3R90_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .90, v = 'PV4') # PV4
LU3R90_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .90, v = 'PV5') # PV5
LU3R95_1 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .95, v = 'PV1') # PV1 95th
LU3R95_2 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .95, v = 'PV2') # PV2
LU3R95_3 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .95, v = 'PV3') # PV3
LU3R95_4 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .95, v = 'PV4') # PV4
LU3R95_5 <- LU3Rfn(d1 = T15_M, d2 = T15_F, q = .95, v = 'PV5') # PV5

LU3R05s <- c(LU3R05_1, LU3R05_2, LU3R05_3, LU3R05_4, LU3R05_5) # 5th
LU3R10s <- c(LU3R10_1, LU3R10_2, LU3R10_3, LU3R10_4, LU3R10_5) # 10th
LU3R15s <- c(LU3R15_1, LU3R15_2, LU3R15_3, LU3R15_4, LU3R15_5) # 15th
LU3R20s <- c(LU3R20_1, LU3R20_2, LU3R20_3, LU3R20_4, LU3R20_5) # 20th
LU3R25s <- c(LU3R25_1, LU3R25_2, LU3R25_3, LU3R25_4, LU3R25_5) # 25th
LU3R30s <- c(LU3R30_1, LU3R30_2, LU3R30_3, LU3R30_4, LU3R30_5) # 30th
LU3R35s <- c(LU3R35_1, LU3R35_2, LU3R35_3, LU3R35_4, LU3R35_5) # 35th
LU3R40s <- c(LU3R40_1, LU3R40_2, LU3R40_3, LU3R40_4, LU3R40_5) # 40th
LU3R45s <- c(LU3R45_1, LU3R45_2, LU3R45_3, LU3R45_4, LU3R45_5) # 45th
LU3R50s <- c(LU3R50_1, LU3R50_2, LU3R50_3, LU3R50_4, LU3R50_5) # 50th
LU3R55s <- c(LU3R55_1, LU3R55_2, LU3R55_3, LU3R55_4, LU3R55_5) # 55th
LU3R60s <- c(LU3R60_1, LU3R60_2, LU3R60_3, LU3R60_4, LU3R60_5) # 60th
LU3R65s <- c(LU3R65_1, LU3R65_2, LU3R65_3, LU3R65_4, LU3R65_5) # 65th
LU3R70s <- c(LU3R70_1, LU3R70_2, LU3R70_3, LU3R70_4, LU3R70_5) # 70th
LU3R75s <- c(LU3R75_1, LU3R75_2, LU3R75_3, LU3R75_4, LU3R75_5) # 75th
LU3R80s <- c(LU3R80_1, LU3R80_2, LU3R80_3, LU3R80_4, LU3R80_5) # 80th
LU3R85s <- c(LU3R85_1, LU3R85_2, LU3R85_3, LU3R85_4, LU3R85_5) # 85th
LU3R90s <- c(LU3R90_1, LU3R90_2, LU3R90_3, LU3R90_4, LU3R90_5) # 90th
LU3R95s <- c(LU3R95_1, LU3R95_2, LU3R95_3, LU3R95_4, LU3R95_5) # 95th

LU3R05 <- mean(LU3R05s) # 5th
LU3R10 <- mean(LU3R10s) # 10th
LU3R15 <- mean(LU3R15s) # 15th
LU3R20 <- mean(LU3R20s) # 20th
LU3R25 <- mean(LU3R25s) # 25th
LU3R30 <- mean(LU3R30s) # 30th
LU3R35 <- mean(LU3R35s) # 35th
LU3R40 <- mean(LU3R40s) # 40th
LU3R45 <- mean(LU3R45s) # 45th
LU3R50 <- mean(LU3R50s) # 50th
LU3R55 <- mean(LU3R55s) # 55th
LU3R60 <- mean(LU3R60s) # 60th
LU3R65 <- mean(LU3R65s) # 65th
LU3R70 <- mean(LU3R70s) # 70th
LU3R75 <- mean(LU3R75s) # 75th
LU3R80 <- mean(LU3R80s) # 80th
LU3R85 <- mean(LU3R85s) # 85th
LU3R90 <- mean(LU3R90s) # 90th
LU3R95 <- mean(LU3R95s) # 95th

U3R05 <- exp(LU3R05) # 5th
U3R10 <- exp(LU3R10) # 10th
U3R15 <- exp(LU3R15) # 15th
U3R20 <- exp(LU3R20) # 20th
U3R25 <- exp(LU3R25) # 25th
U3R30 <- exp(LU3R30) # 30th
U3R35 <- exp(LU3R35) # 35th
U3R40 <- exp(LU3R40) # 40th
U3R45 <- exp(LU3R45) # 45th
U3R50 <- exp(LU3R50) # 50th
U3R55 <- exp(LU3R55) # 55th
U3R60 <- exp(LU3R60) # 60th
U3R65 <- exp(LU3R65) # 65th
U3R70 <- exp(LU3R70) # 70th
U3R75 <- exp(LU3R75) # 75th
U3R80 <- exp(LU3R80) # 80th
U3R85 <- exp(LU3R85) # 85th
U3R90 <- exp(LU3R90) # 90th
U3R95 <- exp(LU3R95) # 95th



### U3R tail-center differences

# LU3R difference between 95th percentile and median
Md95U_1 <- LU3R95_1-LU3R50_1 # PV1
Md95U_2 <- LU3R95_2-LU3R50_2 # PV2
Md95U_3 <- LU3R95_3-LU3R50_3 # PV3
Md95U_4 <- LU3R95_4-LU3R50_4 # PV4
Md95U_5 <- LU3R95_5-LU3R50_5 # PV5

Md95Us <- c(Md95U_1, Md95U_2, Md95U_3, Md95U_4, Md95U_5)
Md95U <- mean(Md95Us)

# LU3R difference between 90th percentile and median
Md90U_1 <- LU3R90_1-LU3R50_1 # PV1
Md90U_2 <- LU3R90_2-LU3R50_2 # PV2
Md90U_3 <- LU3R90_3-LU3R50_3 # PV3
Md90U_4 <- LU3R90_4-LU3R50_4 # PV4
Md90U_5 <- LU3R90_5-LU3R50_5 # PV5

Md90Us <- c(Md90U_1, Md90U_2, Md90U_3, Md90U_4, Md90U_5)
Md90U <- mean(Md90Us)

# LU3R difference between median and 10th percentile
Md10U_1 <- LU3R50_1-LU3R10_1 # PV1
Md10U_2 <- LU3R50_2-LU3R10_2 # PV2
Md10U_3 <- LU3R50_3-LU3R10_3 # PV3
Md10U_4 <- LU3R50_4-LU3R10_4 # PV4
Md10U_5 <- LU3R50_5-LU3R10_5 # PV5

Md10Us <- c(Md10U_1, Md10U_2, Md10U_3, Md10U_4, Md10U_5)
Md10U <- mean(Md10Us)

# LU3R difference between median and 5th percentile
Md05U_1 <- LU3R50_1-LU3R05_1 # PV1
Md05U_2 <- LU3R50_2-LU3R05_2 # PV2
Md05U_3 <- LU3R50_3-LU3R05_3 # PV3
Md05U_4 <- LU3R50_4-LU3R05_4 # PV4
Md05U_5 <- LU3R50_5-LU3R05_5 # PV5

Md05Us <- c(Md05U_1, Md05U_2, Md05U_3, Md05U_4, Md05U_5)
Md05U <- mean(Md05Us)




##############################
##### Other Effect Sizes #####
##############################

### Cohen's d

d1 <- dfn(d1 = T15_M, d2 = T15_F, v = 'PV1') # PV1
d2 <- dfn(d1 = T15_M, d2 = T15_F, v = 'PV2') # PV2
d3 <- dfn(d1 = T15_M, d2 = T15_F, v = 'PV3') # PV3
d4 <- dfn(d1 = T15_M, d2 = T15_F, v = 'PV4') # PV4
d5 <- dfn(d1 = T15_M, d2 = T15_F, v = 'PV5') # PV5

ds <- c(d1, d2, d3, d4, d5)
d <- mean(ds)


### U3

U3_1 <- U3fn(d1 = T15_M, d2 = T15_F, v = 'PV1') # PV1
U3_2 <- U3fn(d1 = T15_M, d2 = T15_F, v = 'PV2') # PV2
U3_3 <- U3fn(d1 = T15_M, d2 = T15_F, v = 'PV3') # PV3
U3_4 <- U3fn(d1 = T15_M, d2 = T15_F, v = 'PV4') # PV4
U3_5 <- U3fn(d1 = T15_M, d2 = T15_F, v = 'PV5') # PV5

U3s <- c(U3_1, U3_2, U3_3, U3_4, U3_5)
U3 <- mean(U3s)


### Probability of superiority (PS)

PS1 <- PSfn(d1 = T15_M, d2 = T15_F, v = 'PV1') # PV1
PS2 <- PSfn(d1 = T15_M, d2 = T15_F, v = 'PV2') # PV2
PS3 <- PSfn(d1 = T15_M, d2 = T15_F, v = 'PV3') # PV3
PS4 <- PSfn(d1 = T15_M, d2 = T15_F, v = 'PV4') # PV4
PS5 <- PSfn(d1 = T15_M, d2 = T15_F, v = 'PV5') # PV5

PSs <- c(PS1, PS2, PS3, PS4, PS5)
PS <- mean(PSs)


### Variance ratio (VR), Log-transformed VR (LVR)

LVR1 <- LVRfn(d1 = T15_M, d2 = T15_F, v = 'PV1') # PV1
LVR2 <- LVRfn(d1 = T15_M, d2 = T15_F, v = 'PV2') # PV2
LVR3 <- LVRfn(d1 = T15_M, d2 = T15_F, v = 'PV3') # PV3
LVR4 <- LVRfn(d1 = T15_M, d2 = T15_F, v = 'PV4') # PV4
LVR5 <- LVRfn(d1 = T15_M, d2 = T15_F, v = 'PV5') # PV5

LVRs <- c(LVR1, LVR2, LVR3, LVR4, LVR5)
LVR <- mean(LVRs)
VR <- exp(LVR)


### Left and right VR and LVR (VR_L, VR_R, LVR_L, LVR_R)

LVR_L1 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV1', t = 'L') # PV1
LVR_L2 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV2', t = 'L') # PV2
LVR_L3 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV3', t = 'L') # PV3
LVR_L4 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV4', t = 'L') # PV4
LVR_L5 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV5', t = 'L') # PV5

LVR_Ls <- c(LVR_L1, LVR_L2, LVR_L3, LVR_L4, LVR_L5)
LVR_L <- mean(LVR_Ls)
VR_L <- exp(LVR_L)

LVR_R1 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV1', t = 'R') # PV1
LVR_R2 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV2', t = 'R') # PV2
LVR_R3 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV3', t = 'R') # PV3
LVR_R4 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV4', t = 'R') # PV4
LVR_R5 <- LVR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV5', t = 'R') # PV5

LVR_Rs <- c(LVR_R1, LVR_R2, LVR_R3, LVR_R4, LVR_R5)
LVR_R <- mean(LVR_Rs)
VR_R <- exp(LVR_R)


### Mean absolute deviation (from the median) ratio (MADR), Log-transformed MADR (LMADR)

LMADR1 <- LMADRfn(d1 = T15_M, d2 = T15_F, v = 'PV1') # PV1
LMADR2 <- LMADRfn(d1 = T15_M, d2 = T15_F, v = 'PV2') # PV2
LMADR3 <- LMADRfn(d1 = T15_M, d2 = T15_F, v = 'PV3') # PV3
LMADR4 <- LMADRfn(d1 = T15_M, d2 = T15_F, v = 'PV4') # PV4
LMADR5 <- LMADRfn(d1 = T15_M, d2 = T15_F, v = 'PV5') # PV5

LMADRs <- c(LMADR1, LMADR2, LMADR3, LMADR4, LMADR5)
LMADR <- mean(LMADRs)
MADR <- exp(LMADR)


### Left and right MADR and LMADR (MADR_L, MADR_R, LMADR_L, LMADR_R)

LMADR_L1 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV1', t = 'L') # PV1
LMADR_L2 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV2', t = 'L') # PV2
LMADR_L3 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV3', t = 'L') # PV3
LMADR_L4 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV4', t = 'L') # PV4
LMADR_L5 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV5', t = 'L') # PV5

LMADR_Ls <- c(LMADR_L1, LMADR_L2, LMADR_L3, LMADR_L4, LMADR_L5)
LMADR_L <- mean(LMADR_Ls)
MADR_L <- exp(LMADR_L)

LMADR_R1 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV1', t = 'R') # PV1
LMADR_R2 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV2', t = 'R') # PV2
LMADR_R3 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV3', t = 'R') # PV3
LMADR_R4 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV4', t = 'R') # PV4
LMADR_R5 <- LMADR_Tfn(d1 = T15_M, d2 = T15_F, v = 'PV5', t = 'R') # PV5

LMADR_Rs <- c(LMADR_R1, LMADR_R2, LMADR_R3, LMADR_R4, LMADR_R5)
LMADR_R <- mean(LMADR_Rs)
MADR_R <- exp(LMADR_R)


### Gini's mean difference ratio (GMDR), Log-transformed GMDR (LGMDR)

LGMDR1 <- LGMDRfn(d1 = T15_M, d2 = T15_F, v = 'PV1') # PV1
LGMDR2 <- LGMDRfn(d1 = T15_M, d2 = T15_F, v = 'PV2') # PV2
LGMDR3 <- LGMDRfn(d1 = T15_M, d2 = T15_F, v = 'PV3') # PV3
LGMDR4 <- LGMDRfn(d1 = T15_M, d2 = T15_F, v = 'PV4') # PV4
LGMDR5 <- LGMDRfn(d1 = T15_M, d2 = T15_F, v = 'PV5') # PV5

LGMDRs <- c(LGMDR1, LGMDR2, LGMDR3, LGMDR4, LGMDR5)
LGMDR <- mean(LGMDRs)
GMDR <- exp(LGMDR)




###############################################
##### Other Effect Sizes Adjusted for Age #####
###############################################

A15 <- T15[!is.na(T15$Age),] # T15 with Age NAs (if any) excluded

A15_F <- A15[which(A15$Sex == 1),] # female subset
A15_M <- A15[which(A15$Sex == 2),] # male subset

# sex differences in age distribution
AgeU3 <- U3fn(d1 = A15_M, d2 = A15_F, v = 'Age')
AgeLMADR <- LMADRfn(d1 = A15_M, d2 = A15_F, v = 'Age')
AgeMADR <- exp(AgeLMADR)

# control for the correlation between age and score
Slope1 <- lm(formula = PV1 ~ Age, data = A15, weights = HWt)$coefficients[2] # PV1
Slope2 <- lm(formula = PV2 ~ Age, data = A15, weights = HWt)$coefficients[2] # PV2
Slope3 <- lm(formula = PV3 ~ Age, data = A15, weights = HWt)$coefficients[2] # PV3
Slope4 <- lm(formula = PV4 ~ Age, data = A15, weights = HWt)$coefficients[2] # PV4
Slope5 <- lm(formula = PV5 ~ Age, data = A15, weights = HWt)$coefficients[2] # PV5

AgeMn <- wt.mn(x = A15$Age, w = A15$HWt) # mean age

# if there are Age NAs, replace with mean age so their nominally age-corrected scores do not change
if (length(T15[is.na(T15$Age),'Age']) != 0) {T15[is.na(T15$Age),'Age'] <- AgeMn}

# new scores linearly corrected for age
PV1A <- T15$PV1+(AgeMn-T15$Age)*Slope1 # PV1
PV2A <- T15$PV2+(AgeMn-T15$Age)*Slope2 # PV2
PV3A <- T15$PV3+(AgeMn-T15$Age)*Slope3 # PV3
PV4A <- T15$PV4+(AgeMn-T15$Age)*Slope4 # PV4
PV5A <- T15$PV5+(AgeMn-T15$Age)*Slope5 # PV5

A15 <- data.frame(T15, PV1A, PV2A, PV3A, PV4A, PV5A) # T15 with age-corrected scores
A15_F <- A15[which(A15$Sex == 1),] # female subset
A15_M <- A15[which(A15$Sex == 2),] # male subset


### Age-corrected Cohen's d (d_A)

d_A <- mean(c(dfn(d1 = A15_M, d2 = A15_F, v = 'PV1A'),  # PV1
              dfn(d1 = A15_M, d2 = A15_F, v = 'PV2A'),  # PV2
              dfn(d1 = A15_M, d2 = A15_F, v = 'PV3A'),  # PV3
              dfn(d1 = A15_M, d2 = A15_F, v = 'PV4A'),  # PV4
              dfn(d1 = A15_M, d2 = A15_F, v = 'PV5A'))) # PV5


### Age-corrected U3 (U3_A)

U3_A <- mean(c(U3fn(d1 = A15_M, d2 = A15_F, v = 'PV1A'),  # PV1
               U3fn(d1 = A15_M, d2 = A15_F, v = 'PV2A'),  # PV2
               U3fn(d1 = A15_M, d2 = A15_F, v = 'PV3A'),  # PV3
               U3fn(d1 = A15_M, d2 = A15_F, v = 'PV4A'),  # PV4
               U3fn(d1 = A15_M, d2 = A15_F, v = 'PV5A'))) # PV5


### Age-corrected Probability of superiority (PS_A)

PS_A <- mean(c(PSfn(d1 = A15_M, d2 = A15_F, v = 'PV1A'),  # PV1
               PSfn(d1 = A15_M, d2 = A15_F, v = 'PV2A'),  # PV2
               PSfn(d1 = A15_M, d2 = A15_F, v = 'PV3A'),  # PV3
               PSfn(d1 = A15_M, d2 = A15_F, v = 'PV4A'),  # PV4
               PSfn(d1 = A15_M, d2 = A15_F, v = 'PV5A'))) # PV5


### Age-corrected Variance ratio (VR_A), Log-transformed VR (LVR_A)

LVR_A <- mean(c(LVRfn(d1 = A15_M, d2 = A15_F, v = 'PV1A'),  # PV1
                LVRfn(d1 = A15_M, d2 = A15_F, v = 'PV2A'),  # PV2
                LVRfn(d1 = A15_M, d2 = A15_F, v = 'PV3A'),  # PV3
                LVRfn(d1 = A15_M, d2 = A15_F, v = 'PV4A'),  # PV4
                LVRfn(d1 = A15_M, d2 = A15_F, v = 'PV5A'))) # PV5

VR_A <- exp(LVR_A)


### Age-corrected Left and right VR and LVR (VR_LA, VR_RA, LVR_LA, LVR_RA)

LVR_LA <- mean(c(LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV1A', t = 'L'),  # PV1
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV2A', t = 'L'),  # PV2
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV3A', t = 'L'),  # PV3
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV4A', t = 'L'),  # PV4
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV5A', t = 'L'))) # PV5

VR_LA <- exp(LVR_LA)

LVR_RA <- mean(c(LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV1A', t = 'R'),  # PV1
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV2A', t = 'R'),  # PV2
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV3A', t = 'R'),  # PV3
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV4A', t = 'R'),  # PV4
                 LVR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV5A', t = 'R'))) # PV5

VR_RA <- exp(LVR_RA)


### Age-corrected Mean absolute deviation ratio (MADR_A), Log-transformed MADR (LMADR_A)

LMADR_A <- mean(c(LMADRfn(d1 = A15_M, d2 = A15_F, v = 'PV1A'),  # PV1
                  LMADRfn(d1 = A15_M, d2 = A15_F, v = 'PV2A'),  # PV2
                  LMADRfn(d1 = A15_M, d2 = A15_F, v = 'PV3A'),  # PV3
                  LMADRfn(d1 = A15_M, d2 = A15_F, v = 'PV4A'),  # PV4
                  LMADRfn(d1 = A15_M, d2 = A15_F, v = 'PV5A'))) # PV5

MADR_A <- exp(LMADR_A)


### Age-corrected Left and right MADR and LMADR (MADR_LA, MADR_RA, LMADR_LA, LMADR_RA)

LMADR_LA <- mean(c(LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV1A', t = 'L'),  # PV1
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV2A', t = 'L'),  # PV2
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV3A', t = 'L'),  # PV3
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV4A', t = 'L'),  # PV4
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV5A', t = 'L'))) # PV5

MADR_LA <- exp(LMADR_LA)

LMADR_RA <- mean(c(LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV1A', t = 'R'),  # PV1
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV2A', t = 'R'),  # PV2
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV3A', t = 'R'),  # PV3
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV4A', t = 'R'),  # PV4
                   LMADR_Tfn(d1 = A15_M, d2 = A15_F, v = 'PV5A', t = 'R'))) # PV5

MADR_RA <- exp(LMADR_RA)


### Age-corrected Gini's mean difference ratio (GMDR), Log-transformed GMDR (LGMDR)

LGMDR_A <- mean(c(LGMDRfn(d1 = A15_M, d2 = A15_F, v = 'PV1A'),  # PV1
                  LGMDRfn(d1 = A15_M, d2 = A15_F, v = 'PV2A'),  # PV2
                  LGMDRfn(d1 = A15_M, d2 = A15_F, v = 'PV3A'),  # PV3
                  LGMDRfn(d1 = A15_M, d2 = A15_F, v = 'PV4A'),  # PV4
                  LGMDRfn(d1 = A15_M, d2 = A15_F, v = 'PV5A'))) # PV5

GMDR_A <- exp(LGMDR_A)




###########################
##### Standard Errors #####
###########################

#### SEs: Means and Medians ####

J1 <- J2 <- J3 <- J4 <- J5 <- numeric(2*L) # empty containers
A1 <- A2 <- A3 <- A4 <- A5 <- numeric(2*L) # empty containers
C1 <- C2 <- C3 <- C4 <- C5 <- numeric(2*L) # empty containers
K1 <- K2 <- K3 <- K4 <- K5 <- numeric(2*L) # empty containers
N1 <- N2 <- N3 <- N4 <- N5 <- numeric(2*L) # empty containers
F1 <- F2 <- F3 <- F4 <- F5 <- numeric(2*L) # empty containers

# perform jackknife resampling of means and mean difference
for (i in 1:L) { # for each JK zone
  T0 <- T15 # create/restore duplicate
  # double weights if JK code is 1 and zero weights if JK code is 0
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 1),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  J1[i] <- wt.mn(T0$PV1, T0$HWt) # PV1 reweighted mean (total)
  J2[i] <- wt.mn(T0$PV2, T0$HWt) # PV2
  J3[i] <- wt.mn(T0$PV3, T0$HWt) # PV3
  J4[i] <- wt.mn(T0$PV4, T0$HWt) # PV4
  J5[i] <- wt.mn(T0$PV5, T0$HWt) # PV5
  A1[i] <- wt.qnt(T0$PV1, T0$HWt, .5) # PV1 reweighted median (total)
  A2[i] <- wt.qnt(T0$PV2, T0$HWt, .5) # PV2
  A3[i] <- wt.qnt(T0$PV3, T0$HWt, .5) # PV3
  A4[i] <- wt.qnt(T0$PV4, T0$HWt, .5) # PV4
  A5[i] <- wt.qnt(T0$PV5, T0$HWt, .5) # PV5
  C1[i] <- wt.mn(T0_F$PV1, T0_F$HWt) # PV1 reweighted mean (females)
  C2[i] <- wt.mn(T0_F$PV2, T0_F$HWt) # PV2
  C3[i] <- wt.mn(T0_F$PV3, T0_F$HWt) # PV3
  C4[i] <- wt.mn(T0_F$PV4, T0_F$HWt) # PV4
  C5[i] <- wt.mn(T0_F$PV5, T0_F$HWt) # PV5
  K1[i] <- wt.qnt(T0_F$PV1, T0_F$HWt, .5) # PV1 reweighted median (females)
  K2[i] <- wt.qnt(T0_F$PV2, T0_F$HWt, .5) # PV2
  K3[i] <- wt.qnt(T0_F$PV3, T0_F$HWt, .5) # PV3
  K4[i] <- wt.qnt(T0_F$PV4, T0_F$HWt, .5) # PV4
  K5[i] <- wt.qnt(T0_F$PV5, T0_F$HWt, .5) # PV5
  N1[i] <- wt.mn(T0_M$PV1, T0_M$HWt) # PV1 reweighted mean (males)
  N2[i] <- wt.mn(T0_M$PV2, T0_M$HWt) # PV2
  N3[i] <- wt.mn(T0_M$PV3, T0_M$HWt) # PV3
  N4[i] <- wt.mn(T0_M$PV4, T0_M$HWt) # PV4
  N5[i] <- wt.mn(T0_M$PV5, T0_M$HWt) # PV5
  F1[i] <- wt.qnt(T0_M$PV1, T0_M$HWt, .5) # PV1 reweighted median (males)
  F2[i] <- wt.qnt(T0_M$PV2, T0_M$HWt, .5) # PV2
  F3[i] <- wt.qnt(T0_M$PV3, T0_M$HWt, .5) # PV3
  F4[i] <- wt.qnt(T0_M$PV4, T0_M$HWt, .5) # PV4
  F5[i] <- wt.qnt(T0_M$PV5, T0_M$HWt, .5) # PV5
  
  T0 <- T15 # restore duplicate
  # double weights if JK code is 0 and zero weights if JK code is 1
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 0),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  J1[i+L] <- wt.mn(T0$PV1, T0$HWt) # PV1 reweighted mean (total)
  J2[i+L] <- wt.mn(T0$PV2, T0$HWt) # PV2
  J3[i+L] <- wt.mn(T0$PV3, T0$HWt) # PV3
  J4[i+L] <- wt.mn(T0$PV4, T0$HWt) # PV4
  J5[i+L] <- wt.mn(T0$PV5, T0$HWt) # PV5
  A1[i+L] <- wt.qnt(T0$PV1, T0$HWt, .5) # PV1 reweighted median (total)
  A2[i+L] <- wt.qnt(T0$PV2, T0$HWt, .5) # PV2
  A3[i+L] <- wt.qnt(T0$PV3, T0$HWt, .5) # PV3
  A4[i+L] <- wt.qnt(T0$PV4, T0$HWt, .5) # PV4
  A5[i+L] <- wt.qnt(T0$PV5, T0$HWt, .5) # PV5
  C1[i+L] <- wt.mn(T0_F$PV1, T0_F$HWt) # PV1 reweighted mean (females)
  C2[i+L] <- wt.mn(T0_F$PV2, T0_F$HWt) # PV2
  C3[i+L] <- wt.mn(T0_F$PV3, T0_F$HWt) # PV3
  C4[i+L] <- wt.mn(T0_F$PV4, T0_F$HWt) # PV4
  C5[i+L] <- wt.mn(T0_F$PV5, T0_F$HWt) # PV5
  K1[i+L] <- wt.qnt(T0_F$PV1, T0_F$HWt, .5) # PV1 reweighted median (females)
  K2[i+L] <- wt.qnt(T0_F$PV2, T0_F$HWt, .5) # PV2
  K3[i+L] <- wt.qnt(T0_F$PV3, T0_F$HWt, .5) # PV3
  K4[i+L] <- wt.qnt(T0_F$PV4, T0_F$HWt, .5) # PV4
  K5[i+L] <- wt.qnt(T0_F$PV5, T0_F$HWt, .5) # PV5
  N1[i+L] <- wt.mn(T0_M$PV1, T0_M$HWt) # PV1 reweighted mean (males)
  N2[i+L] <- wt.mn(T0_M$PV2, T0_M$HWt) # PV2
  N3[i+L] <- wt.mn(T0_M$PV3, T0_M$HWt) # PV3
  N4[i+L] <- wt.mn(T0_M$PV4, T0_M$HWt) # PV4
  N5[i+L] <- wt.mn(T0_M$PV5, T0_M$HWt) # PV5
  F1[i+L] <- wt.qnt(T0_M$PV1, T0_M$HWt, .5) # PV1 reweighted median (males)
  F2[i+L] <- wt.qnt(T0_M$PV2, T0_M$HWt, .5) # PV2
  F3[i+L] <- wt.qnt(T0_M$PV3, T0_M$HWt, .5) # PV3
  F4[i+L] <- wt.qnt(T0_M$PV4, T0_M$HWt, .5) # PV4
  F5[i+L] <- wt.qnt(T0_M$PV5, T0_M$HWt, .5) # PV5
}

# jackknife sampling variance
JSV_MnT <- mean(c(sum((J1-Mn1)^2), # mean (total)
                  sum((J2-Mn2)^2),
                  sum((J3-Mn3)^2),
                  sum((J4-Mn4)^2),
                  sum((J5-Mn5)^2)))/2
JSV_MdT <- mean(c(sum((A1-Md1)^2), # median (total)
                  sum((A2-Md2)^2),
                  sum((A3-Md3)^2),
                  sum((A4-Md4)^2),
                  sum((A5-Md5)^2)))/2
JSV_MnF <- mean(c(sum((C1-Mn1_F)^2), # mean (females)
                  sum((C2-Mn2_F)^2),
                  sum((C3-Mn3_F)^2),
                  sum((C4-Mn4_F)^2),
                  sum((C5-Mn5_F)^2)))/2
JSV_MdF <- mean(c(sum((K1-Md1_F)^2), # median (females)
                  sum((K2-Md2_F)^2),
                  sum((K3-Md3_F)^2),
                  sum((K4-Md4_F)^2),
                  sum((K5-Md5_F)^2)))/2
JSV_MnM <- mean(c(sum((N1-Mn1_M)^2), # mean (males)
                  sum((N2-Mn2_M)^2),
                  sum((N3-Mn3_M)^2),
                  sum((N4-Mn4_M)^2),
                  sum((N5-Mn5_M)^2)))/2
JSV_MdM <- mean(c(sum((F1-Md1_M)^2), # median (males)
                  sum((F2-Md2_M)^2),
                  sum((F3-Md3_M)^2),
                  sum((F4-Md4_M)^2),
                  sum((F5-Md5_M)^2)))/2
JSV_MnDf <- mean(c(sum((N1-C1-MnDf1)^2), # mean difference
                   sum((N2-C2-MnDf2)^2),
                   sum((N3-C3-MnDf3)^2),
                   sum((N4-C4-MnDf4)^2),
                   sum((N5-C5-MnDf5)^2)))/2
JSV_MdDf <- mean(c(sum((F1-K1-MdDf1)^2), # median difference
                   sum((F2-K2-MdDf2)^2),
                   sum((F3-K3-MdDf3)^2),
                   sum((F4-K4-MdDf4)^2),
                   sum((F5-K5-MdDf5)^2)))/2

# imputation variance
IV_MnT <- .3*sum((Mns-Mn)^2) # mean (total)
IV_MdT <- .3*sum((Mds-Md)^2) # median (total)
IV_MnF <- .3*sum((Mns_F-Mn_F)^2) # mean (females)
IV_MdF <- .3*sum((Mds_F-Md_F)^2) # median (females)
IV_MnM <- .3*sum((Mns_M-Mn_M)^2) # mean (males)
IV_MdM <- .3*sum((Mds_M-Md_M)^2) # median (males)
IV_MnDf <- .3*sum((MnDfs-MnDf)^2) # mean difference
IV_MdDf <- .3*sum((MdDfs-MdDf)^2) # median difference


# total variance = jackknife + imputation
TV_MnT <- JSV_MnT + IV_MnT # total variance of mean (total)
SE_MnT <- sqrt(TV_MnT) # standard error of mean (total)

TV_MdT <- JSV_MdT + IV_MdT # total variance of median (total)
SE_MdT <- sqrt(TV_MdT) # standard error of median (total)

TV_MnF <- JSV_MnF + IV_MnF # total variance of mean (females)
SE_MnF <- sqrt(TV_MnF) # standard error of mean (females)

TV_MdF <- JSV_MdF + IV_MdF # total variance of median (females)
SE_MdF <- sqrt(TV_MdF) # standard error of median (females)

TV_MnM <- JSV_MnM + IV_MnM # total variance of mean (males)
SE_MnM <- sqrt(TV_MnM) # standard error of mean (males)

TV_MdM <- JSV_MdM + IV_MdM # total variance of median (males)
SE_MdM <- sqrt(TV_MdM) # standard error of median (males)

TV_MnDf <- JSV_MnDf + IV_MnDf # total variance of mean difference
SE_MnDf <- sqrt(TV_MnDf) # standard error of mean difference

TV_MdDf <- JSV_MdDf + IV_MdDf # total variance of median difference
SE_MdDf <- sqrt(TV_MdDf) # standard error of median difference



#### SEs: LTPRs and LTPR tail-center differences ####

A1 <- A2 <- A3 <- A4 <- A5 <- numeric(2*L) # empty containers
B1 <- B2 <- B3 <- B4 <- B5 <- numeric(2*L) # empty containers
C1 <- C2 <- C3 <- C4 <- C5 <- numeric(2*L) # empty containers
E1 <- E2 <- E3 <- E4 <- E5 <- numeric(2*L) # empty containers
F1 <- F2 <- F3 <- F4 <- F5 <- numeric(2*L) # empty containers
G1 <- G2 <- G3 <- G4 <- G5 <- numeric(2*L) # empty containers
H1 <- H2 <- H3 <- H4 <- H5 <- numeric(2*L) # empty containers
I1 <- I2 <- I3 <- I4 <- I5 <- numeric(2*L) # empty containers
J1 <- J2 <- J3 <- J4 <- J5 <- numeric(2*L) # empty containers
K1 <- K2 <- K3 <- K4 <- K5 <- numeric(2*L) # empty containers
M1 <- M2 <- M3 <- M4 <- M5 <- numeric(2*L) # empty containers
N1 <- N2 <- N3 <- N4 <- N5 <- numeric(2*L) # empty containers
O1 <- O2 <- O3 <- O4 <- O5 <- numeric(2*L) # empty containers
P1 <- P2 <- P3 <- P4 <- P5 <- numeric(2*L) # empty containers
Q1 <- Q2 <- Q3 <- Q4 <- Q5 <- numeric(2*L) # empty containers
R1 <- R2 <- R3 <- R4 <- R5 <- numeric(2*L) # empty containers
S1 <- S2 <- S3 <- S4 <- S5 <- numeric(2*L) # empty containers
T1 <- T2 <- T3 <- T4 <- T5 <- numeric(2*L) # empty containers
V1 <- V2 <- V3 <- V4 <- V5 <- numeric(2*L) # empty containers
W1 <- W2 <- W3 <- W4 <- W5 <- numeric(2*L) # empty containers

# perform jackknife resampling of LTPRs and LTPR tail-center differences
for (i in 1:L) { # for each JK zone
  T0 <- T15 # create/restore duplicate
  # double weights if JK code is 1 and zero weights if JK code is 0
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 1),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  A1[i] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV1, T0$HWt), 'PV1') # PV1 reweighted mean
  A2[i] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV2, T0$HWt), 'PV2') # PV2
  A3[i] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV3, T0$HWt), 'PV3') # PV3
  A4[i] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV4, T0$HWt), 'PV4') # PV4
  A5[i] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV5, T0$HWt), 'PV5') # PV5
  B1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .05), 'PV1') # PV1 reweighted 5th LTPR
  B2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .05), 'PV2') # PV2
  B3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .05), 'PV3') # PV3
  B4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .05), 'PV4') # PV4
  B5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .05), 'PV5') # PV5
  C1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .10), 'PV1') # PV1 reweighted 10th LTPR
  C2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .10), 'PV2') # PV2
  C3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .10), 'PV3') # PV3
  C4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .10), 'PV4') # PV4
  C5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .10), 'PV5') # PV5
  E1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .15), 'PV1') # PV1 reweighted 15th LTPR
  E2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .15), 'PV2') # PV2
  E3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .15), 'PV3') # PV3
  E4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .15), 'PV4') # PV4
  E5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .15), 'PV5') # PV5
  F1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .20), 'PV1') # PV1 reweighted 20th LTPR
  F2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .20), 'PV2') # PV2
  F3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .20), 'PV3') # PV3
  F4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .20), 'PV4') # PV4
  F5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .20), 'PV5') # PV5
  G1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .25), 'PV1') # PV1 reweighted 25th LTPR
  G2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .25), 'PV2') # PV2
  G3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .25), 'PV3') # PV3
  G4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .25), 'PV4') # PV4
  G5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .25), 'PV5') # PV5
  H1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .30), 'PV1') # PV1 reweighted 30th LTPR
  H2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .30), 'PV2') # PV2
  H3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .30), 'PV3') # PV3
  H4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .30), 'PV4') # PV4
  H5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .30), 'PV5') # PV5
  I1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .35), 'PV1') # PV1 reweighted 35th LTPR
  I2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .35), 'PV2') # PV2
  I3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .35), 'PV3') # PV3
  I4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .35), 'PV4') # PV4
  I5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .35), 'PV5') # PV5
  J1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .40), 'PV1') # PV1 reweighted 40th LTPR
  J2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .40), 'PV2') # PV2
  J3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .40), 'PV3') # PV3
  J4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .40), 'PV4') # PV4
  J5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .40), 'PV5') # PV5
  K1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .45), 'PV1') # PV1 reweighted 45th LTPR
  K2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .45), 'PV2') # PV2
  K3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .45), 'PV3') # PV3
  K4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .45), 'PV4') # PV4
  K5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .45), 'PV5') # PV5
  M1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .50), 'PV1') # PV1 reweighted 50th LTPR
  M2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .50), 'PV2') # PV2
  M3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .50), 'PV3') # PV3
  M4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .50), 'PV4') # PV4
  M5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .50), 'PV5') # PV5
  N1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .55), 'PV1') # PV1 reweighted 55th LTPR
  N2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .55), 'PV2') # PV2
  N3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .55), 'PV3') # PV3
  N4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .55), 'PV4') # PV4
  N5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .55), 'PV5') # PV5
  O1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .60), 'PV1') # PV1 reweighted 60th LTPR
  O2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .60), 'PV2') # PV2
  O3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .60), 'PV3') # PV3
  O4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .60), 'PV4') # PV4
  O5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .60), 'PV5') # PV5
  P1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .65), 'PV1') # PV1 reweighted 65th LTPR
  P2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .65), 'PV2') # PV2
  P3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .65), 'PV3') # PV3
  P4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .65), 'PV4') # PV4
  P5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .65), 'PV5') # PV5
  Q1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .70), 'PV1') # PV1 reweighted 70th LTPR
  Q2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .70), 'PV2') # PV2
  Q3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .70), 'PV3') # PV3
  Q4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .70), 'PV4') # PV4
  Q5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .70), 'PV5') # PV5
  R1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .75), 'PV1') # PV1 reweighted 75th LTPR
  R2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .75), 'PV2') # PV2
  R3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .75), 'PV3') # PV3
  R4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .75), 'PV4') # PV4
  R5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .75), 'PV5') # PV5
  S1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .80), 'PV1') # PV1 reweighted 80th LTPR
  S2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .80), 'PV2') # PV2
  S3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .80), 'PV3') # PV3
  S4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .80), 'PV4') # PV4
  S5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .80), 'PV5') # PV5
  T1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .85), 'PV1') # PV1 reweighted 85th LTPR
  T2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .85), 'PV2') # PV2
  T3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .85), 'PV3') # PV3
  T4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .85), 'PV4') # PV4
  T5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .85), 'PV5') # PV5
  V1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .90), 'PV1') # PV1 reweighted 90th LTPR
  V2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .90), 'PV2') # PV2
  V3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .90), 'PV3') # PV3
  V4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .90), 'PV4') # PV4
  V5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .90), 'PV5') # PV5
  W1[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .95), 'PV1') # PV1 reweighted 95th LTPR
  W2[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .95), 'PV2') # PV2
  W3[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .95), 'PV3') # PV3
  W4[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .95), 'PV4') # PV4
  W5[i] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .95), 'PV5') # PV5
  
  T0 <- T15 # restore duplicate
  # double weights if JK code is 0 and zero weights if JK code is 1
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 0),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  A1[i+L] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV1, T0$HWt), 'PV1') # PV1 reweighted mean
  A2[i+L] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV2, T0$HWt), 'PV2') # PV2
  A3[i+L] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV3, T0$HWt), 'PV3') # PV3
  A4[i+L] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV4, T0$HWt), 'PV4') # PV4
  A5[i+L] <- LTPRfn(T0_M, T0_F, wt.mn(T0$PV5, T0$HWt), 'PV5') # PV5
  B1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .05), 'PV1') # PV1 reweighted 5th LTPR
  B2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .05), 'PV2') # PV2
  B3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .05), 'PV3') # PV3
  B4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .05), 'PV4') # PV4
  B5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .05), 'PV5') # PV5
  C1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .10), 'PV1') # PV1 reweighted 10th LTPR
  C2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .10), 'PV2') # PV2
  C3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .10), 'PV3') # PV3
  C4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .10), 'PV4') # PV4
  C5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .10), 'PV5') # PV5
  E1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .15), 'PV1') # PV1 reweighted 15th LTPR
  E2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .15), 'PV2') # PV2
  E3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .15), 'PV3') # PV3
  E4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .15), 'PV4') # PV4
  E5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .15), 'PV5') # PV5
  F1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .20), 'PV1') # PV1 reweighted 20th LTPR
  F2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .20), 'PV2') # PV2
  F3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .20), 'PV3') # PV3
  F4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .20), 'PV4') # PV4
  F5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .20), 'PV5') # PV5
  G1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .25), 'PV1') # PV1 reweighted 25th LTPR
  G2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .25), 'PV2') # PV2
  G3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .25), 'PV3') # PV3
  G4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .25), 'PV4') # PV4
  G5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .25), 'PV5') # PV5
  H1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .30), 'PV1') # PV1 reweighted 30th LTPR
  H2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .30), 'PV2') # PV2
  H3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .30), 'PV3') # PV3
  H4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .30), 'PV4') # PV4
  H5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .30), 'PV5') # PV5
  I1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .35), 'PV1') # PV1 reweighted 35th LTPR
  I2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .35), 'PV2') # PV2
  I3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .35), 'PV3') # PV3
  I4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .35), 'PV4') # PV4
  I5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .35), 'PV5') # PV5
  J1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .40), 'PV1') # PV1 reweighted 40th LTPR
  J2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .40), 'PV2') # PV2
  J3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .40), 'PV3') # PV3
  J4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .40), 'PV4') # PV4
  J5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .40), 'PV5') # PV5
  K1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .45), 'PV1') # PV1 reweighted 45th LTPR
  K2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .45), 'PV2') # PV2
  K3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .45), 'PV3') # PV3
  K4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .45), 'PV4') # PV4
  K5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .45), 'PV5') # PV5
  M1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .50), 'PV1') # PV1 reweighted 50th LTPR
  M2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .50), 'PV2') # PV2
  M3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .50), 'PV3') # PV3
  M4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .50), 'PV4') # PV4
  M5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .50), 'PV5') # PV5
  N1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .55), 'PV1') # PV1 reweighted 55th LTPR
  N2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .55), 'PV2') # PV2
  N3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .55), 'PV3') # PV3
  N4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .55), 'PV4') # PV4
  N5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .55), 'PV5') # PV5
  O1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .60), 'PV1') # PV1 reweighted 60th LTPR
  O2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .60), 'PV2') # PV2
  O3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .60), 'PV3') # PV3
  O4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .60), 'PV4') # PV4
  O5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .60), 'PV5') # PV5
  P1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .65), 'PV1') # PV1 reweighted 65th LTPR
  P2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .65), 'PV2') # PV2
  P3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .65), 'PV3') # PV3
  P4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .65), 'PV4') # PV4
  P5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .65), 'PV5') # PV5
  Q1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .70), 'PV1') # PV1 reweighted 70th LTPR
  Q2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .70), 'PV2') # PV2
  Q3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .70), 'PV3') # PV3
  Q4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .70), 'PV4') # PV4
  Q5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .70), 'PV5') # PV5
  R1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .75), 'PV1') # PV1 reweighted 75th LTPR
  R2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .75), 'PV2') # PV2
  R3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .75), 'PV3') # PV3
  R4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .75), 'PV4') # PV4
  R5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .75), 'PV5') # PV5
  S1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .80), 'PV1') # PV1 reweighted 80th LTPR
  S2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .80), 'PV2') # PV2
  S3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .80), 'PV3') # PV3
  S4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .80), 'PV4') # PV4
  S5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .80), 'PV5') # PV5
  T1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .85), 'PV1') # PV1 reweighted 85th LTPR
  T2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .85), 'PV2') # PV2
  T3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .85), 'PV3') # PV3
  T4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .85), 'PV4') # PV4
  T5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .85), 'PV5') # PV5
  V1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .90), 'PV1') # PV1 reweighted 90th LTPR
  V2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .90), 'PV2') # PV2
  V3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .90), 'PV3') # PV3
  V4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .90), 'PV4') # PV4
  V5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .90), 'PV5') # PV5
  W1[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV1, T0$HWt, .95), 'PV1') # PV1 reweighted 95th LTPR
  W2[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV2, T0$HWt, .95), 'PV2') # PV2
  W3[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV3, T0$HWt, .95), 'PV3') # PV3
  W4[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV4, T0$HWt, .95), 'PV4') # PV4
  W5[i+L] <- LTPRfn(T0_M, T0_F, wt.qnt(T0$PV5, T0$HWt, .95), 'PV5') # PV5
  
  if (i %% 5 == 0) {print(paste0(i, '/', L, ' at ', Sys.time()), quote = F)} # print updates
}

# jackknife sampling variance
JSV_Mn <- mean(c(sum((A1-LTPRMn_1)^2), # mean LTPR
                 sum((A2-LTPRMn_2)^2),
                 sum((A3-LTPRMn_3)^2),
                 sum((A4-LTPRMn_4)^2),
                 sum((A5-LTPRMn_5)^2)))/2
JSV_05 <- mean(c(sum((B1-LTPR05_1)^2), # 5th LTPR
                 sum((B2-LTPR05_2)^2),
                 sum((B3-LTPR05_3)^2),
                 sum((B4-LTPR05_4)^2),
                 sum((B5-LTPR05_5)^2)))/2
JSV_10 <- mean(c(sum((C1-LTPR10_1)^2), # 10th LTPR
                 sum((C2-LTPR10_2)^2),
                 sum((C3-LTPR10_3)^2),
                 sum((C4-LTPR10_4)^2),
                 sum((C5-LTPR10_5)^2)))/2
JSV_15 <- mean(c(sum((E1-LTPR15_1)^2), # 15th LTPR
                 sum((E2-LTPR15_2)^2),
                 sum((E3-LTPR15_3)^2),
                 sum((E4-LTPR15_4)^2),
                 sum((E5-LTPR15_5)^2)))/2
JSV_20 <- mean(c(sum((F1-LTPR20_1)^2), # 20th LTPR
                 sum((F2-LTPR20_2)^2),
                 sum((F3-LTPR20_3)^2),
                 sum((F4-LTPR20_4)^2),
                 sum((F5-LTPR20_5)^2)))/2
JSV_25 <- mean(c(sum((G1-LTPR25_1)^2), # 25th LTPR
                 sum((G2-LTPR25_2)^2),
                 sum((G3-LTPR25_3)^2),
                 sum((G4-LTPR25_4)^2),
                 sum((G5-LTPR25_5)^2)))/2
JSV_30 <- mean(c(sum((H1-LTPR30_1)^2), # 30th LTPR
                 sum((H2-LTPR30_2)^2),
                 sum((H3-LTPR30_3)^2),
                 sum((H4-LTPR30_4)^2),
                 sum((H5-LTPR30_5)^2)))/2
JSV_35 <- mean(c(sum((I1-LTPR35_1)^2), # 35th LTPR
                 sum((I2-LTPR35_2)^2),
                 sum((I3-LTPR35_3)^2),
                 sum((I4-LTPR35_4)^2),
                 sum((I5-LTPR35_5)^2)))/2
JSV_40 <- mean(c(sum((J1-LTPR40_1)^2), # 40th LTPR
                 sum((J2-LTPR40_2)^2),
                 sum((J3-LTPR40_3)^2),
                 sum((J4-LTPR40_4)^2),
                 sum((J5-LTPR40_5)^2)))/2
JSV_45 <- mean(c(sum((K1-LTPR45_1)^2), # 45th LTPR
                 sum((K2-LTPR45_2)^2),
                 sum((K3-LTPR45_3)^2),
                 sum((K4-LTPR45_4)^2),
                 sum((K5-LTPR45_5)^2)))/2
JSV_50 <- mean(c(sum((M1-LTPR50_1)^2), # 50th LTPR
                 sum((M2-LTPR50_2)^2),
                 sum((M3-LTPR50_3)^2),
                 sum((M4-LTPR50_4)^2),
                 sum((M5-LTPR50_5)^2)))/2
JSV_55 <- mean(c(sum((N1-LTPR55_1)^2), # 55th LTPR
                 sum((N2-LTPR55_2)^2),
                 sum((N3-LTPR55_3)^2),
                 sum((N4-LTPR55_4)^2),
                 sum((N5-LTPR55_5)^2)))/2
JSV_60 <- mean(c(sum((O1-LTPR60_1)^2), # 60th LTPR
                 sum((O2-LTPR60_2)^2),
                 sum((O3-LTPR60_3)^2),
                 sum((O4-LTPR60_4)^2),
                 sum((O5-LTPR60_5)^2)))/2
JSV_65 <- mean(c(sum((P1-LTPR65_1)^2), # 65th LTPR
                 sum((P2-LTPR65_2)^2),
                 sum((P3-LTPR65_3)^2),
                 sum((P4-LTPR65_4)^2),
                 sum((P5-LTPR65_5)^2)))/2
JSV_70 <- mean(c(sum((Q1-LTPR70_1)^2), # 70th LTPR
                 sum((Q2-LTPR70_2)^2),
                 sum((Q3-LTPR70_3)^2),
                 sum((Q4-LTPR70_4)^2),
                 sum((Q5-LTPR70_5)^2)))/2
JSV_75 <- mean(c(sum((R1-LTPR75_1)^2), # 75th LTPR
                 sum((R2-LTPR75_2)^2),
                 sum((R3-LTPR75_3)^2),
                 sum((R4-LTPR75_4)^2),
                 sum((R5-LTPR75_5)^2)))/2
JSV_80 <- mean(c(sum((S1-LTPR80_1)^2), # 80th LTPR
                 sum((S2-LTPR80_2)^2),
                 sum((S3-LTPR80_3)^2),
                 sum((S4-LTPR80_4)^2),
                 sum((S5-LTPR80_5)^2)))/2
JSV_85 <- mean(c(sum((T1-LTPR85_1)^2), # 85th LTPR
                 sum((T2-LTPR85_2)^2),
                 sum((T3-LTPR85_3)^2),
                 sum((T4-LTPR85_4)^2),
                 sum((T5-LTPR85_5)^2)))/2
JSV_90 <- mean(c(sum((V1-LTPR90_1)^2), # 90th LTPR
                 sum((V2-LTPR90_2)^2),
                 sum((V3-LTPR90_3)^2),
                 sum((V4-LTPR90_4)^2),
                 sum((V5-LTPR90_5)^2)))/2
JSV_95 <- mean(c(sum((W1-LTPR95_1)^2), # 95th LTPR
                 sum((W2-LTPR95_2)^2),
                 sum((W3-LTPR95_3)^2),
                 sum((W4-LTPR95_4)^2),
                 sum((W5-LTPR95_5)^2)))/2
JSV_Md95T <- mean(c(sum((W1-M1-Md95T_1)^2), # Md95T
                    sum((W2-M2-Md95T_2)^2),
                    sum((W3-M3-Md95T_3)^2),
                    sum((W4-M4-Md95T_4)^2),
                    sum((W5-M5-Md95T_5)^2)))/2
JSV_Md90T <- mean(c(sum((V1-M1-Md90T_1)^2), # Md90T
                    sum((V2-M2-Md90T_2)^2),
                    sum((V3-M3-Md90T_3)^2),
                    sum((V4-M4-Md90T_4)^2),
                    sum((V5-M5-Md90T_5)^2)))/2
JSV_Md10T <- mean(c(sum((M1-C1-Md10T_1)^2), # Md10T
                    sum((M2-C2-Md10T_2)^2),
                    sum((M3-C3-Md10T_3)^2),
                    sum((M4-C4-Md10T_4)^2),
                    sum((M5-C5-Md10T_5)^2)))/2
JSV_Md05T <- mean(c(sum((M1-B1-Md05T_1)^2), # Md05T
                    sum((M2-B2-Md05T_2)^2),
                    sum((M3-B3-Md05T_3)^2),
                    sum((M4-B4-Md05T_4)^2),
                    sum((M5-B5-Md05T_5)^2)))/2
JSV_Mn95T <- mean(c(sum((W1-A1-Mn95T_1)^2), # Mn95T
                    sum((W2-A2-Mn95T_2)^2),
                    sum((W3-A3-Mn95T_3)^2),
                    sum((W4-A4-Mn95T_4)^2),
                    sum((W5-A5-Mn95T_5)^2)))/2
JSV_Mn90T <- mean(c(sum((V1-A1-Mn90T_1)^2), # Mn90T
                    sum((V2-A2-Mn90T_2)^2),
                    sum((V3-A3-Mn90T_3)^2),
                    sum((V4-A4-Mn90T_4)^2),
                    sum((V5-A5-Mn90T_5)^2)))/2
JSV_Mn10T <- mean(c(sum((A1-C1-Mn10T_1)^2), # Mn10T
                    sum((A2-C2-Mn10T_2)^2),
                    sum((A3-C3-Mn10T_3)^2),
                    sum((A4-C4-Mn10T_4)^2),
                    sum((A5-C5-Mn10T_5)^2)))/2
JSV_Mn05T <- mean(c(sum((A1-B1-Mn05T_1)^2), # Mn05T
                    sum((A2-B2-Mn05T_2)^2),
                    sum((A3-B3-Mn05T_3)^2),
                    sum((A4-B4-Mn05T_4)^2),
                    sum((A5-B5-Mn05T_5)^2)))/2

# imputation variance
IV_Mn <- .3*sum((LTPRMns-LTPRMn)^2) # mean LTPR
IV_05 <- .3*sum((LTPR05s-LTPR05)^2) # 5th LTPR
IV_10 <- .3*sum((LTPR10s-LTPR10)^2) # 10th LTPR
IV_15 <- .3*sum((LTPR15s-LTPR15)^2) # 15th LTPR
IV_20 <- .3*sum((LTPR20s-LTPR20)^2) # 20th LTPR
IV_25 <- .3*sum((LTPR25s-LTPR25)^2) # 25th LTPR
IV_30 <- .3*sum((LTPR30s-LTPR30)^2) # 30th LTPR
IV_35 <- .3*sum((LTPR35s-LTPR35)^2) # 35th LTPR
IV_40 <- .3*sum((LTPR40s-LTPR40)^2) # 40th LTPR
IV_45 <- .3*sum((LTPR45s-LTPR45)^2) # 45th LTPR
IV_50 <- .3*sum((LTPR50s-LTPR50)^2) # 50th LTPR
IV_55 <- .3*sum((LTPR55s-LTPR55)^2) # 55th LTPR
IV_60 <- .3*sum((LTPR60s-LTPR60)^2) # 60th LTPR
IV_65 <- .3*sum((LTPR65s-LTPR65)^2) # 65th LTPR
IV_70 <- .3*sum((LTPR70s-LTPR70)^2) # 70th LTPR
IV_75 <- .3*sum((LTPR75s-LTPR75)^2) # 75th LTPR
IV_80 <- .3*sum((LTPR80s-LTPR80)^2) # 80th LTPR
IV_85 <- .3*sum((LTPR85s-LTPR85)^2) # 85th LTPR
IV_90 <- .3*sum((LTPR90s-LTPR90)^2) # 90th LTPR
IV_95 <- .3*sum((LTPR95s-LTPR95)^2) # 95th LTPR
IV_Md95T <- .3*sum((Md95Ts-Md95T)^2) # Md95T
IV_Md90T <- .3*sum((Md90Ts-Md90T)^2) # Md90T
IV_Md10T <- .3*sum((Md10Ts-Md10T)^2) # Md10T
IV_Md05T <- .3*sum((Md05Ts-Md05T)^2) # Md05T
IV_Mn95T <- .3*sum((Mn95Ts-Mn95T)^2) # Mn95T
IV_Mn90T <- .3*sum((Mn90Ts-Mn90T)^2) # Mn90T
IV_Mn10T <- .3*sum((Mn10Ts-Mn10T)^2) # Mn10T
IV_Mn05T <- .3*sum((Mn05Ts-Mn05T)^2) # Mn05T


# total variance = jackknife + imputation
TV_Mn <- JSV_Mn + IV_Mn # total variance of mean LTPR
SE_Mn <- sqrt(TV_Mn) # standard error of mean LTPR

TV_05T <- JSV_05 + IV_05 # total variance of 5th percentile LTPR
SE_05T <- sqrt(TV_05T) # standard error of 5th percentile LTPR

TV_10T <- JSV_10 + IV_10 # total variance of 10th percentile LTPR
SE_10T <- sqrt(TV_10T) # standard error of 10th percentile LTPR

TV_15T <- JSV_15 + IV_15 # total variance of 15th percentile LTPR
SE_15T <- sqrt(TV_15T) # standard error of 15th percentile LTPR

TV_20T <- JSV_20 + IV_20 # total variance of 20th percentile LTPR
SE_20T <- sqrt(TV_20T) # standard error of 20th percentile LTPR

TV_25T <- JSV_25 + IV_25 # total variance of 25th percentile LTPR
SE_25T <- sqrt(TV_25T) # standard error of 25th percentile LTPR

TV_30T <- JSV_30 + IV_30 # total variance of 30th percentile LTPR
SE_30T <- sqrt(TV_30T) # standard error of 30th percentile LTPR

TV_35T <- JSV_35 + IV_35 # total variance of 35th percentile LTPR
SE_35T <- sqrt(TV_35T) # standard error of 35th percentile LTPR

TV_40T <- JSV_40 + IV_40 # total variance of 40th percentile LTPR
SE_40T <- sqrt(TV_40T) # standard error of 40th percentile LTPR

TV_45T <- JSV_45 + IV_45 # total variance of 45th percentile LTPR
SE_45T <- sqrt(TV_45T) # standard error of 45th percentile LTPR

TV_50T <- JSV_50 + IV_50 # total variance of 50th percentile LTPR
SE_50T <- sqrt(TV_50T) # standard error of 50th percentile LTPR

TV_55T <- JSV_55 + IV_55 # total variance of 55th percentile LTPR
SE_55T <- sqrt(TV_55T) # standard error of 55th percentile LTPR

TV_60T <- JSV_60 + IV_60 # total variance of 60th percentile LTPR
SE_60T <- sqrt(TV_60T) # standard error of 60th percentile LTPR

TV_65T <- JSV_65 + IV_65 # total variance of 65th percentile LTPR
SE_65T <- sqrt(TV_65T) # standard error of 65th percentile LTPR

TV_70T <- JSV_70 + IV_70 # total variance of 70th percentile LTPR
SE_70T <- sqrt(TV_70T) # standard error of 70th percentile LTPR

TV_75T <- JSV_75 + IV_75 # total variance of 75th percentile LTPR
SE_75T <- sqrt(TV_75T) # standard error of 75th percentile LTPR

TV_80T <- JSV_80 + IV_80 # total variance of 80th percentile LTPR
SE_80T <- sqrt(TV_80T) # standard error of 80th percentile LTPR

TV_85T <- JSV_85 + IV_85 # total variance of 85th percentile LTPR
SE_85T <- sqrt(TV_85T) # standard error of 85th percentile LTPR

TV_90T <- JSV_90 + IV_90 # total variance of 90th percentile LTPR
SE_90T <- sqrt(TV_90T) # standard error of 90th percentile LTPR

TV_95T <- JSV_95 + IV_95 # total variance of 95th percentile LTPR
SE_95T <- sqrt(TV_95T) # standard error of 95th percentile LTPR

TV_Md95T <- JSV_Md95T + IV_Md95T # total variance of Md95T
SE_Md95T <- sqrt(TV_Md95T) # standard error of Md95T

TV_Md90T <- JSV_Md90T + IV_Md90T # total variance of Md90T
SE_Md90T <- sqrt(TV_Md90T) # standard error of Md90T

TV_Md10T <- JSV_Md10T + IV_Md10T # total variance of Md10T
SE_Md10T <- sqrt(TV_Md10T) # standard error of Md10T

TV_Md05T <- JSV_Md05T + IV_Md05T # total variance of Md05T
SE_Md05T <- sqrt(TV_Md05T) # standard error of Md05T

TV_Mn95T <- JSV_Mn95T + IV_Mn95T # total variance of Mn95T
SE_Mn95T <- sqrt(TV_Mn95T) # standard error of Mn95T

TV_Mn90T <- JSV_Mn90T + IV_Mn90T # total variance of Mn90T
SE_Mn90T <- sqrt(TV_Mn90T) # standard error of Mn90T

TV_Mn10T <- JSV_Mn10T + IV_Mn10T # total variance of Mn10T
SE_Mn10T <- sqrt(TV_Mn10T) # standard error of Mn10T

TV_Mn05T <- JSV_Mn05T + IV_Mn05T # total variance of Mn05T
SE_Mn05T <- sqrt(TV_Mn05T) # standard error of Mn05T



#### SEs: LU3Rs and LU3R tail-center differences ####

B1 <- B2 <- B3 <- B4 <- B5 <- numeric(2*L) # empty containers
C1 <- C2 <- C3 <- C4 <- C5 <- numeric(2*L) # empty containers
E1 <- E2 <- E3 <- E4 <- E5 <- numeric(2*L) # empty containers
F1 <- F2 <- F3 <- F4 <- F5 <- numeric(2*L) # empty containers
G1 <- G2 <- G3 <- G4 <- G5 <- numeric(2*L) # empty containers
H1 <- H2 <- H3 <- H4 <- H5 <- numeric(2*L) # empty containers
I1 <- I2 <- I3 <- I4 <- I5 <- numeric(2*L) # empty containers
J1 <- J2 <- J3 <- J4 <- J5 <- numeric(2*L) # empty containers
K1 <- K2 <- K3 <- K4 <- K5 <- numeric(2*L) # empty containers
M1 <- M2 <- M3 <- M4 <- M5 <- numeric(2*L) # empty containers
N1 <- N2 <- N3 <- N4 <- N5 <- numeric(2*L) # empty containers
O1 <- O2 <- O3 <- O4 <- O5 <- numeric(2*L) # empty containers
P1 <- P2 <- P3 <- P4 <- P5 <- numeric(2*L) # empty containers
Q1 <- Q2 <- Q3 <- Q4 <- Q5 <- numeric(2*L) # empty containers
R1 <- R2 <- R3 <- R4 <- R5 <- numeric(2*L) # empty containers
S1 <- S2 <- S3 <- S4 <- S5 <- numeric(2*L) # empty containers
T1 <- T2 <- T3 <- T4 <- T5 <- numeric(2*L) # empty containers
V1 <- V2 <- V3 <- V4 <- V5 <- numeric(2*L) # empty containers
W1 <- W2 <- W3 <- W4 <- W5 <- numeric(2*L) # empty containers

# perform jackknife resampling of LU3Rs and LU3R tail-center differences
for (i in 1:L) { # for each JK zone
  T0 <- T15 # create/restore duplicate
  # double weights if JK code is 1 and zero weights if JK code is 0
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 1),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  B1[i] <- LU3Rfn(T0_M, T0_F, .05, 'PV1') # PV1 reweighted 5th LU3R
  B2[i] <- LU3Rfn(T0_M, T0_F, .05, 'PV2') # PV2
  B3[i] <- LU3Rfn(T0_M, T0_F, .05, 'PV3') # PV3
  B4[i] <- LU3Rfn(T0_M, T0_F, .05, 'PV4') # PV4
  B5[i] <- LU3Rfn(T0_M, T0_F, .05, 'PV5') # PV5
  C1[i] <- LU3Rfn(T0_M, T0_F, .10, 'PV1') # PV1 reweighted 10th LU3R
  C2[i] <- LU3Rfn(T0_M, T0_F, .10, 'PV2') # PV2
  C3[i] <- LU3Rfn(T0_M, T0_F, .10, 'PV3') # PV3
  C4[i] <- LU3Rfn(T0_M, T0_F, .10, 'PV4') # PV4
  C5[i] <- LU3Rfn(T0_M, T0_F, .10, 'PV5') # PV5
  E1[i] <- LU3Rfn(T0_M, T0_F, .15, 'PV1') # PV1 reweighted 15th LU3R
  E2[i] <- LU3Rfn(T0_M, T0_F, .15, 'PV2') # PV2
  E3[i] <- LU3Rfn(T0_M, T0_F, .15, 'PV3') # PV3
  E4[i] <- LU3Rfn(T0_M, T0_F, .15, 'PV4') # PV4
  E5[i] <- LU3Rfn(T0_M, T0_F, .15, 'PV5') # PV5
  F1[i] <- LU3Rfn(T0_M, T0_F, .20, 'PV1') # PV1 reweighted 20th LU3R
  F2[i] <- LU3Rfn(T0_M, T0_F, .20, 'PV2') # PV2
  F3[i] <- LU3Rfn(T0_M, T0_F, .20, 'PV3') # PV3
  F4[i] <- LU3Rfn(T0_M, T0_F, .20, 'PV4') # PV4
  F5[i] <- LU3Rfn(T0_M, T0_F, .20, 'PV5') # PV5
  G1[i] <- LU3Rfn(T0_M, T0_F, .25, 'PV1') # PV1 reweighted 25th LU3R
  G2[i] <- LU3Rfn(T0_M, T0_F, .25, 'PV2') # PV2
  G3[i] <- LU3Rfn(T0_M, T0_F, .25, 'PV3') # PV3
  G4[i] <- LU3Rfn(T0_M, T0_F, .25, 'PV4') # PV4
  G5[i] <- LU3Rfn(T0_M, T0_F, .25, 'PV5') # PV5
  H1[i] <- LU3Rfn(T0_M, T0_F, .30, 'PV1') # PV1 reweighted 30th LU3R
  H2[i] <- LU3Rfn(T0_M, T0_F, .30, 'PV2') # PV2
  H3[i] <- LU3Rfn(T0_M, T0_F, .30, 'PV3') # PV3
  H4[i] <- LU3Rfn(T0_M, T0_F, .30, 'PV4') # PV4
  H5[i] <- LU3Rfn(T0_M, T0_F, .30, 'PV5') # PV5
  I1[i] <- LU3Rfn(T0_M, T0_F, .35, 'PV1') # PV1 reweighted 35th LU3R
  I2[i] <- LU3Rfn(T0_M, T0_F, .35, 'PV2') # PV2
  I3[i] <- LU3Rfn(T0_M, T0_F, .35, 'PV3') # PV3
  I4[i] <- LU3Rfn(T0_M, T0_F, .35, 'PV4') # PV4
  I5[i] <- LU3Rfn(T0_M, T0_F, .35, 'PV5') # PV5
  J1[i] <- LU3Rfn(T0_M, T0_F, .40, 'PV1') # PV1 reweighted 40th LU3R
  J2[i] <- LU3Rfn(T0_M, T0_F, .40, 'PV2') # PV2
  J3[i] <- LU3Rfn(T0_M, T0_F, .40, 'PV3') # PV3
  J4[i] <- LU3Rfn(T0_M, T0_F, .40, 'PV4') # PV4
  J5[i] <- LU3Rfn(T0_M, T0_F, .40, 'PV5') # PV5
  K1[i] <- LU3Rfn(T0_M, T0_F, .45, 'PV1') # PV1 reweighted 45th LU3R
  K2[i] <- LU3Rfn(T0_M, T0_F, .45, 'PV2') # PV2
  K3[i] <- LU3Rfn(T0_M, T0_F, .45, 'PV3') # PV3
  K4[i] <- LU3Rfn(T0_M, T0_F, .45, 'PV4') # PV4
  K5[i] <- LU3Rfn(T0_M, T0_F, .45, 'PV5') # PV5
  M1[i] <- LU3Rfn(T0_M, T0_F, .50, 'PV1') # PV1 reweighted 50th LU3R
  M2[i] <- LU3Rfn(T0_M, T0_F, .50, 'PV2') # PV2
  M3[i] <- LU3Rfn(T0_M, T0_F, .50, 'PV3') # PV3
  M4[i] <- LU3Rfn(T0_M, T0_F, .50, 'PV4') # PV4
  M5[i] <- LU3Rfn(T0_M, T0_F, .50, 'PV5') # PV5
  N1[i] <- LU3Rfn(T0_M, T0_F, .55, 'PV1') # PV1 reweighted 55th LU3R
  N2[i] <- LU3Rfn(T0_M, T0_F, .55, 'PV2') # PV2
  N3[i] <- LU3Rfn(T0_M, T0_F, .55, 'PV3') # PV3
  N4[i] <- LU3Rfn(T0_M, T0_F, .55, 'PV4') # PV4
  N5[i] <- LU3Rfn(T0_M, T0_F, .55, 'PV5') # PV5
  O1[i] <- LU3Rfn(T0_M, T0_F, .60, 'PV1') # PV1 reweighted 60th LU3R
  O2[i] <- LU3Rfn(T0_M, T0_F, .60, 'PV2') # PV2
  O3[i] <- LU3Rfn(T0_M, T0_F, .60, 'PV3') # PV3
  O4[i] <- LU3Rfn(T0_M, T0_F, .60, 'PV4') # PV4
  O5[i] <- LU3Rfn(T0_M, T0_F, .60, 'PV5') # PV5
  P1[i] <- LU3Rfn(T0_M, T0_F, .65, 'PV1') # PV1 reweighted 65th LU3R
  P2[i] <- LU3Rfn(T0_M, T0_F, .65, 'PV2') # PV2
  P3[i] <- LU3Rfn(T0_M, T0_F, .65, 'PV3') # PV3
  P4[i] <- LU3Rfn(T0_M, T0_F, .65, 'PV4') # PV4
  P5[i] <- LU3Rfn(T0_M, T0_F, .65, 'PV5') # PV5
  Q1[i] <- LU3Rfn(T0_M, T0_F, .70, 'PV1') # PV1 reweighted 70th LU3R
  Q2[i] <- LU3Rfn(T0_M, T0_F, .70, 'PV2') # PV2
  Q3[i] <- LU3Rfn(T0_M, T0_F, .70, 'PV3') # PV3
  Q4[i] <- LU3Rfn(T0_M, T0_F, .70, 'PV4') # PV4
  Q5[i] <- LU3Rfn(T0_M, T0_F, .70, 'PV5') # PV5
  R1[i] <- LU3Rfn(T0_M, T0_F, .75, 'PV1') # PV1 reweighted 75th LU3R
  R2[i] <- LU3Rfn(T0_M, T0_F, .75, 'PV2') # PV2
  R3[i] <- LU3Rfn(T0_M, T0_F, .75, 'PV3') # PV3
  R4[i] <- LU3Rfn(T0_M, T0_F, .75, 'PV4') # PV4
  R5[i] <- LU3Rfn(T0_M, T0_F, .75, 'PV5') # PV5
  S1[i] <- LU3Rfn(T0_M, T0_F, .80, 'PV1') # PV1 reweighted 80th LU3R
  S2[i] <- LU3Rfn(T0_M, T0_F, .80, 'PV2') # PV2
  S3[i] <- LU3Rfn(T0_M, T0_F, .80, 'PV3') # PV3
  S4[i] <- LU3Rfn(T0_M, T0_F, .80, 'PV4') # PV4
  S5[i] <- LU3Rfn(T0_M, T0_F, .80, 'PV5') # PV5
  T1[i] <- LU3Rfn(T0_M, T0_F, .85, 'PV1') # PV1 reweighted 85th LU3R
  T2[i] <- LU3Rfn(T0_M, T0_F, .85, 'PV2') # PV2
  T3[i] <- LU3Rfn(T0_M, T0_F, .85, 'PV3') # PV3
  T4[i] <- LU3Rfn(T0_M, T0_F, .85, 'PV4') # PV4
  T5[i] <- LU3Rfn(T0_M, T0_F, .85, 'PV5') # PV5
  V1[i] <- LU3Rfn(T0_M, T0_F, .90, 'PV1') # PV1 reweighted 90th LU3R
  V2[i] <- LU3Rfn(T0_M, T0_F, .90, 'PV2') # PV2
  V3[i] <- LU3Rfn(T0_M, T0_F, .90, 'PV3') # PV3
  V4[i] <- LU3Rfn(T0_M, T0_F, .90, 'PV4') # PV4
  V5[i] <- LU3Rfn(T0_M, T0_F, .90, 'PV5') # PV5
  W1[i] <- LU3Rfn(T0_M, T0_F, .95, 'PV1') # PV1 reweighted 95th LU3R
  W2[i] <- LU3Rfn(T0_M, T0_F, .95, 'PV2') # PV2
  W3[i] <- LU3Rfn(T0_M, T0_F, .95, 'PV3') # PV3
  W4[i] <- LU3Rfn(T0_M, T0_F, .95, 'PV4') # PV4
  W5[i] <- LU3Rfn(T0_M, T0_F, .95, 'PV5') # PV5
  
  T0 <- T15 # restore duplicate
  # double weights if JK code is 0 and zero weights if JK code is 1
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 0),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  B1[i+L] <- LU3Rfn(T0_M, T0_F, .05, 'PV1') # PV1 reweighted 5th LU3R
  B2[i+L] <- LU3Rfn(T0_M, T0_F, .05, 'PV2') # PV2
  B3[i+L] <- LU3Rfn(T0_M, T0_F, .05, 'PV3') # PV3
  B4[i+L] <- LU3Rfn(T0_M, T0_F, .05, 'PV4') # PV4
  B5[i+L] <- LU3Rfn(T0_M, T0_F, .05, 'PV5') # PV5
  C1[i+L] <- LU3Rfn(T0_M, T0_F, .10, 'PV1') # PV1 reweighted 10th LU3R
  C2[i+L] <- LU3Rfn(T0_M, T0_F, .10, 'PV2') # PV2
  C3[i+L] <- LU3Rfn(T0_M, T0_F, .10, 'PV3') # PV3
  C4[i+L] <- LU3Rfn(T0_M, T0_F, .10, 'PV4') # PV4
  C5[i+L] <- LU3Rfn(T0_M, T0_F, .10, 'PV5') # PV5
  E1[i+L] <- LU3Rfn(T0_M, T0_F, .15, 'PV1') # PV1 reweighted 15th LU3R
  E2[i+L] <- LU3Rfn(T0_M, T0_F, .15, 'PV2') # PV2
  E3[i+L] <- LU3Rfn(T0_M, T0_F, .15, 'PV3') # PV3
  E4[i+L] <- LU3Rfn(T0_M, T0_F, .15, 'PV4') # PV4
  E5[i+L] <- LU3Rfn(T0_M, T0_F, .15, 'PV5') # PV5
  F1[i+L] <- LU3Rfn(T0_M, T0_F, .20, 'PV1') # PV1 reweighted 20th LU3R
  F2[i+L] <- LU3Rfn(T0_M, T0_F, .20, 'PV2') # PV2
  F3[i+L] <- LU3Rfn(T0_M, T0_F, .20, 'PV3') # PV3
  F4[i+L] <- LU3Rfn(T0_M, T0_F, .20, 'PV4') # PV4
  F5[i+L] <- LU3Rfn(T0_M, T0_F, .20, 'PV5') # PV5
  G1[i+L] <- LU3Rfn(T0_M, T0_F, .25, 'PV1') # PV1 reweighted 25th LU3R
  G2[i+L] <- LU3Rfn(T0_M, T0_F, .25, 'PV2') # PV2
  G3[i+L] <- LU3Rfn(T0_M, T0_F, .25, 'PV3') # PV3
  G4[i+L] <- LU3Rfn(T0_M, T0_F, .25, 'PV4') # PV4
  G5[i+L] <- LU3Rfn(T0_M, T0_F, .25, 'PV5') # PV5
  H1[i+L] <- LU3Rfn(T0_M, T0_F, .30, 'PV1') # PV1 reweighted 30th LU3R
  H2[i+L] <- LU3Rfn(T0_M, T0_F, .30, 'PV2') # PV2
  H3[i+L] <- LU3Rfn(T0_M, T0_F, .30, 'PV3') # PV3
  H4[i+L] <- LU3Rfn(T0_M, T0_F, .30, 'PV4') # PV4
  H5[i+L] <- LU3Rfn(T0_M, T0_F, .30, 'PV5') # PV5
  I1[i+L] <- LU3Rfn(T0_M, T0_F, .35, 'PV1') # PV1 reweighted 35th LU3R
  I2[i+L] <- LU3Rfn(T0_M, T0_F, .35, 'PV2') # PV2
  I3[i+L] <- LU3Rfn(T0_M, T0_F, .35, 'PV3') # PV3
  I4[i+L] <- LU3Rfn(T0_M, T0_F, .35, 'PV4') # PV4
  I5[i+L] <- LU3Rfn(T0_M, T0_F, .35, 'PV5') # PV5
  J1[i+L] <- LU3Rfn(T0_M, T0_F, .40, 'PV1') # PV1 reweighted 40th LU3R
  J2[i+L] <- LU3Rfn(T0_M, T0_F, .40, 'PV2') # PV2
  J3[i+L] <- LU3Rfn(T0_M, T0_F, .40, 'PV3') # PV3
  J4[i+L] <- LU3Rfn(T0_M, T0_F, .40, 'PV4') # PV4
  J5[i+L] <- LU3Rfn(T0_M, T0_F, .40, 'PV5') # PV5
  K1[i+L] <- LU3Rfn(T0_M, T0_F, .45, 'PV1') # PV1 reweighted 45th LU3R
  K2[i+L] <- LU3Rfn(T0_M, T0_F, .45, 'PV2') # PV2
  K3[i+L] <- LU3Rfn(T0_M, T0_F, .45, 'PV3') # PV3
  K4[i+L] <- LU3Rfn(T0_M, T0_F, .45, 'PV4') # PV4
  K5[i+L] <- LU3Rfn(T0_M, T0_F, .45, 'PV5') # PV5
  M1[i+L] <- LU3Rfn(T0_M, T0_F, .50, 'PV1') # PV1 reweighted 50th LU3R
  M2[i+L] <- LU3Rfn(T0_M, T0_F, .50, 'PV2') # PV2
  M3[i+L] <- LU3Rfn(T0_M, T0_F, .50, 'PV3') # PV3
  M4[i+L] <- LU3Rfn(T0_M, T0_F, .50, 'PV4') # PV4
  M5[i+L] <- LU3Rfn(T0_M, T0_F, .50, 'PV5') # PV5
  N1[i+L] <- LU3Rfn(T0_M, T0_F, .55, 'PV1') # PV1 reweighted 55th LU3R
  N2[i+L] <- LU3Rfn(T0_M, T0_F, .55, 'PV2') # PV2
  N3[i+L] <- LU3Rfn(T0_M, T0_F, .55, 'PV3') # PV3
  N4[i+L] <- LU3Rfn(T0_M, T0_F, .55, 'PV4') # PV4
  N5[i+L] <- LU3Rfn(T0_M, T0_F, .55, 'PV5') # PV5
  O1[i+L] <- LU3Rfn(T0_M, T0_F, .60, 'PV1') # PV1 reweighted 60th LU3R
  O2[i+L] <- LU3Rfn(T0_M, T0_F, .60, 'PV2') # PV2
  O3[i+L] <- LU3Rfn(T0_M, T0_F, .60, 'PV3') # PV3
  O4[i+L] <- LU3Rfn(T0_M, T0_F, .60, 'PV4') # PV4
  O5[i+L] <- LU3Rfn(T0_M, T0_F, .60, 'PV5') # PV5
  P1[i+L] <- LU3Rfn(T0_M, T0_F, .65, 'PV1') # PV1 reweighted 65th LU3R
  P2[i+L] <- LU3Rfn(T0_M, T0_F, .65, 'PV2') # PV2
  P3[i+L] <- LU3Rfn(T0_M, T0_F, .65, 'PV3') # PV3
  P4[i+L] <- LU3Rfn(T0_M, T0_F, .65, 'PV4') # PV4
  P5[i+L] <- LU3Rfn(T0_M, T0_F, .65, 'PV5') # PV5
  Q1[i+L] <- LU3Rfn(T0_M, T0_F, .70, 'PV1') # PV1 reweighted 70th LU3R
  Q2[i+L] <- LU3Rfn(T0_M, T0_F, .70, 'PV2') # PV2
  Q3[i+L] <- LU3Rfn(T0_M, T0_F, .70, 'PV3') # PV3
  Q4[i+L] <- LU3Rfn(T0_M, T0_F, .70, 'PV4') # PV4
  Q5[i+L] <- LU3Rfn(T0_M, T0_F, .70, 'PV5') # PV5
  R1[i+L] <- LU3Rfn(T0_M, T0_F, .75, 'PV1') # PV1 reweighted 75th LU3R
  R2[i+L] <- LU3Rfn(T0_M, T0_F, .75, 'PV2') # PV2
  R3[i+L] <- LU3Rfn(T0_M, T0_F, .75, 'PV3') # PV3
  R4[i+L] <- LU3Rfn(T0_M, T0_F, .75, 'PV4') # PV4
  R5[i+L] <- LU3Rfn(T0_M, T0_F, .75, 'PV5') # PV5
  S1[i+L] <- LU3Rfn(T0_M, T0_F, .80, 'PV1') # PV1 reweighted 80th LU3R
  S2[i+L] <- LU3Rfn(T0_M, T0_F, .80, 'PV2') # PV2
  S3[i+L] <- LU3Rfn(T0_M, T0_F, .80, 'PV3') # PV3
  S4[i+L] <- LU3Rfn(T0_M, T0_F, .80, 'PV4') # PV4
  S5[i+L] <- LU3Rfn(T0_M, T0_F, .80, 'PV5') # PV5
  T1[i+L] <- LU3Rfn(T0_M, T0_F, .85, 'PV1') # PV1 reweighted 85th LU3R
  T2[i+L] <- LU3Rfn(T0_M, T0_F, .85, 'PV2') # PV2
  T3[i+L] <- LU3Rfn(T0_M, T0_F, .85, 'PV3') # PV3
  T4[i+L] <- LU3Rfn(T0_M, T0_F, .85, 'PV4') # PV4
  T5[i+L] <- LU3Rfn(T0_M, T0_F, .85, 'PV5') # PV5
  V1[i+L] <- LU3Rfn(T0_M, T0_F, .90, 'PV1') # PV1 reweighted 90th LU3R
  V2[i+L] <- LU3Rfn(T0_M, T0_F, .90, 'PV2') # PV2
  V3[i+L] <- LU3Rfn(T0_M, T0_F, .90, 'PV3') # PV3
  V4[i+L] <- LU3Rfn(T0_M, T0_F, .90, 'PV4') # PV4
  V5[i+L] <- LU3Rfn(T0_M, T0_F, .90, 'PV5') # PV5
  W1[i+L] <- LU3Rfn(T0_M, T0_F, .95, 'PV1') # PV1 reweighted 95th LU3R
  W2[i+L] <- LU3Rfn(T0_M, T0_F, .95, 'PV2') # PV2
  W3[i+L] <- LU3Rfn(T0_M, T0_F, .95, 'PV3') # PV3
  W4[i+L] <- LU3Rfn(T0_M, T0_F, .95, 'PV4') # PV4
  W5[i+L] <- LU3Rfn(T0_M, T0_F, .95, 'PV5') # PV5
}

# jackknife sampling variance
JSV_05 <- mean(c(sum((B1-LU3R05_1)^2), # 5th LU3R
                 sum((B2-LU3R05_2)^2),
                 sum((B3-LU3R05_3)^2),
                 sum((B4-LU3R05_4)^2),
                 sum((B5-LU3R05_5)^2)))/2
JSV_10 <- mean(c(sum((C1-LU3R10_1)^2), # 10th LU3R
                 sum((C2-LU3R10_2)^2),
                 sum((C3-LU3R10_3)^2),
                 sum((C4-LU3R10_4)^2),
                 sum((C5-LU3R10_5)^2)))/2
JSV_15 <- mean(c(sum((E1-LU3R15_1)^2), # 15th LU3R
                 sum((E2-LU3R15_2)^2),
                 sum((E3-LU3R15_3)^2),
                 sum((E4-LU3R15_4)^2),
                 sum((E5-LU3R15_5)^2)))/2
JSV_20 <- mean(c(sum((F1-LU3R20_1)^2), # 20th LU3R
                 sum((F2-LU3R20_2)^2),
                 sum((F3-LU3R20_3)^2),
                 sum((F4-LU3R20_4)^2),
                 sum((F5-LU3R20_5)^2)))/2
JSV_25 <- mean(c(sum((G1-LU3R25_1)^2), # 25th LU3R
                 sum((G2-LU3R25_2)^2),
                 sum((G3-LU3R25_3)^2),
                 sum((G4-LU3R25_4)^2),
                 sum((G5-LU3R25_5)^2)))/2
JSV_30 <- mean(c(sum((H1-LU3R30_1)^2), # 30th LU3R
                 sum((H2-LU3R30_2)^2),
                 sum((H3-LU3R30_3)^2),
                 sum((H4-LU3R30_4)^2),
                 sum((H5-LU3R30_5)^2)))/2
JSV_35 <- mean(c(sum((I1-LU3R35_1)^2), # 35th LU3R
                 sum((I2-LU3R35_2)^2),
                 sum((I3-LU3R35_3)^2),
                 sum((I4-LU3R35_4)^2),
                 sum((I5-LU3R35_5)^2)))/2
JSV_40 <- mean(c(sum((J1-LU3R40_1)^2), # 40th LU3R
                 sum((J2-LU3R40_2)^2),
                 sum((J3-LU3R40_3)^2),
                 sum((J4-LU3R40_4)^2),
                 sum((J5-LU3R40_5)^2)))/2
JSV_45 <- mean(c(sum((K1-LU3R45_1)^2), # 45th LU3R
                 sum((K2-LU3R45_2)^2),
                 sum((K3-LU3R45_3)^2),
                 sum((K4-LU3R45_4)^2),
                 sum((K5-LU3R45_5)^2)))/2
JSV_50 <- mean(c(sum((M1-LU3R50_1)^2), # 50th LU3R
                 sum((M2-LU3R50_2)^2),
                 sum((M3-LU3R50_3)^2),
                 sum((M4-LU3R50_4)^2),
                 sum((M5-LU3R50_5)^2)))/2
JSV_55 <- mean(c(sum((N1-LU3R55_1)^2), # 55th LU3R
                 sum((N2-LU3R55_2)^2),
                 sum((N3-LU3R55_3)^2),
                 sum((N4-LU3R55_4)^2),
                 sum((N5-LU3R55_5)^2)))/2
JSV_60 <- mean(c(sum((O1-LU3R60_1)^2), # 60th LU3R
                 sum((O2-LU3R60_2)^2),
                 sum((O3-LU3R60_3)^2),
                 sum((O4-LU3R60_4)^2),
                 sum((O5-LU3R60_5)^2)))/2
JSV_65 <- mean(c(sum((P1-LU3R65_1)^2), # 65th LU3R
                 sum((P2-LU3R65_2)^2),
                 sum((P3-LU3R65_3)^2),
                 sum((P4-LU3R65_4)^2),
                 sum((P5-LU3R65_5)^2)))/2
JSV_70 <- mean(c(sum((Q1-LU3R70_1)^2), # 70th LU3R
                 sum((Q2-LU3R70_2)^2),
                 sum((Q3-LU3R70_3)^2),
                 sum((Q4-LU3R70_4)^2),
                 sum((Q5-LU3R70_5)^2)))/2
JSV_75 <- mean(c(sum((R1-LU3R75_1)^2), # 75th LU3R
                 sum((R2-LU3R75_2)^2),
                 sum((R3-LU3R75_3)^2),
                 sum((R4-LU3R75_4)^2),
                 sum((R5-LU3R75_5)^2)))/2
JSV_80 <- mean(c(sum((S1-LU3R80_1)^2), # 80th LU3R
                 sum((S2-LU3R80_2)^2),
                 sum((S3-LU3R80_3)^2),
                 sum((S4-LU3R80_4)^2),
                 sum((S5-LU3R80_5)^2)))/2
JSV_85 <- mean(c(sum((T1-LU3R85_1)^2), # 85th LU3R
                 sum((T2-LU3R85_2)^2),
                 sum((T3-LU3R85_3)^2),
                 sum((T4-LU3R85_4)^2),
                 sum((T5-LU3R85_5)^2)))/2
JSV_90 <- mean(c(sum((V1-LU3R90_1)^2), # 90th LU3R
                 sum((V2-LU3R90_2)^2),
                 sum((V3-LU3R90_3)^2),
                 sum((V4-LU3R90_4)^2),
                 sum((V5-LU3R90_5)^2)))/2
JSV_95 <- mean(c(sum((W1-LU3R95_1)^2), # 95th LU3R
                 sum((W2-LU3R95_2)^2),
                 sum((W3-LU3R95_3)^2),
                 sum((W4-LU3R95_4)^2),
                 sum((W5-LU3R95_5)^2)))/2
JSV_Md95U <- mean(c(sum((W1-M1-Md95U_1)^2), # Md95U
                    sum((W2-M2-Md95U_2)^2),
                    sum((W3-M3-Md95U_3)^2),
                    sum((W4-M4-Md95U_4)^2),
                    sum((W5-M5-Md95U_5)^2)))/2
JSV_Md90U <- mean(c(sum((V1-M1-Md90U_1)^2), # Md90U
                    sum((V2-M2-Md90U_2)^2),
                    sum((V3-M3-Md90U_3)^2),
                    sum((V4-M4-Md90U_4)^2),
                    sum((V5-M5-Md90U_5)^2)))/2
JSV_Md10U <- mean(c(sum((M1-C1-Md10U_1)^2), # Md10U
                    sum((M2-C2-Md10U_2)^2),
                    sum((M3-C3-Md10U_3)^2),
                    sum((M4-C4-Md10U_4)^2),
                    sum((M5-C5-Md10U_5)^2)))/2
JSV_Md05U <- mean(c(sum((M1-B1-Md05U_1)^2), # Md05U
                    sum((M2-B2-Md05U_2)^2),
                    sum((M3-B3-Md05U_3)^2),
                    sum((M4-B4-Md05U_4)^2),
                    sum((M5-B5-Md05U_5)^2)))/2

# imputation variance
IV_05 <- .3*sum((LU3R05s-LU3R05)^2) # 5th LU3R
IV_10 <- .3*sum((LU3R10s-LU3R10)^2) # 10th LU3R
IV_15 <- .3*sum((LU3R15s-LU3R15)^2) # 15th LU3R
IV_20 <- .3*sum((LU3R20s-LU3R20)^2) # 20th LU3R
IV_25 <- .3*sum((LU3R25s-LU3R25)^2) # 25th LU3R
IV_30 <- .3*sum((LU3R30s-LU3R30)^2) # 30th LU3R
IV_35 <- .3*sum((LU3R35s-LU3R35)^2) # 35th LU3R
IV_40 <- .3*sum((LU3R40s-LU3R40)^2) # 40th LU3R
IV_45 <- .3*sum((LU3R45s-LU3R45)^2) # 45th LU3R
IV_50 <- .3*sum((LU3R50s-LU3R50)^2) # 50th LU3R
IV_55 <- .3*sum((LU3R55s-LU3R55)^2) # 55th LU3R
IV_60 <- .3*sum((LU3R60s-LU3R60)^2) # 60th LU3R
IV_65 <- .3*sum((LU3R65s-LU3R65)^2) # 65th LU3R
IV_70 <- .3*sum((LU3R70s-LU3R70)^2) # 70th LU3R
IV_75 <- .3*sum((LU3R75s-LU3R75)^2) # 75th LU3R
IV_80 <- .3*sum((LU3R80s-LU3R80)^2) # 80th LU3R
IV_85 <- .3*sum((LU3R85s-LU3R85)^2) # 85th LU3R
IV_90 <- .3*sum((LU3R90s-LU3R90)^2) # 90th LU3R
IV_95 <- .3*sum((LU3R95s-LU3R95)^2) # 95th LU3R
IV_Md95U <- .3*sum((Md95Us-Md95U)^2) # Md95U
IV_Md90U <- .3*sum((Md90Us-Md90U)^2) # Md90U
IV_Md10U <- .3*sum((Md10Us-Md10U)^2) # Md10U
IV_Md05U <- .3*sum((Md05Us-Md05U)^2) # Md05U


# total variance = jackknife + imputation
TV_05U <- JSV_05 + IV_05 # total variance of 5th percentile LU3R
SE_05U <- sqrt(TV_05U) # standard error of 5th percentile LU3R

TV_10U <- JSV_10 + IV_10 # total variance of 10th percentile LU3R
SE_10U <- sqrt(TV_10U) # standard error of 10th percentile LU3R

TV_15U <- JSV_15 + IV_15 # total variance of 15th percentile LU3R
SE_15U <- sqrt(TV_15U) # standard error of 15th percentile LU3R

TV_20U <- JSV_20 + IV_20 # total variance of 20th percentile LU3R
SE_20U <- sqrt(TV_20U) # standard error of 20th percentile LU3R

TV_25U <- JSV_25 + IV_25 # total variance of 25th percentile LU3R
SE_25U <- sqrt(TV_25U) # standard error of 25th percentile LU3R

TV_30U <- JSV_30 + IV_30 # total variance of 30th percentile LU3R
SE_30U <- sqrt(TV_30U) # standard error of 30th percentile LU3R

TV_35U <- JSV_35 + IV_35 # total variance of 35th percentile LU3R
SE_35U <- sqrt(TV_35U) # standard error of 35th percentile LU3R

TV_40U <- JSV_40 + IV_40 # total variance of 40th percentile LU3R
SE_40U <- sqrt(TV_40U) # standard error of 40th percentile LU3R

TV_45U <- JSV_45 + IV_45 # total variance of 45th percentile LU3R
SE_45U <- sqrt(TV_45U) # standard error of 45th percentile LU3R

TV_50U <- JSV_50 + IV_50 # total variance of 50th percentile LU3R
SE_50U <- sqrt(TV_50U) # standard error of 50th percentile LU3R

TV_55U <- JSV_55 + IV_55 # total variance of 55th percentile LU3R
SE_55U <- sqrt(TV_55U) # standard error of 55th percentile LU3R

TV_60U <- JSV_60 + IV_60 # total variance of 60th percentile LU3R
SE_60U <- sqrt(TV_60U) # standard error of 60th percentile LU3R

TV_65U <- JSV_65 + IV_65 # total variance of 65th percentile LU3R
SE_65U <- sqrt(TV_65U) # standard error of 65th percentile LU3R

TV_70U <- JSV_70 + IV_70 # total variance of 70th percentile LU3R
SE_70U <- sqrt(TV_70U) # standard error of 70th percentile LU3R

TV_75U <- JSV_75 + IV_75 # total variance of 75th percentile LU3R
SE_75U <- sqrt(TV_75U) # standard error of 75th percentile LU3R

TV_80U <- JSV_80 + IV_80 # total variance of 80th percentile LU3R
SE_80U <- sqrt(TV_80U) # standard error of 80th percentile LU3R

TV_85U <- JSV_85 + IV_85 # total variance of 85th percentile LU3R
SE_85U <- sqrt(TV_85U) # standard error of 85th percentile LU3R

TV_90U <- JSV_90 + IV_90 # total variance of 90th percentile LU3R
SE_90U <- sqrt(TV_90U) # standard error of 90th percentile LU3R

TV_95U <- JSV_95 + IV_95 # total variance of 95th percentile LU3R
SE_95U <- sqrt(TV_95U) # standard error of 95th percentile LU3R

TV_Md95U <- JSV_Md95U + IV_Md95U # total variance of Md95U
SE_Md95U <- sqrt(TV_Md95U) # standard error of Md95U

TV_Md90U <- JSV_Md90U + IV_Md90U # total variance of Md90U
SE_Md90U <- sqrt(TV_Md90U) # standard error of Md90U

TV_Md10U <- JSV_Md10U + IV_Md10U # total variance of Md10U
SE_Md10U <- sqrt(TV_Md10U) # standard error of Md10U

TV_Md05U <- JSV_Md05U + IV_Md05U # total variance of Md05U
SE_Md05U <- sqrt(TV_Md05U) # standard error of Md05U



#### SEs: d, U3, LVR, LVR_L, LVR_R, LMADR, LMADR_L, LMADR_R, LGMDR ####

J1 <- J2 <- J3 <- J4 <- J5 <- numeric(2*L) # empty containers
A1 <- A2 <- A3 <- A4 <- A5 <- numeric(2*L) # empty containers
C1 <- C2 <- C3 <- C4 <- C5 <- numeric(2*L) # empty containers
K1 <- K2 <- K3 <- K4 <- K5 <- numeric(2*L) # empty containers
N1 <- N2 <- N3 <- N4 <- N5 <- numeric(2*L) # empty containers
I1 <- I2 <- I3 <- I4 <- I5 <- numeric(2*L) # empty containers
F1 <- F2 <- F3 <- F4 <- F5 <- numeric(2*L) # empty containers
E1 <- E2 <- E3 <- E4 <- E5 <- numeric(2*L) # empty containers
G1 <- G2 <- G3 <- G4 <- G5 <- numeric(2*L) # empty containers

# perform jackknife resampling of d, U3, LVR, LVR_L, LVR_R, LMADR, LMADR_L, LMADR_R
for (i in 1:L) { # for each JK zone
  T0 <- T15 # create/restore duplicate
  # double weights if JK code is 1 and zero weights if JK code is 0
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 1),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  J1[i] <- dfn(T0_M, T0_F, 'PV1') # PV1 reweighted d
  J2[i] <- dfn(T0_M, T0_F, 'PV2') # PV2
  J3[i] <- dfn(T0_M, T0_F, 'PV3') # PV3
  J4[i] <- dfn(T0_M, T0_F, 'PV4') # PV4
  J5[i] <- dfn(T0_M, T0_F, 'PV5') # PV5
  A1[i] <- U3fn(T0_M, T0_F, 'PV1') # PV1 reweighted U3
  A2[i] <- U3fn(T0_M, T0_F, 'PV2') # PV2
  A3[i] <- U3fn(T0_M, T0_F, 'PV3') # PV3
  A4[i] <- U3fn(T0_M, T0_F, 'PV4') # PV4
  A5[i] <- U3fn(T0_M, T0_F, 'PV5') # PV5
  C1[i] <- LVRfn(T0_M, T0_F, 'PV1') # PV1 reweighted LVR
  C2[i] <- LVRfn(T0_M, T0_F, 'PV2') # PV2
  C3[i] <- LVRfn(T0_M, T0_F, 'PV3') # PV3
  C4[i] <- LVRfn(T0_M, T0_F, 'PV4') # PV4
  C5[i] <- LVRfn(T0_M, T0_F, 'PV5') # PV5
  K1[i] <- LVR_Tfn(T0_M, T0_F, 'PV1', 'L') # PV1 reweighted LVR_L
  K2[i] <- LVR_Tfn(T0_M, T0_F, 'PV2', 'L') # PV2
  K3[i] <- LVR_Tfn(T0_M, T0_F, 'PV3', 'L') # PV3
  K4[i] <- LVR_Tfn(T0_M, T0_F, 'PV4', 'L') # PV4
  K5[i] <- LVR_Tfn(T0_M, T0_F, 'PV5', 'L') # PV5
  N1[i] <- LVR_Tfn(T0_M, T0_F, 'PV1', 'R') # PV1 reweighted LVR_R
  N2[i] <- LVR_Tfn(T0_M, T0_F, 'PV2', 'R') # PV2
  N3[i] <- LVR_Tfn(T0_M, T0_F, 'PV3', 'R') # PV3
  N4[i] <- LVR_Tfn(T0_M, T0_F, 'PV4', 'R') # PV4
  N5[i] <- LVR_Tfn(T0_M, T0_F, 'PV5', 'R') # PV5
  I1[i] <- LMADRfn(T0_M, T0_F, 'PV1') # PV1 reweighted LMADR
  I2[i] <- LMADRfn(T0_M, T0_F, 'PV2') # PV2
  I3[i] <- LMADRfn(T0_M, T0_F, 'PV3') # PV3
  I4[i] <- LMADRfn(T0_M, T0_F, 'PV4') # PV4
  I5[i] <- LMADRfn(T0_M, T0_F, 'PV5') # PV5
  F1[i] <- LMADR_Tfn(T0_M, T0_F, 'PV1', 'L') # PV1 reweighted LMADR_L
  F2[i] <- LMADR_Tfn(T0_M, T0_F, 'PV2', 'L') # PV2
  F3[i] <- LMADR_Tfn(T0_M, T0_F, 'PV3', 'L') # PV3
  F4[i] <- LMADR_Tfn(T0_M, T0_F, 'PV4', 'L') # PV4
  F5[i] <- LMADR_Tfn(T0_M, T0_F, 'PV5', 'L') # PV5
  E1[i] <- LMADR_Tfn(T0_M, T0_F, 'PV1', 'R') # PV1 reweighted LMADR_R
  E2[i] <- LMADR_Tfn(T0_M, T0_F, 'PV2', 'R') # PV2
  E3[i] <- LMADR_Tfn(T0_M, T0_F, 'PV3', 'R') # PV3
  E4[i] <- LMADR_Tfn(T0_M, T0_F, 'PV4', 'R') # PV4
  E5[i] <- LMADR_Tfn(T0_M, T0_F, 'PV5', 'R') # PV5
  G1[i] <- LGMDRfn(T0_M, T0_F, 'PV1') # PV1 reweighted LGMDR
  G2[i] <- LGMDRfn(T0_M, T0_F, 'PV2') # PV2
  G3[i] <- LGMDRfn(T0_M, T0_F, 'PV3') # PV3
  G4[i] <- LGMDRfn(T0_M, T0_F, 'PV4') # PV4
  G5[i] <- LGMDRfn(T0_M, T0_F, 'PV5') # PV5
  
  T0 <- T15 # restore duplicate
  # double weights if JK code is 0 and zero weights if JK code is 1
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 0),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  J1[i+L] <- dfn(T0_M, T0_F, 'PV1') # PV1 reweighted d
  J2[i+L] <- dfn(T0_M, T0_F, 'PV2') # PV2
  J3[i+L] <- dfn(T0_M, T0_F, 'PV3') # PV3
  J4[i+L] <- dfn(T0_M, T0_F, 'PV4') # PV4
  J5[i+L] <- dfn(T0_M, T0_F, 'PV5') # PV5
  A1[i+L] <- U3fn(T0_M, T0_F, 'PV1') # PV1 reweighted U3
  A2[i+L] <- U3fn(T0_M, T0_F, 'PV2') # PV2
  A3[i+L] <- U3fn(T0_M, T0_F, 'PV3') # PV3
  A4[i+L] <- U3fn(T0_M, T0_F, 'PV4') # PV4
  A5[i+L] <- U3fn(T0_M, T0_F, 'PV5') # PV5
  C1[i+L] <- LVRfn(T0_M, T0_F, 'PV1') # PV1 reweighted LVR
  C2[i+L] <- LVRfn(T0_M, T0_F, 'PV2') # PV2
  C3[i+L] <- LVRfn(T0_M, T0_F, 'PV3') # PV3
  C4[i+L] <- LVRfn(T0_M, T0_F, 'PV4') # PV4
  C5[i+L] <- LVRfn(T0_M, T0_F, 'PV5') # PV5
  K1[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV1', 'L') # PV1 reweighted LVR_L
  K2[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV2', 'L') # PV2
  K3[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV3', 'L') # PV3
  K4[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV4', 'L') # PV4
  K5[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV5', 'L') # PV5
  N1[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV1', 'R') # PV1 reweighted LVR_R
  N2[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV2', 'R') # PV2
  N3[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV3', 'R') # PV3
  N4[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV4', 'R') # PV4
  N5[i+L] <- LVR_Tfn(T0_M, T0_F, 'PV5', 'R') # PV5
  I1[i+L] <- LMADRfn(T0_M, T0_F, 'PV1') # PV1 reweighted LMADR
  I2[i+L] <- LMADRfn(T0_M, T0_F, 'PV2') # PV2
  I3[i+L] <- LMADRfn(T0_M, T0_F, 'PV3') # PV3
  I4[i+L] <- LMADRfn(T0_M, T0_F, 'PV4') # PV4
  I5[i+L] <- LMADRfn(T0_M, T0_F, 'PV5') # PV5
  F1[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV1', 'L') # PV1 reweighted LMADR_L
  F2[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV2', 'L') # PV2
  F3[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV3', 'L') # PV3
  F4[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV4', 'L') # PV4
  F5[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV5', 'L') # PV5
  E1[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV1', 'R') # PV1 reweighted LMADR_R
  E2[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV2', 'R') # PV2
  E3[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV3', 'R') # PV3
  E4[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV4', 'R') # PV4
  E5[i+L] <- LMADR_Tfn(T0_M, T0_F, 'PV5', 'R') # PV5
  G1[i+L] <- LGMDRfn(T0_M, T0_F, 'PV1') # PV1 reweighted LGMDR
  G2[i+L] <- LGMDRfn(T0_M, T0_F, 'PV2') # PV2
  G3[i+L] <- LGMDRfn(T0_M, T0_F, 'PV3') # PV3
  G4[i+L] <- LGMDRfn(T0_M, T0_F, 'PV4') # PV4
  G5[i+L] <- LGMDRfn(T0_M, T0_F, 'PV5') # PV5
}

# jackknife sampling variance
JSV_d <- mean(c(sum((J1-d1)^2), # d
                sum((J2-d2)^2),
                sum((J3-d3)^2),
                sum((J4-d4)^2),
                sum((J5-d5)^2)))/2
JSV_U3 <- mean(c(sum((A1-U3_1)^2), # U3
                 sum((A2-U3_2)^2),
                 sum((A3-U3_3)^2),
                 sum((A4-U3_4)^2),
                 sum((A5-U3_5)^2)))/2
JSV_LVR <- mean(c(sum((C1-LVR1)^2), # LVR
                  sum((C2-LVR2)^2),
                  sum((C3-LVR3)^2),
                  sum((C4-LVR4)^2),
                  sum((C5-LVR5)^2)))/2
JSV_LVR_L <- mean(c(sum((K1-LVR_L1)^2), # LVR_L
                    sum((K2-LVR_L2)^2),
                    sum((K3-LVR_L3)^2),
                    sum((K4-LVR_L4)^2),
                    sum((K5-LVR_L5)^2)))/2
JSV_LVR_R <- mean(c(sum((N1-LVR_R1)^2), # LVR_R
                    sum((N2-LVR_R2)^2),
                    sum((N3-LVR_R3)^2),
                    sum((N4-LVR_R4)^2),
                    sum((N5-LVR_R5)^2)))/2
JSV_LMADR <- mean(c(sum((I1-LMADR1)^2), # LMADR
                    sum((I2-LMADR2)^2),
                    sum((I3-LMADR3)^2),
                    sum((I4-LMADR4)^2),
                    sum((I5-LMADR5)^2)))/2
JSV_LMADR_L <- mean(c(sum((F1-LMADR_L1)^2), # LMADR_L
                      sum((F2-LMADR_L2)^2),
                      sum((F3-LMADR_L3)^2),
                      sum((F4-LMADR_L4)^2),
                      sum((F5-LMADR_L5)^2)))/2
JSV_LMADR_R <- mean(c(sum((E1-LMADR_R1)^2), # LMADR_R
                      sum((E2-LMADR_R2)^2),
                      sum((E3-LMADR_R3)^2),
                      sum((E4-LMADR_R4)^2),
                      sum((E5-LMADR_R5)^2)))/2
JSV_LGMDR <- mean(c(sum((G1-LGMDR1)^2), # LGMDR
                    sum((G2-LGMDR2)^2),
                    sum((G3-LGMDR3)^2),
                    sum((G4-LGMDR4)^2),
                    sum((G5-LGMDR5)^2)))/2

# imputation variance
IV_d <- .3*sum((ds-d)^2) # d
IV_U3 <- .3*sum((U3s-U3)^2) # U3
IV_LVR <- .3*sum((LVRs-LVR)^2) # LVR
IV_LVR_L <- .3*sum((LVR_Ls-LVR_L)^2) # LVR_L
IV_LVR_R <- .3*sum((LVR_Rs-LVR_R)^2) # LVR_R
IV_LMADR <- .3*sum((LMADRs-LMADR)^2) # LMADR
IV_LMADR_L <- .3*sum((LMADR_Ls-LMADR_L)^2) # LMADR_L
IV_LMADR_R <- .3*sum((LMADR_Rs-LMADR_R)^2) # LMADR_R
IV_LGMDR <- .3*sum((LGMDRs-LGMDR)^2) # LGMDR


# total variance = jackknife + imputation
TV_d <- JSV_d + IV_d # total variance of d
SE_d <- sqrt(TV_d) # standard error of d

TV_U3 <- JSV_U3 + IV_U3 # total variance of U3
SE_U3 <- sqrt(TV_U3) # standard error of U3

TV_LVR <- JSV_LVR + IV_LVR # total variance of LVR
SE_LVR <- sqrt(TV_LVR) # standard error of LVR

TV_LVR_L <- JSV_LVR_L + IV_LVR_L # total variance of LVR_L
SE_LVR_L <- sqrt(TV_LVR_L) # standard error of LVR_L

TV_LVR_R <- JSV_LVR_R + IV_LVR_R # total variance of LVR_R
SE_LVR_R <- sqrt(TV_LVR_R) # standard error of LVR_R

TV_LMADR <- JSV_LMADR + IV_LMADR # total variance of LMADR
SE_LMADR <- sqrt(TV_LMADR) # standard error of LMADR

TV_LMADR_L <- JSV_LMADR_L + IV_LMADR_L # total variance of LMADR_L
SE_LMADR_L <- sqrt(TV_LMADR_L) # standard error of LMADR_L

TV_LMADR_R <- JSV_LMADR_R + IV_LMADR_R # total variance of LMADR_R
SE_LMADR_R <- sqrt(TV_LMADR_R) # standard error of LMADR_R

TV_LGMDR <- JSV_LGMDR + IV_LGMDR # total variance of LGMDR
SE_LGMDR <- sqrt(TV_LGMDR) # standard error of LGMDR



#### SEs: Probability of superiority ####

J1 <- J2 <- J3 <- J4 <- J5 <- numeric(2*L) # empty containers

# perform jackknife resampling of PSs
for (i in 1:L) { # for each JK zone
  T0 <- T15 # create/restore duplicate
  # double weights if JK code is 1 and zero weights if JK code is 0
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 1),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  J1[i] <- PSfn(T0_M, T0_F, 'PV1') # PV1 reweighted PS
  J2[i] <- PSfn(T0_M, T0_F, 'PV2') # PV2
  J3[i] <- PSfn(T0_M, T0_F, 'PV3') # PV3
  J4[i] <- PSfn(T0_M, T0_F, 'PV4') # PV4
  J5[i] <- PSfn(T0_M, T0_F, 'PV5') # PV5
  
  T0 <- T15 # restore duplicate
  # double weights if JK code is 0 and zero weights if JK code is 1
  T0[which(T0$JKZ == i & T0$JKR == 0),'HWt'] <- 2*T0[which(T0$JKZ == i & T0$JKR == 0),'HWt']
  T0[which(T0$JKZ == i & T0$JKR == 1),'HWt'] <- 0
  T0_F <- T0[which(T0$Sex == 1),] # female subset
  T0_M <- T0[which(T0$Sex == 2),] # male subset
  J1[i+L] <- PSfn(T0_M, T0_F, 'PV1') # PV1 reweighted PS
  J2[i+L] <- PSfn(T0_M, T0_F, 'PV2') # PV2
  J3[i+L] <- PSfn(T0_M, T0_F, 'PV3') # PV3
  J4[i+L] <- PSfn(T0_M, T0_F, 'PV4') # PV4
  J5[i+L] <- PSfn(T0_M, T0_F, 'PV5') # PV5
  
  print(paste0(i, '/', L, ' at ', Sys.time()), quote = F) # print updates
}

# jackknife sampling variance
JSV_PS <- mean(c(sum((J1-PS1)^2),
                 sum((J2-PS2)^2),
                 sum((J3-PS3)^2),
                 sum((J4-PS4)^2),
                 sum((J5-PS5)^2)))/2

# imputation variance
IV_PS <- .3*sum((PSs-PS)^2)


# total variance = jackknife + imputation
TV_PS <- JSV_PS + IV_PS # total variance of PS
SE_PS <- sqrt(TV_PS) # standard error of PS



####################################
##### 95% Confidence Intervals #####
####################################

# 95% CIs computed under the assumption of normality
CI95 <- qnorm(.975)


### Means and Medians
Mn_lo <- Mn-SE_MnT*CI95
Mn_up <- Mn+SE_MnT*CI95

Md_lo <- Md-SE_MdT*CI95
Md_up <- Md+SE_MdT*CI95

Mn_F_lo <- Mn_F-SE_MnF*CI95
Mn_F_up <- Mn_F+SE_MnF*CI95

Md_F_lo <- Md_F-SE_MdF*CI95
Md_F_up <- Md_F+SE_MdF*CI95

Mn_M_lo <- Mn_M-SE_MnM*CI95
Mn_M_up <- Mn_M+SE_MnM*CI95

Md_M_lo <- Md_M-SE_MdM*CI95
Md_M_up <- Md_M+SE_MdM*CI95

MnDf_lo <- MnDf-SE_MnDf*CI95
MnDf_up <- MnDf+SE_MnDf*CI95

MdDf_lo <- MdDf-SE_MdDf*CI95
MdDf_up <- MdDf+SE_MdDf*CI95


### TPRs
TPRMn_lo <- exp(LTPRMn-SE_MnT*CI95)
TPRMn_up <- exp(LTPRMn+SE_MnT*CI95)

TPR05_lo <- exp(LTPR05-SE_05T*CI95)
TPR05_up <- exp(LTPR05+SE_05T*CI95)

TPR10_lo <- exp(LTPR10-SE_10T*CI95)
TPR10_up <- exp(LTPR10+SE_10T*CI95)

TPR15_lo <- exp(LTPR15-SE_15T*CI95)
TPR15_up <- exp(LTPR15+SE_15T*CI95)

TPR20_lo <- exp(LTPR20-SE_20T*CI95)
TPR20_up <- exp(LTPR20+SE_20T*CI95)

TPR25_lo <- exp(LTPR25-SE_25T*CI95)
TPR25_up <- exp(LTPR25+SE_25T*CI95)

TPR30_lo <- exp(LTPR30-SE_30T*CI95)
TPR30_up <- exp(LTPR30+SE_30T*CI95)

TPR35_lo <- exp(LTPR35-SE_35T*CI95)
TPR35_up <- exp(LTPR35+SE_35T*CI95)

TPR40_lo <- exp(LTPR40-SE_40T*CI95)
TPR40_up <- exp(LTPR40+SE_40T*CI95)

TPR45_lo <- exp(LTPR45-SE_45T*CI95)
TPR45_up <- exp(LTPR45+SE_45T*CI95)

TPR50_lo <- exp(LTPR50-SE_50T*CI95)
TPR50_up <- exp(LTPR50+SE_50T*CI95)

TPR55_lo <- exp(LTPR55-SE_55T*CI95)
TPR55_up <- exp(LTPR55+SE_55T*CI95)

TPR60_lo <- exp(LTPR60-SE_60T*CI95)
TPR60_up <- exp(LTPR60+SE_60T*CI95)

TPR65_lo <- exp(LTPR65-SE_65T*CI95)
TPR65_up <- exp(LTPR65+SE_65T*CI95)

TPR70_lo <- exp(LTPR70-SE_70T*CI95)
TPR70_up <- exp(LTPR70+SE_70T*CI95)

TPR75_lo <- exp(LTPR75-SE_75T*CI95)
TPR75_up <- exp(LTPR75+SE_75T*CI95)

TPR80_lo <- exp(LTPR80-SE_80T*CI95)
TPR80_up <- exp(LTPR80+SE_80T*CI95)

TPR85_lo <- exp(LTPR85-SE_85T*CI95)
TPR85_up <- exp(LTPR85+SE_85T*CI95)

TPR90_lo <- exp(LTPR90-SE_90T*CI95)
TPR90_up <- exp(LTPR90+SE_90T*CI95)

TPR95_lo <- exp(LTPR95-SE_95T*CI95)
TPR95_up <- exp(LTPR95+SE_95T*CI95)

Md95T_lo <- Md95T-SE_Md95T*CI95
Md95T_up <- Md95T+SE_Md95T*CI95
Md90T_lo <- Md90T-SE_Md90T*CI95
Md90T_up <- Md90T+SE_Md90T*CI95
Md10T_lo <- Md10T-SE_Md10T*CI95
Md10T_up <- Md10T+SE_Md10T*CI95
Md05T_lo <- Md05T-SE_Md05T*CI95
Md05T_up <- Md05T+SE_Md05T*CI95

Mn95T_lo <- Mn95T-SE_Mn95T*CI95
Mn95T_up <- Mn95T+SE_Mn95T*CI95
Mn90T_lo <- Mn90T-SE_Mn90T*CI95
Mn90T_up <- Mn90T+SE_Mn90T*CI95
Mn10T_lo <- Mn10T-SE_Mn10T*CI95
Mn10T_up <- Mn10T+SE_Mn10T*CI95
Mn05T_lo <- Mn05T-SE_Mn05T*CI95
Mn05T_up <- Mn05T+SE_Mn05T*CI95


### U3Rs
U3R05_lo <- exp(LU3R05-SE_05U*CI95)
U3R05_up <- exp(LU3R05+SE_05U*CI95)

U3R10_lo <- exp(LU3R10-SE_10U*CI95)
U3R10_up <- exp(LU3R10+SE_10U*CI95)

U3R15_lo <- exp(LU3R15-SE_15U*CI95)
U3R15_up <- exp(LU3R15+SE_15U*CI95)

U3R20_lo <- exp(LU3R20-SE_20U*CI95)
U3R20_up <- exp(LU3R20+SE_20U*CI95)

U3R25_lo <- exp(LU3R25-SE_25U*CI95)
U3R25_up <- exp(LU3R25+SE_25U*CI95)

U3R30_lo <- exp(LU3R30-SE_30U*CI95)
U3R30_up <- exp(LU3R30+SE_30U*CI95)

U3R35_lo <- exp(LU3R35-SE_35U*CI95)
U3R35_up <- exp(LU3R35+SE_35U*CI95)

U3R40_lo <- exp(LU3R40-SE_40U*CI95)
U3R40_up <- exp(LU3R40+SE_40U*CI95)

U3R45_lo <- exp(LU3R45-SE_45U*CI95)
U3R45_up <- exp(LU3R45+SE_45U*CI95)

U3R50_lo <- exp(LU3R50-SE_50U*CI95)
U3R50_up <- exp(LU3R50+SE_50U*CI95)

U3R55_lo <- exp(LU3R55-SE_55U*CI95)
U3R55_up <- exp(LU3R55+SE_55U*CI95)

U3R60_lo <- exp(LU3R60-SE_60U*CI95)
U3R60_up <- exp(LU3R60+SE_60U*CI95)

U3R65_lo <- exp(LU3R65-SE_65U*CI95)
U3R65_up <- exp(LU3R65+SE_65U*CI95)

U3R70_lo <- exp(LU3R70-SE_70U*CI95)
U3R70_up <- exp(LU3R70+SE_70U*CI95)

U3R75_lo <- exp(LU3R75-SE_75U*CI95)
U3R75_up <- exp(LU3R75+SE_75U*CI95)

U3R80_lo <- exp(LU3R80-SE_80U*CI95)
U3R80_up <- exp(LU3R80+SE_80U*CI95)

U3R85_lo <- exp(LU3R85-SE_85U*CI95)
U3R85_up <- exp(LU3R85+SE_85U*CI95)

U3R90_lo <- exp(LU3R90-SE_90U*CI95)
U3R90_up <- exp(LU3R90+SE_90U*CI95)

U3R95_lo <- exp(LU3R95-SE_95U*CI95)
U3R95_up <- exp(LU3R95+SE_95U*CI95)

Md95U_lo <- Md95U-SE_Md95U*CI95
Md95U_up <- Md95U+SE_Md95U*CI95
Md90U_lo <- Md90U-SE_Md90U*CI95
Md90U_up <- Md90U+SE_Md90U*CI95
Md10U_lo <- Md10U-SE_Md10U*CI95
Md10U_up <- Md10U+SE_Md10U*CI95
Md05U_lo <- Md05U-SE_Md05U*CI95
Md05U_up <- Md05U+SE_Md05U*CI95



### Other effect sizes

# Cohen's d
d_lo <- d-SE_d*CI95
d_up <- d+SE_d*CI95

# U3
U3_lo <- U3-SE_U3*CI95
U3_up <- U3+SE_U3*CI95

# Probability of superiority
PS_lo <- PS-SE_PS*CI95
PS_up <- PS+SE_PS*CI95

# Variance ratio
VR_lo <- exp(LVR-SE_LVR*CI95)
VR_up <- exp(LVR+SE_LVR*CI95)

# Variance ratio (left tail)
VR_L_lo <- exp(LVR_L-SE_LVR_L*CI95)
VR_L_up <- exp(LVR_L+SE_LVR_L*CI95)

# Variance ratio (right tail)
VR_R_lo <- exp(LVR_R-SE_LVR_R*CI95)
VR_R_up <- exp(LVR_R+SE_LVR_R*CI95)

# Mean absolute deviation ratio
MADR_lo <- exp(LMADR-SE_LMADR*CI95)
MADR_up <- exp(LMADR+SE_LMADR*CI95)

# Mean absolute deviation ratio (left tail)
MADR_L_lo <- exp(LMADR_L-SE_LMADR_L*CI95)
MADR_L_up <- exp(LMADR_L+SE_LMADR_L*CI95)

# Mean absolute deviation ratio (right tail)
MADR_R_lo <- exp(LMADR_R-SE_LMADR_R*CI95)
MADR_R_up <- exp(LMADR_R+SE_LMADR_R*CI95)

# Gini's mean difference ratio
GMDR_lo <- exp(LGMDR-SE_LGMDR*CI95)
GMDR_up <- exp(LGMDR+SE_LGMDR*CI95)




##################
##### Output #####
##################

# further variables

WtRatio <- sum(T15_M$HWt)/sum(T15_F$HWt)

Low <- sum(T15[which(T15$Low == 1),'HWt'])/sum(T15$HWt)*100 # percent too low for estimation

# Standard deviation ratio (SDR), as VR is not comparable to MADR and GMDR
SDR <- sqrt(VR)
SDR_lo <- sqrt(VR_lo)
SDR_up <- sqrt(VR_up)
LSDR <- log(sqrt(VR))
TV_LSDR <- (log(SDR_up/SDR)/CI95)^2

SDR_L <- sqrt(VR_L)
SDR_L_lo <- sqrt(VR_L_lo)
SDR_L_up <- sqrt(VR_L_up)
LSDR_L <- log(sqrt(VR_L))
TV_LSDR_L <- (log(SDR_L_up/SDR_L)/CI95)^2

SDR_R <- sqrt(VR_R)
SDR_R_lo <- sqrt(VR_R_lo)
SDR_R_up <- sqrt(VR_R_up)
LSDR_R <- log(sqrt(VR_R))
TV_LSDR_R <- (log(SDR_R_up/SDR_R)/CI95)^2


# summary table
Labels <- c('CNT', 'Grade', 'Size', 'FSize', 'MSize', 'M/F Wt Ratio', 'Low %',
            'Mean', 'Mean Low', 'Mean Upp', 'Median', 'Median Low', 'Median Upp',
            'F Mean', 'F Mean Low', 'F Mean Upp', 'F Median', 'F Median Low', 'F Median Upp',
            'M Mean', 'M Mean Low', 'M Mean Upp', 'M Median', 'M Median Low', 'M Median Upp',
            'Mean Diff', 'Mean Diff Low', 'Mean Diff Upp',
            'Med Diff', 'Med Diff Low', 'Med Diff Upp',
            'TPRMn', 'TPRMn Low', 'TPRMn Upp', 'TPR05', 'TPR05 Low', 'TPR05 Upp',
            'TPR10', 'TPR10 Low', 'TPR10 Upp', 'TPR15', 'TPR15 Low', 'TPR15 Upp',
            'TPR20', 'TPR20 Low', 'TPR20 Upp', 'TPR25', 'TPR25 Low', 'TPR25 Upp',
            'TPR30', 'TPR30 Low', 'TPR30 Upp', 'TPR35', 'TPR35 Low', 'TPR35 Upp',
            'TPR40', 'TPR40 Low', 'TPR40 Upp', 'TPR45', 'TPR45 Low', 'TPR45 Upp',
            'TPR50', 'TPR50 Low', 'TPR50 Upp', 'TPR55', 'TPR55 Low', 'TPR55 Upp',
            'TPR60', 'TPR60 Low', 'TPR60 Upp', 'TPR65', 'TPR65 Low', 'TPR65 Upp',
            'TPR70', 'TPR70 Low', 'TPR70 Upp', 'TPR75', 'TPR75 Low', 'TPR75 Upp',
            'TPR80', 'TPR80 Low', 'TPR80 Upp', 'TPR85', 'TPR85 Low', 'TPR85 Upp',
            'TPR90', 'TPR90 Low', 'TPR90 Upp', 'TPR95', 'TPR95 Low', 'TPR95 Upp',
            'LTPRMn', 'LTPRMn TV', 'LTPR05', 'LTPR05 TV',
            'LTPR10', 'LTPR10 TV', 'LTPR15', 'LTPR15 TV',
            'LTPR20', 'LTPR20 TV', 'LTPR25', 'LTPR25 TV',
            'LTPR30', 'LTPR30 TV', 'LTPR35', 'LTPR35 TV',
            'LTPR40', 'LTPR40 TV', 'LTPR45', 'LTPR45 TV',
            'LTPR50', 'LTPR50 TV', 'LTPR55', 'LTPR55 TV',
            'LTPR60', 'LTPR60 TV', 'LTPR65', 'LTPR65 TV',
            'LTPR70', 'LTPR70 TV', 'LTPR75', 'LTPR75 TV',
            'LTPR80', 'LTPR80 TV', 'LTPR85', 'LTPR85 TV',
            'LTPR90', 'LTPR90 TV', 'LTPR95', 'LTPR95 TV',
            'Med95T', 'Med95T Low', 'Med95T Upp', 'Med95T TV',
            'Med90T', 'Med90T Low', 'Med90T Upp', 'Med90T TV',
            'Med10T', 'Med10T Low', 'Med10T Upp', 'Med10T TV',
            'Med05T', 'Med05T Low', 'Med05T Upp', 'Med05T TV',
            'Mn95T', 'Mn95T Low', 'Mn95T Upp', 'Mn95T TV',
            'Mn90T', 'Mn90T Low', 'Mn90T Upp', 'Mn90T TV',
            'Mn10T', 'Mn10T Low', 'Mn10T Upp', 'Mn10T TV',
            'Mn05T', 'Mn05T Low', 'Mn05T Upp', 'Mn05T TV',
            'U3R05', 'U3R05 Low', 'U3R05 Upp',
            'U3R10', 'U3R10 Low', 'U3R10 Upp', 'U3R15', 'U3R15 Low', 'U3R15 Upp',
            'U3R20', 'U3R20 Low', 'U3R20 Upp', 'U3R25', 'U3R25 Low', 'U3R25 Upp',
            'U3R30', 'U3R30 Low', 'U3R30 Upp', 'U3R35', 'U3R35 Low', 'U3R35 Upp',
            'U3R40', 'U3R40 Low', 'U3R40 Upp', 'U3R45', 'U3R45 Low', 'U3R45 Upp',
            'U3R50', 'U3R50 Low', 'U3R50 Upp', 'U3R55', 'U3R55 Low', 'U3R55 Upp',
            'U3R60', 'U3R60 Low', 'U3R60 Upp', 'U3R65', 'U3R65 Low', 'U3R65 Upp',
            'U3R70', 'U3R70 Low', 'U3R70 Upp', 'U3R75', 'U3R75 Low', 'U3R75 Upp',
            'U3R80', 'U3R80 Low', 'U3R80 Upp', 'U3R85', 'U3R85 Low', 'U3R85 Upp',
            'U3R90', 'U3R90 Low', 'U3R90 Upp', 'U3R95', 'U3R95 Low', 'U3R95 Upp',
            'LU3R05', 'LU3R05 TV',
            'LU3R10', 'LU3R10 TV', 'LU3R15', 'LU3R15 TV',
            'LU3R20', 'LU3R20 TV', 'LU3R25', 'LU3R25 TV',
            'LU3R30', 'LU3R30 TV', 'LU3R35', 'LU3R35 TV',
            'LU3R40', 'LU3R40 TV', 'LU3R45', 'LU3R45 TV',
            'LU3R50', 'LU3R50 TV', 'LU3R55', 'LU3R55 TV',
            'LU3R60', 'LU3R60 TV', 'LU3R65', 'LU3R65 TV',
            'LU3R70', 'LU3R70 TV', 'LU3R75', 'LU3R75 TV',
            'LU3R80', 'LU3R80 TV', 'LU3R85', 'LU3R85 TV',
            'LU3R90', 'LU3R90 TV', 'LU3R95', 'LU3R95 TV',
            'Med95U', 'Med95U Low', 'Med95U Upp', 'Med95U TV',
            'Med90U', 'Med90U Low', 'Med90U Upp', 'Med90U TV',
            'Med10U', 'Med10U Low', 'Med10U Upp', 'Med10U TV',
            'Med05U', 'Med05U Low', 'Med05U Upp', 'Med05U TV',
            'd', 'd Low', 'd Upp', 'd TV',
            'U3', 'U3 Low', 'U3 Upp', 'U3 TV',
            'PS', 'PS Low', 'PS Upp', 'PS TV',
            'VR', 'VR Low', 'VR Upp', 'LVR', 'LVR TV',
            'VR_L', 'VR_L Low', 'VR_L Upp', 'LVR_L', 'LVR_L TV',
            'VR_R', 'VR_R Low', 'VR_R Upp', 'LVR_R', 'LVR_R TV',
            'SDR', 'SDR Low', 'SDR Upp', 'LSDR', 'LSDR TV',
            'SDR_L', 'SDR_L Low', 'SDR_L Upp', 'LSDR_L', 'LSDR_L TV',
            'SDR_R', 'SDR_R Low', 'SDR_R Upp', 'LSDR_R', 'LSDR_R TV',
            'MADR', 'MADR Low', 'MADR Upp', 'LMADR', 'LMADR TV',
            'MADR_L', 'MADR_L Low', 'MADR_L Upp', 'LMADR_L', 'LMADR_L TV',
            'MADR_R', 'MADR_R Low', 'MADR_R Upp', 'LMADR_R', 'LMADR_R TV',
            'GMDR', 'GMDR Low', 'GMDR Upp', 'LGMDR', 'LGMDR TV',
            'Age U3', 'Age MADR', 'Age LMADR', 'd AgeCor', 'U3 AgeCor',
            'PS AgeCor', 'VR AgeCor', 'LVR AgeCor', 'VR_L AgeCor',
            'LVR_L AgeCor', 'VR_R AgeCor', 'LVR_R AgeCor',
            'MADR AgeCor', 'LMADR AgeCor', 'MADR_L AgeCor', 'LMADR_L AgeCor',
            'MADR_R AgeCor', 'LMADR_R AgeCor', 'GMDR AgeCor', 'LGMDR AgeCor')

Variables <- c(CNT, Grade, Size, FSize, MSize, WtRatio, Low,
               Mn, Mn_lo, Mn_up, Md, Md_lo, Md_up,
               Mn_F, Mn_F_lo, Mn_F_up, Md_F, Md_F_lo, Md_F_up,
               Mn_M, Mn_M_lo, Mn_M_up, Md_M, Md_M_lo, Md_M_up,
               MnDf, MnDf_lo, MnDf_up,
               MdDf, MdDf_lo, MdDf_up,
               TPRMn, TPRMn_lo, TPRMn_up, TPR05, TPR05_lo, TPR05_up,
               TPR10, TPR10_lo, TPR10_up, TPR15, TPR15_lo, TPR15_up,
               TPR20, TPR20_lo, TPR20_up, TPR25, TPR25_lo, TPR25_up,
               TPR30, TPR30_lo, TPR30_up, TPR35, TPR35_lo, TPR35_up,
               TPR40, TPR40_lo, TPR40_up, TPR45, TPR45_lo, TPR45_up,
               TPR50, TPR50_lo, TPR50_up, TPR55, TPR55_lo, TPR55_up,
               TPR60, TPR60_lo, TPR60_up, TPR65, TPR65_lo, TPR65_up,
               TPR70, TPR70_lo, TPR70_up, TPR75, TPR75_lo, TPR75_up,
               TPR80, TPR80_lo, TPR80_up, TPR85, TPR85_lo, TPR85_up,
               TPR90, TPR90_lo, TPR90_up, TPR95, TPR95_lo, TPR95_up,
               LTPRMn, TV_Mn, LTPR05, TV_05T,
               LTPR10, TV_10T, LTPR15, TV_15T,
               LTPR20, TV_20T, LTPR25, TV_25T,
               LTPR30, TV_30T, LTPR35, TV_35T,
               LTPR40, TV_40T, LTPR45, TV_45T,
               LTPR50, TV_50T, LTPR55, TV_55T,
               LTPR60, TV_60T, LTPR65, TV_65T,
               LTPR70, TV_70T, LTPR75, TV_75T,
               LTPR80, TV_80T, LTPR85, TV_85T,
               LTPR90, TV_90T, LTPR95, TV_95T,
               Md95T, Md95T_lo, Md95T_up, TV_Md95T,
               Md90T, Md90T_lo, Md90T_up, TV_Md90T,
               Md10T, Md10T_lo, Md10T_up, TV_Md10T,
               Md05T, Md05T_lo, Md05T_up, TV_Md05T,
               Mn95T, Mn95T_lo, Mn95T_up, TV_Mn95T,
               Mn90T, Mn90T_lo, Mn90T_up, TV_Mn90T,
               Mn10T, Mn10T_lo, Mn10T_up, TV_Mn10T,
               Mn05T, Mn05T_lo, Mn05T_up, TV_Mn05T,
               U3R05, U3R05_lo, U3R05_up,
               U3R10, U3R10_lo, U3R10_up, U3R15, U3R15_lo, U3R15_up,
               U3R20, U3R20_lo, U3R20_up, U3R25, U3R25_lo, U3R25_up,
               U3R30, U3R30_lo, U3R30_up, U3R35, U3R35_lo, U3R35_up,
               U3R40, U3R40_lo, U3R40_up, U3R45, U3R45_lo, U3R45_up,
               U3R50, U3R50_lo, U3R50_up, U3R55, U3R55_lo, U3R55_up,
               U3R60, U3R60_lo, U3R60_up, U3R65, U3R65_lo, U3R65_up,
               U3R70, U3R70_lo, U3R70_up, U3R75, U3R75_lo, U3R75_up,
               U3R80, U3R80_lo, U3R80_up, U3R85, U3R85_lo, U3R85_up,
               U3R90, U3R90_lo, U3R90_up, U3R95, U3R95_lo, U3R95_up,
               LU3R05, TV_05U,
               LU3R10, TV_10U, LU3R15, TV_15U,
               LU3R20, TV_20U, LU3R25, TV_25U,
               LU3R30, TV_30U, LU3R35, TV_35U,
               LU3R40, TV_40U, LU3R45, TV_45U,
               LU3R50, TV_50U, LU3R55, TV_55U,
               LU3R60, TV_60U, LU3R65, TV_65U,
               LU3R70, TV_70U, LU3R75, TV_75U,
               LU3R80, TV_80U, LU3R85, TV_85U,
               LU3R90, TV_90U, LU3R95, TV_95U,
               Md95U, Md95U_lo, Md95U_up, TV_Md95U,
               Md90U, Md90U_lo, Md90U_up, TV_Md90U,
               Md10U, Md10U_lo, Md10U_up, TV_Md10U,
               Md05U, Md05U_lo, Md05U_up, TV_Md05U,
               d, d_lo, d_up, TV_d,
               U3, U3_lo, U3_up, TV_U3,
               PS, PS_lo, PS_up, TV_PS,
               VR, VR_lo, VR_up, LVR, TV_LVR,
               VR_L, VR_L_lo, VR_L_up, LVR_L, TV_LVR_L,
               VR_R, VR_R_lo, VR_R_up, LVR_R, TV_LVR_R,
               SDR, SDR_lo, SDR_up, LSDR, TV_LSDR,
               SDR_L, SDR_L_lo, SDR_L_up, LSDR_L, TV_LSDR_L,
               SDR_R, SDR_R_lo, SDR_R_up, LSDR_R, TV_LSDR_R,
               MADR, MADR_lo, MADR_up, LMADR, TV_LMADR,
               MADR_L, MADR_L_lo, MADR_L_up, LMADR_L, TV_LMADR_L,
               MADR_R, MADR_R_lo, MADR_R_up, LMADR_R, TV_LMADR_R,
               GMDR, GMDR_lo, GMDR_up, LGMDR, TV_LGMDR, 
               AgeU3, AgeMADR, AgeLMADR, d_A, U3_A,
               PS_A, VR_A, LVR_A, VR_LA,
               LVR_LA, VR_RA, LVR_RA,
               MADR_A, LMADR_A, MADR_LA, LMADR_LA,
               MADR_RA, LMADR_RA, GMDR_A, LGMDR_A)

Output <- format(data.frame(Labels, Variables), scientific = F) # put everything in this

write.csv(x = Output, file = 'TIMSS output/2015 4/AUS.csv') # select file name, store as a csv



















