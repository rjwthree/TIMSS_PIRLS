######################################################################
######################################################################
############### SEX DIFFERENCES IN THE PIRLS 2016 DATA ###############
######################################################################
######################################################################


#### Outline ####

# [1] Read and Format Data
# select country
# objects needed throughout the script
# [2] Weighted Functions
# functions needed throughout the script, with descriptions
# [3] Effect Sizes
# various effect sizes for sex differences throughout the distribution
# [4] Some Effect Sizes Adjusted for Age
# dataframe with scores linearly controlled for age
# age-adjusted effect sizes
# [5] Standard Errors and Confidence Intervals
# jackknife resampling
# total variance and CI bounds for effect sizes
# [6] Output
# write dataframe with all needed variables

# ratios are log-transformed to place them on a linear scale
# log-transformed ratios are signified by an 'L' preceding the acronym




################################
##### Read and Format Data #####
################################

#### Select Country ####

# PIRLS Literacy (L) is an easier version of the standard PIRLS assessment (R)
# ePIRLS (E) is a computer assessment; all polities that did ePIRLS also did PIRLS
# Norway administered the standard assessment to 5th graders,
# and tested a 4th grade benchmark sample of similar size


# manually subset a test and a polity (benchmarks at the end)

test <- 'R' # PIRLS — 47/50 countries, 10/11 benchmarks
test <- 'E' # ePIRLS — 14/50 countries, 2/11 benchmarks
test <- 'L' # PIRLS Literacy — 5/50 countries, 1/11 benchmarks

CNT <- 'AUS' # Australia | R
CNT <- 'AUT' # Austria | R
CNT <- 'AZE' # Azerbaijan | R
CNT <- 'BHR' # Bahrain | R
CNT <- 'BFL' # Belgium (Flemish) | R
CNT <- 'BFR' # Belgium (French) | R
CNT <- 'BGR' # Bulgaria | R
CNT <- 'CAN' # Canada | R E
CNT <- 'CHL' # Chile | R
CNT <- 'CZE' # Czech Republic | R
CNT <- 'DNK' # Denmark | R E
CNT <- 'EGY' # Egypt | L
CNT <- 'ENG' # England | R
CNT <- 'FIN' # Finland | R
CNT <- 'FRA' # France | R
CNT <- 'GEO' # Georgia | R E
CNT <- 'DEU' # Germany | R
CNT <- 'HKG' # Hong Kong | R
CNT <- 'HUN' # Hungary | R
CNT <- 'IRN' # Iran | R L
CNT <- 'IRL' # Ireland | R E
CNT <- 'ISR' # Israel | R E
CNT <- 'ITA' # Italy | R E
CNT <- 'KAZ' # Kazakhstan | R
CNT <- 'KWT' # Kuwait | L
CNT <- 'LVA' # Latvia | R
CNT <- 'LTU' # Lithuania | R
CNT <- 'MAC' # Macao, China | R
CNT <- 'MLT' # Malta | R
CNT <- 'MAR' # Morocco | R L
CNT <- 'NLD' # Netherlands | R
CNT <- 'NZL' # New Zealand | R
CNT <- 'NIR' # Northern Ireland | R
CNT <- 'NOR' # Norway | R E
CNT <- 'OMN' # Oman | R
CNT <- 'POL' # Poland | R
CNT <- 'PRT' # Portugal | R E
CNT <- 'QAT' # Qatar | R
CNT <- 'RUS' # Russia | R
CNT <- 'SAU' # Saudi Arabia | R
CNT <- 'SGP' # Singapore | R E
CNT <- 'SVK' # Slovakia | R
CNT <- 'SVN' # Slovenia | R E
CNT <- 'ZAF' # South Africa | L
CNT <- 'ESP' # Spain | R
CNT <- 'SWE' # Sweden | R E
CNT <- 'TWN' # Taiwan | R E
CNT <- 'TTO' # Trinidad and Tobago | R
CNT <- 'ARE' # United Arab Emirates | R E
CNT <- 'USA' # United States | R E

CNT <- 'ABA' # Buenos Aires, Argentina | R
CNT <- 'COT' # Ontario, Canada | R
CNT <- 'CQU' # Quebec, Canada | R
CNT <- 'DN3' # Denmark (3rd grade) | L
CNT <- 'NO4' # Norway (4th grade) | R
CNT <- 'RMO' # Moscow city, Russia | R
CNT <- 'ZA5' # South Africa (5th grade) | R
CNT <- 'EAN' # Andalusia, Spain | R
CNT <- 'EMA' # Madrid, Spain | R
CNT <- 'AAD' # Abu Dhabi, United Arab Emirates | R E
CNT <- 'ADU' # Dubai, United Arab Emirates | R E



#### Create Objects ####

library(haven) # read SPSS

# lookup table
CT <- data.frame('Australia', 'Austria', 'Azerbaijan', 'Bahrain', 'Belgium (Fl)',
                 'Belgium (Fr)', 'Bulgaria', 'Canada', 'Chile', 'Czech Republic',
                 'Denmark', 'Egypt', 'England', 'Finland', 'France',
                 'Georgia', 'Germany', 'Hong Kong', 'Hungary', 'Iran',
                 'Ireland', 'Israel', 'Italy', 'Kazakhstan', 'Kuwait',
                 'Latvia', 'Lithuania', 'Macao', 'Malta', 'Morocco',
                 'Netherlands', 'New Zealand', 'Northern Ireland', 'Norway', 'Oman',
                 'Poland', 'Portugal', 'Qatar', 'Russia', 'Saudi Arabia',
                 'Singapore', 'Slovakia', 'Slovenia', 'South Africa', 'Spain',
                 'Sweden', 'Taiwan', 'Trinidad and Tobago', 'United Arab Emirates', 'United States',
                 'Buenos Aires', 'Ontario', 'Quebec', 'Denmark (3)', 'Norway (4)',
                 'Moscow', 'South Africa (5)', 'Andalusia', 'Madrid', 'Abu Dhabi', 'Dubai')
names(CT) <- c('AUS', 'AUT', 'AZE', 'BHR', 'BFL', 'BFR', 'BGR', 'CAN', 'CHL', 'CZE',
               'DNK', 'EGY', 'ENG', 'FIN', 'FRA', 'GEO', 'DEU', 'HKG', 'HUN', 'IRN',
               'IRL', 'ISR', 'ITA', 'KAZ', 'KWT', 'LVA', 'LTU', 'MAC', 'MLT', 'MAR',
               'NLD', 'NZL', 'NIR', 'NOR', 'OMN', 'POL', 'PRT', 'QAT', 'RUS', 'SAU',
               'SGP', 'SVK', 'SVN', 'ZAF', 'ESP', 'SWE', 'TWN', 'TTO', 'ARE', 'USA',
               'ABA', 'COT', 'CQU', 'DN3', 'NO4', 'RMO', 'ZA5', 'EAN', 'EMA', 'AAD', 'ADU')

Country <- as.character(CT[,CNT])

if (test == 'E') { # ePIRLS
  PIRLS16 <- read_spss(paste0('PIRLS/2016/4E/ASG', CNT, 'E1.sav'))
  P16 <- PIRLS16[c('IDSTUD', 'ITSEX', 'ASDAGE', 'HOUWGT', 'JKZONE', 'JKREP', 'ASDELOWP',
                   'ASEREA01', 'ASEREA02', 'ASEREA03', 'ASEREA04', 'ASEREA05')] # subset columns
} else { # PIRLS and PIRLS Literacy have the same column names
  if (test == 'R') { # PIRLS
    PIRLS16 <- read_spss(paste0('PIRLS/2016/4/ASG', CNT, 'R4.sav'))
  } else if (test == 'L') {PIRLS16 <- read_spss(paste0('PIRLS/2016/4L/ASG', CNT, 'L2.sav'))} # Literacy
  P16 <- PIRLS16[c('IDSTUD', 'ITSEX', 'ASDAGE', 'HOUWGT', 'JKZONE', 'JKREP', 'ASDRLOWP',
                   'ASRREA01', 'ASRREA02', 'ASRREA03', 'ASRREA04', 'ASRREA05')] # subset columns
}

names(P16) <- c('Student', 'Sex', 'Age', 'HWt', 'JKZ', 'JKR', 'Low',
                'PV1', 'PV2', 'PV3', 'PV4', 'PV5') # rename columns

for (i in names(P16)) {attributes(P16[[i]])$label <- NULL} # delete column labels
P16 <- data.frame(zap_labels(P16)) # remove all labels and convert to dataframe

Size <- nrow(P16) # sample size
P16$HWt <- Size/sum(P16$HWt)*P16$HWt # equalize sample size and sum of HWt

P16_F <- P16[which(P16$Sex == 1),] # female subset
P16_M <- P16[which(P16$Sex == 2),] # male subset

FSize <- nrow(P16_F) # female sample size
MSize <- nrow(P16_M) # male sample size




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
  i <- wt <- 0
  if (q > 0 & q <= .5) {
    ascend <- order(x)
    a <- x[ascend]
    b <- w[ascend]
    if (k <= b[1]) {return(min(x))}
    while (wt < k) {
      i <- i + 1
      wt <- wt + b[i]
    }
    if (i == length(a)) {
      if (a[i] == a[i-1]) {i <- -1} else {return(a[i]-(wt-k)/b[i]*(a[i]-a[i-1]))}
    } else if (a[i] == a[i-1] | a[i] == a[i+1]) {
      i <- -1
    } else {return(a[i]-(wt-k)/b[i]*(a[i]-a[i-1]))}
  } else if (q > .5 & q < 1) { 
    descend <- order(x, decreasing = T)
    a <- x[descend]
    b <- w[descend]
    if (k <= b[length(b)]) {return(min(x))}
    wt <- sum(w)
    while (wt > k) {
      i <- i + 1
      wt <- wt - b[i]
    }
    if (i+1 == length(a)) {
      if (a[i+1] == a[i]) {i <- -1} else {return(a[i+1]+(k-wt)/b[i]*(a[i]-a[i+1]))}
    } else if (a[i+1] == a[i] | a[i+1] == a[i+2]) {
      i <- -1
    } else {return(a[i+1]+(k-wt)/b[i]*(a[i]-a[i+1]))}
  } else if (q == 0) {return(min(x))}
  
  if (i < 0) {
    d <- data.frame(x, w)
    W <- data.frame(table(d$x))
    X <- W[W$Freq > 1,]
    Y <- as.numeric(as.character(X$Var1))
    Z <- numeric(length(Y))
    for (i in 1:length(Y)) {Z[i] <- sum(d[d$x == Y[i],'w'])}
    d <- d[! d$x %in% Y,]
    d <- rbind(d, data.frame(x = Y, w = Z))
    
    i <- wt <- 0
    ascend <- order(d$x)
    a <- d$x[ascend]
    b <- d$w[ascend]
    if (k <= b[1]) {return(min(x))}
    while (wt < k) {
      i <- i + 1
      wt <- wt + b[i]
    }
    return(a[i]-(wt-k)/b[i]*(a[i]-a[i-1]))
  }
}
# sort the scores and weights in ascending order of scores (a and b)
# add weights until the fraction of total weight exceeds q
# perform linear interpolation between the closest scores
# for efficiency, reverse the process if q exceeds .5
# if needed, divert to an alternate algorithm by assigning -1 to i (impossible otherwise)
# this sums the weights of students with identical scores before proceeding
# it is a less efficient algorithm that remedies the slight distortion which occurs in
# rare cases when a[i] (if q <= .5) or a[i+1] (if q > .5) is one of multiple equal scores
# other weighted quantile functions using linear interpolation would produce a different
# result if, e.g., one simply duplicated all values and weights, which should not happen
# edge cases are tested for and addressed; one of the tests is necessary only
# because of the adjustments made to check whether the alternate algorithm is needed


# empirical cumulative distribution function
wt.ecdf <- function(x, w, p) {
  ascend <- order(x)
  a <- x[ascend]
  b <- w[ascend]
  n <- length(x)
  z <- sum(a <= p)
  
  if (z == 0) {
    return(0)
  } else if (z+1 == n) {
    s <- (p-a[z]) / (a[z+1]-a[z]) * b[z+1]
    return((s + sum(b[1:z])) / sum(b))
  } else if (z == n) {return(1)}
  
  if (a[z+1] == a[z+2]) {
    d <- data.frame(a, b)
    wt <- sum(d[d$a == a[z+1],'b'])
  } else {wt <- b[z+1]}
  
  s <- (p-a[z]) / (a[z+1]-a[z]) * wt
  return((s + sum(b[1:z])) / sum(b))
}
# sort the scores and weights in ascending order of scores (a and b)
# find the number of scores less than or equal to p (z); their weight is sum(b[1:z])
# perform linear interpolation between the closest scores to find the additional weight (s)
# the fraction of total weight gives the weighted quantile of p
# the rare case in which two or more students with identical scores have
# the score immediately above p is addressed by summing their weights


# Cohen's d
dfn <- function(d1, d2, v) {
  (wt.mn(d1[,v], d1$HWt) - wt.mn(d2[,v], d2$HWt))/sqrt((wt.var(d1[,v], d1$HWt) + wt.var(d2[,v], d2$HWt))/2)
}
# raw mean difference divided by quadratic mean of standard deviations


# U3
U3fn <- function(d1, d2, v) {1 - wt.ecdf(d1[,v], d1$HWt, wt.qnt(d2[,v], d2$HWt, .5))}
# the precise share of male weight above the female median


# probability of superiority (PS)
PSfn <- function(d1, d2, v) {
  d1 <- d1[order(d1[,v]),]
  SuperM <- numeric(FSize)
  for (i in 1:FSize) {SuperM[i] <- sum(d1[which(d1[,v] > d2[i,v]),'HWt'])}
  
  equal <- intersect(d1[,v], d2[,v])
  EWtsM <- EWtsF <- numeric(length(equal))
  for (i in 1:length(equal)) {
    EWtsM[i] <- sum(d1[which(d1[,v] %in% equal[i]),'HWt'])
    EWtsF[i] <- sum(d2[which(d2[,v] %in% equal[i]),'HWt'])
  }
  
  return(sum(c(SuperM*d2$HWt, .5*EWtsM*EWtsF)) / (sum(d1$HWt)*sum(d2$HWt)))
}
# for efficiency, sort male scores and weights in order of scores (d1)
# compute the male weight above each female score (SuperM)
# store the small subset of scores that are equal to a score in the other subgroup (equal)
# find the associated weights in males and females (EWtsM, EWtsF)
# sum weights of superior males multiplied by associated female weights
# and half of the sum of equal male and female weights multiplied
# then divide by the total weight of male-female pairs
# this gives the probability that a random male has a higher score than a random female


# log-transformed standard deviation ratio (LSDR)
LSDRfn <- function(d1, d2, v) {log(sqrt(wt.var(d1[,v], d1$HWt) / wt.var(d2[,v], d2$HWt)))}
# M/F ratio of standard deviation
# log-transform the ratio for linear scale


# log-transformed tail SDR (LSDR_T)
LSDR_Tfn <- function(d1, d2, v, t) {
  q1 <- wt.mn(d1[,v], d1$HWt)
  q2 <- wt.mn(d2[,v], d2$HWt)
  if (t == 'L') {
    Mtail <- d1[which(d1[,v] < q1),]
    Ftail <- d2[which(d2[,v] < q2),]
  } else if (t == 'R') {
    Mtail <- d1[which(d1[,v] > q1),]
    Ftail <- d2[which(d2[,v] > q2),]
  }
  m <- sum(Mtail$HWt*(Mtail[,v]-q1)^2) / sum(Mtail$HWt)
  f <- sum(Ftail$HWt*(Ftail[,v]-q2)^2) / sum(Ftail$HWt)
  return(log(sqrt(m/f)))
}
# subset the males and females below (t = 'L') or above (t = 'R') the subgroup mean
# compute sqrt of the M/F ratio of mean squared deviation from the mean in the tail
# log-transform the ratio for linear scale
# Bessel's correction not applicable because in this case the sample mean's
# deviation from the population mean causes random not systematic error


# log-transformed median absolute deviation ratio (LMADR)
LMADRfn <- function(d1, d2, v) {
  log(wt.qnt(abs(d1[,v]-wt.qnt(d1[,v], d1$HWt, .5)), d1$HWt, .5)/
        wt.qnt(abs(d2[,v]-wt.qnt(d2[,v], d2$HWt, .5)), d2$HWt, .5))
}
# M/F ratio of median absolute deviation from the median
# log-transform the ratio for linear scale


# log-transformed tail MADR (LMADR_T)
LMADR_Tfn <- function(d1, d2, v, t) {
  q1 <- wt.qnt(d1[,v], d1$HWt, .5)
  q2 <- wt.qnt(d2[,v], d2$HWt, .5)
  if (t == 'L') {
    Mtail <- d1[which(d1[,v] < q1),]
    Ftail <- d2[which(d2[,v] < q2),]
    return(log(wt.qnt(q1-Mtail[,v], Mtail$HWt, .5)/wt.qnt(q2-Ftail[,v], Ftail$HWt, .5)))
  } else if (t == 'R') {
    Mtail <- d1[which(d1[,v] > q1),]
    Ftail <- d2[which(d2[,v] > q2),]
    return(log(wt.qnt(Mtail[,v]-q1, Mtail$HWt, .5)/wt.qnt(Ftail[,v]-q2, Ftail$HWt, .5)))
  }
}
# subset the males and females below (t = 'L') or above (t = 'R') the subgroup median
# compute the M/F ratio of mean absolute deviation from the median in the left or right tail
# log-transform the ratio for linear scale


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
  return(log(Gini1*wt.mn(d1[,v], d1$HWt)/(Gini2*wt.mn(d2[,v], d2$HWt))))
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
# log-transform the ratio for linear scale


# log-transformed U3 ratio (LU3R)
LU3Rfn <- function(d1, d2, q, v) {
  m <- wt.ecdf(d1[,v], d1$HWt, wt.qnt(d2[,v], d2$HWt, q))
  if (q < .5) {return(log(m/q))} else {return(log((1-m)/(1-q)))}
}
# share of male weight below a female subgroup quantile (m)
# if left tail, find the ratio of weight share below the female quantile; if right tail, above it
# log-transform the ratio for linear scale


# log-transformed tail proportion ratio (LTPR)
LTPRfn <- function(d1, d2, q, v) {
  d <- rbind(d1, d2)
  k <- wt.qnt(d[,v], d$HWt, q)
  m <- wt.ecdf(d1[,v], d1$HWt, k)
  f <- wt.ecdf(d2[,v], d2$HWt, k)
  if (q < .5) {return(log(m/f))} else {return(log((1-m)/(1-f)))}
}
# compute proportions of male and female weight below a threshold k (m and f)
# if left tail, find the TPR below k; if right tail, above k 
# log-transform the ratio for linear scale


# log-transformed median-aligned U3 ratio (LMU3R)
LMU3Rfn <- function(d1, d2, q, v) {
  d1[,v] <- wt.qnt(d2[,v], d2$HWt, .5) - wt.qnt(d1[,v], d1$HWt, .5) + d1[,v]
  m <- wt.ecdf(d1[,v], d1$HWt, wt.qnt(d2[,v], d2$HWt, q))
  if (q < .5) {return(log(m/q))} else {return(log((1-m)/(1-q)))}
}
# move the male median to the female median (arbitrary method of alignment)
# compute the LMU3R (see LU3R description)


# quantile difference (QD)
QDfn <- function(d1, d2, q, v) {wt.qnt(d1[,v], d1$HWt, q) - wt.qnt(d2[,v], d2$HWt, q)}
# raw quantile M-F difference


# mean MAD (MMAD)
MMADfn <- function(d1, d2, v) {
  (wt.qnt(abs(d1[,v]-wt.qnt(d1[,v], d1$HWt, .5)), d1$HWt, .5) +
     wt.qnt(abs(d2[,v]-wt.qnt(d2[,v], d2$HWt, .5)), d2$HWt, .5)) / 2
}
# mean of male and female MADs, computed once outside of loop for efficiency


# standardized quantile difference (SQD)
SQDfn <- function(d1, d2, q, v) {
  QD <- wt.qnt(d1[,v], d1$HWt, q) - wt.qnt(d2[,v], d2$HWt, q)
  m <- wt.qnt(abs(d1[,v]-wt.qnt(d1[,v], d1$HWt, .5)), d1$HWt, .5)
  f <- wt.qnt(abs(d2[,v]-wt.qnt(d2[,v], d2$HWt, .5)), d2$HWt, .5)
  return(QD / ((m+f)/2) * 100)
}
# raw quantile M-F difference as a percentage of mean MAD: QDfn / MMADfn * 100




########################
##### Effect Sizes #####
########################

### means and medians (M), some effect sizes (E)

Ms <- numeric(30) # empty container
M <- numeric(8) # empty container
Es <- numeric(50) # empty container
E <- numeric(10) # empty container

for (i in 1:5) {
  PV <- paste0('PV',i)
  Ms[i] <- wt.mn(P16[,PV], P16$HWt) # means (total)
  Ms[i+5] <- wt.qnt(P16[,PV], P16$HWt, .5) # medians (total)
  Ms[i+10] <- wt.mn(P16_F[,PV], P16_F$HWt) # means (female)
  Ms[i+15] <- wt.qnt(P16_F[,PV], P16_F$HWt, .5) # medians (female)
  Ms[i+20] <- wt.mn(P16_M[,PV], P16_M$HWt) # means (male)
  Ms[i+25] <- wt.qnt(P16_M[,PV], P16_M$HWt, .5) # medians (male)
  Es[i] <- dfn(P16_M, P16_F, PV) # Cohen's ds
  Es[i+5] <- U3fn(P16_M, P16_F, PV) # U3s
  Es[i+10] <- PSfn(P16_M, P16_F, PV) # PSs
  Es[i+15] <- LSDRfn(P16_M, P16_F, PV) # LSDRs
  Es[i+20] <- LSDR_Tfn(P16_M, P16_F, PV, 'L') # LSDR_Ls
  Es[i+25] <- LSDR_Tfn(P16_M, P16_F, PV, 'R') # LSDR_Rs
  Es[i+30] <- LMADRfn(P16_M, P16_F, PV) # LMADRs
  Es[i+35] <- LMADR_Tfn(P16_M, P16_F, PV, 'L') # LMADR_Ls
  Es[i+40] <- LMADR_Tfn(P16_M, P16_F, PV, 'R') # LMADR_Rs
  Es[i+45] <- LGMDRfn(P16_M, P16_F, PV) # LGMDRs
}

Ms <- c(Ms, Ms[21:25]-Ms[11:15], Ms[26:30]-Ms[16:20]) # mean and median differences

for (i in 1:8) {M[i] <- mean(Ms[(5*i-4):(5*i)])} # means and medians
for (i in 1:10) {E[i] <- mean(Es[(5*i-4):(5*i)])} # effect sizes



### U3 Ratios (U3Rs) and Tail Proportion Ratios (TPRs)

Rs <- numeric(40) # empty container
R <- numeric(8) # empty container
P <- c(.05, .1, .9, .95) # quantiles

for (i in 1:4) { # for percentiles 5, 10, 90, 95
  for (s in 1:5) { # for each PV
    Rs[s+(i-1)*5] <- LU3Rfn(P16_M, P16_F, P[i], paste0('PV',s)) # LU3Rs
    Rs[s+(i+3)*5] <- LTPRfn(P16_M, P16_F, P[i], paste0('PV',s)) # LTPRs
  }
  R[i] <- mean(Rs[(5*i-4):(5*i)]) # LU3R for each percentile
  R[i+4] <- mean(Rs[(5*i+16):(5*i+20)]) # LTPR for each percentile
}



### Median-aligned U3 Ratios (MU3Rs) and Standardized Quantile Differences (SQDs)

LMU3Rs <- SQDs <- numeric(495) # empty containers
LMU3R <- SQD <- numeric(99) # empty containers
MMADs <- numeric(5) # empty container

for (i in 1:5) {MMADs[i] <- MMADfn(P16_M, P16_F, paste0('PV',i))} # mean MADs

for (i in 1:99) { # for each percentile
  for (s in 1:5) { # for each PV
    LMU3Rs[s+(i-1)*5] <- LMU3Rfn(P16_M, P16_F, i/100, paste0('PV',s)) # LMU3Rs
    SQDs[s+(i-1)*5] <- QDfn(P16_M, P16_F, i/100, paste0('PV',s)) / MMADs[s] * 100 # SQDs
  }
  LMU3R[i] <- mean(LMU3Rs[(5*i-4):(5*i)]) # LMU3R for each percentile
  SQD[i] <- mean(SQDs[(5*i-4):(5*i)]) # SQD for each percentile
}

# SQD tail-center shifts (SQDTCs): Med-5, Med-10, 90-Med, 95-Med
Qs <- c(SQDs[246:250]-SQDs[21:25], SQDs[246:250]-SQDs[46:50],
        SQDs[446:450]-SQDs[246:250], SQDs[471:475]-SQDs[246:250])
Q <- c(SQD[50]-SQD[5], SQD[50]-SQD[10], SQD[90]-SQD[50], SQD[95]-SQD[50])

# ratios for standard errors: all LU3Rs and LTPRs, LMU3Rs at 5, 10, 90, 95
Rs <- c(Rs, LMU3Rs[c(21:25, 46:50, 446:450, 471:475)])
R <- c(R, LMU3R[c(5, 10, 90, 95)])




##############################################
##### Some Effect Sizes Adjusted for Age #####
##############################################

A16 <- P16[!is.na(P16$Age),] # P16 with Age NAs (if any) excluded
A16_F <- A16[which(A16$Sex == 1),] # female subset
A16_M <- A16[which(A16$Sex == 2),] # male subset

AgeU3 <- U3fn(A16_M, A16_F, 'Age') # sex differences in age
AgeMADR <- exp(LMADRfn(A16_M, A16_F, 'Age'))

Slope <- numeric(5) # control for the correlation between age and score
S <- c(PV1 ~ Age, PV2 ~ Age, PV3 ~ Age, PV4 ~ Age, PV5 ~ Age)
for (i in 1:5) {Slope[i] <- lm(formula = S[[i]], data = A16, weights = HWt)$coefficients[2]}

AgeMn <- wt.mn(A16$Age, A16$HWt) # mean age

# if there are Age NAs, replace with mean age so their nominally age-corrected scores do not change
if (length(P16[is.na(P16$Age),'Age']) != 0) {P16[is.na(P16$Age),'Age'] <- AgeMn}

A16 <- data.frame(P16, # P16 with age-corrected scores
                  PV1A = P16$PV1+(AgeMn-P16$Age)*Slope[1],
                  PV2A = P16$PV2+(AgeMn-P16$Age)*Slope[2],
                  PV3A = P16$PV3+(AgeMn-P16$Age)*Slope[3],
                  PV4A = P16$PV4+(AgeMn-P16$Age)*Slope[4],
                  PV5A = P16$PV5+(AgeMn-P16$Age)*Slope[5])
A16_F <- A16[which(A16$Sex == 1),] # female subset
A16_M <- A16[which(A16$Sex == 2),] # male subset

EAs <- numeric(50) # empty container
EA <- numeric(10) # empty container

for (i in 1:5) {
  PV <- paste0('PV',i,'A')
  EAs[i] <- dfn(A16_M, A16_F, PV) # Cohen's ds
  EAs[i+5] <- U3fn(A16_M, A16_F, PV) # U3s
  EAs[i+10] <- PSfn(A16_M, A16_F, PV) # PSs
  EAs[i+15] <- LSDRfn(A16_M, A16_F, PV) # LSDRs
  EAs[i+20] <- LSDR_Tfn(A16_M, A16_F, PV, 'L') # LSDR_Ls
  EAs[i+25] <- LSDR_Tfn(A16_M, A16_F, PV, 'R') # LSDR_Rs
  EAs[i+30] <- LMADRfn(A16_M, A16_F, PV) # LMADRs
  EAs[i+35] <- LMADR_Tfn(A16_M, A16_F, PV, 'L') # LMADR_Ls
  EAs[i+40] <- LMADR_Tfn(A16_M, A16_F, PV, 'R') # LMADR_Rs
  EAs[i+45] <- LGMDRfn(A16_M, A16_F, PV) # LGMDRs
}

for (i in 1:10) {EA[i] <- mean(EAs[(5*i-4):(5*i)])} # age-adjusted effect sizes




####################################################
##### Standard Errors and Confidence Intervals #####
####################################################

L <- length(unique(P16$JKZ)) # number of jackknife zones
L1 <- rep(1:L, times = 2) # vector to select jackknife zones
L2 <- rep(c(1,0), each = L) # vector to select jackknife replicate codes
L3 <- rev(L2) # vector to select jackknife replicate codes
L <- 2*L # double L for subsequent use

Y <- c(.05, .1, .5, .9, .95) # quantiles for SQDTCs

# empty containers
MJ <- numeric(30*L)
EJ <- numeric(50*L)
QJ <- numeric(25*L)
RJ <- numeric(60*L)
Msum <- numeric(40)
Esum <- numeric(50)
Qsum <- numeric(20)
Rsum <- numeric(60)
TVM <- numeric(8)
TVE <- numeric(10)
TVQ <- numeric(4)
TVR <- numeric(12)
MCI <- numeric(16)
ECI <- numeric(20)
QCI <- numeric(8)
RCI <- numeric(24)

# perform jackknife resampling
for (i in 1:L) { # for each JK zone, twice
  P0 <- P16 # create/restore duplicate
  P0[which(P0$JKZ == L1[i] & P0$JKR == L2[i]),'HWt'] <- 2*P0[which(P0$JKZ == L1[i] & P0$JKR == L2[i]),'HWt']
  P0[which(P0$JKZ == L1[i] & P0$JKR == L3[i]),'HWt'] <- 0 # reweight according to JK replicate code
  P0_F <- P0[which(P0$Sex == 1),] # female subset
  P0_M <- P0[which(P0$Sex == 2),] # male subset
  
  for (s in 1:5) { # for each PV
    PV <- paste0('PV',s)
    
    MJ[i+(s-1)*L] <- wt.mn(P0[,PV], P0$HWt) # reweighted means (total)
    MJ[i+(s+4)*L] <- wt.qnt(P0[,PV], P0$HWt, .5) # reweighted medians (total)
    MJ[i+(s+9)*L] <- wt.mn(P0_F[,PV], P0_F$HWt) # reweighted means (female)
    MJ[i+(s+14)*L] <- wt.qnt(P0_F[,PV], P0_F$HWt, .5) # reweighted medians (female)
    MJ[i+(s+19)*L] <- wt.mn(P0_M[,PV], P0_M$HWt) # reweighted means (male)
    MJ[i+(s+24)*L] <- wt.qnt(P0_M[,PV], P0_M$HWt, .5) # reweighted medians (male)
    
    EJ[i+(s-1)*L] <- dfn(P0_M, P0_F, PV) # reweighted Cohen's ds
    EJ[i+(s+4)*L] <- U3fn(P0_M, P0_F, PV) # reweighted U3s
    EJ[i+(s+9)*L] <- PSfn(P0_M, P0_F, PV) # reweighted PSs
    EJ[i+(s+14)*L] <- LSDRfn(P0_M, P0_F, PV) # reweighted LSDRs
    EJ[i+(s+19)*L] <- LSDR_Tfn(P0_M, P0_F, PV, 'L') # reweighted LSDR_Ls
    EJ[i+(s+24)*L] <- LSDR_Tfn(P0_M, P0_F, PV, 'R') # reweighted LSDR_Rs
    EJ[i+(s+29)*L] <- LMADRfn(P0_M, P0_F, PV) # reweighted LMADRs
    EJ[i+(s+34)*L] <- LMADR_Tfn(P0_M, P0_F, PV, 'L') # reweighted LMADR_Ls
    EJ[i+(s+39)*L] <- LMADR_Tfn(P0_M, P0_F, PV, 'R') # reweighted LMADR_Rs
    EJ[i+(s+44)*L] <- LGMDRfn(P0_M, P0_F, PV) # reweighted LGMDRs
    
    for (c in 1:5) {QJ[i+(s+5*c-6)*L] <- SQDfn(P0_M, P0_F, Y[c], PV)} # reweighted SQDs
    
    for (c in 1:4) {
      RJ[i+(s+5*c-6)*L] <- LU3Rfn(P0_M, P0_F, P[c], PV) # reweighted LU3Rs
      RJ[i+(s+5*c+14)*L] <- LTPRfn(P0_M, P0_F, P[c], PV) # reweighted LTPRs
      RJ[i+(s+5*c+34)*L] <- LMU3Rfn(P0_M, P0_F, P[c], PV) # reweighted LMU3Rs
    }
  }
  print(paste0(i, '/', L, ' at ', Sys.time()), quote = F) # print updates
}

# append reweighted mean/median differences to MJ, replace reweighted SQDs in QJ with reweighted SQDTCs
MJ <- c(MJ, MJ[(20*L+1):(25*L)]-MJ[(10*L+1):(15*L)], MJ[(25*L+1):(30*L)]-MJ[(15*L+1):(20*L)])
QJ <- c(QJ[(10*L+1):(15*L)]-QJ[1:(5*L)], QJ[(10*L+1):(15*L)]-QJ[(5*L+1):(10*L)],
        QJ[(15*L+1):(20*L)]-QJ[(10*L+1):(15*L)], QJ[(20*L+1):(25*L)]-QJ[(10*L+1):(15*L)])

CI95 <- qnorm(.975) # ratio of 95% CIs to SEs

# total variance = sampling variance + imputation variance
for (i in 1:12) { # for each effect size, up to the appropriate i
  for (s in 1:5) { # for each PV
    if (i < 9) {Msum[s+(i-1)*5] <- sum((MJ[((s+5*i-6)*L+1):((s+(i-1)*5)*L)]-Ms[s+(i-1)*5])^2)}
    if (i < 11) {Esum[s+(i-1)*5] <- sum((EJ[((s+5*i-6)*L+1):((s+(i-1)*5)*L)]-Es[s+(i-1)*5])^2)}
    if (i < 5) {Qsum[s+(i-1)*5] <- sum((QJ[((s+5*i-6)*L+1):((s+(i-1)*5)*L)]-Qs[s+(i-1)*5])^2)}
    Rsum[s+(i-1)*5] <- sum((RJ[((s+5*i-6)*L+1):((s+(i-1)*5)*L)]-Rs[s+(i-1)*5])^2)
  }
  if (i < 9) { # total variance and 95% CI bounds of means and medians
    TVM[i] <- mean(Msum[(5*i-4):(5*i)])/2 + .3*sum((Ms[(5*i-4):(5*i)]-M[i])^2)
    MCI[2*i-1] <- M[i]-sqrt(TVM[i])*CI95
    MCI[2*i] <- M[i]+sqrt(TVM[i])*CI95
  }
  if (i < 11) { # total variance and 95% CI bounds of some effect sizes
    TVE[i] <- mean(Esum[(5*i-4):(5*i)])/2 + .3*sum((Es[(5*i-4):(5*i)]-E[i])^2)
    ECI[2*i-1] <- E[i]-sqrt(TVE[i])*CI95
    ECI[2*i] <- E[i]+sqrt(TVE[i])*CI95
  }
  if (i < 5) { # total variance and 95% CI bounds of SQDTCs
    TVQ[i] <- mean(Qsum[(5*i-4):(5*i)])/2 + .3*sum((Qs[(5*i-4):(5*i)]-Q[i])^2)
    QCI[2*i-1] <- Q[i]-sqrt(TVQ[i])*CI95
    QCI[2*i] <- Q[i]+sqrt(TVQ[i])*CI95
  }            # total variance and 95% CI bounds of U3Rs, TPRs, MU3Rs
  TVR[i] <- mean(Rsum[(5*i-4):(5*i)])/2 + .3*sum((Rs[(5*i-4):(5*i)]-R[i])^2)
  RCI[2*i-1] <- R[i]-sqrt(TVR[i])*CI95
  RCI[2*i] <- R[i]+sqrt(TVR[i])*CI95
}




##################
##### Output #####
##################

WtRatio <- sum(P16_M$HWt)/sum(P16_F$HWt) # M/F weight ratio

MMAD <- mean(MMADs) # mean MAD

Low <- sum(P16[which(P16$Low == 1),'HWt'])/sum(P16$HWt)*100 # percent too low for estimation

LMU3Rnm <- SQDnm <- numeric(99)
for (i in 1:99) { # compressed names for LMU3Rs and SQDs
  LMU3Rnm[i] <- paste0('LMU3R',i)
  SQDnm[i] <- paste0('SQD',i)
}

# summary table
Names <- c('Country', 'CNT', 'Size', 'FSize', 'MSize', 'M/F Wt Ratio', 'Low %',
           'Mean', 'Median', 'F Mean', 'F Median', 'M Mean', 'M Median', 'Mean Diff', 'Median Diff', 
           'd', 'U3', 'PS', 'SDR', 'SDR_L', 'SDR_R', 'MADR', 'MADR_L', 'MADR_R', 'GMDR',
           'U3R05', 'U3R10', 'U3R90', 'U3R95', 'TPR05', 'TPR10', 'TPR90', 'TPR95', 'MU3R05', 'MU3R10',
           'MU3R90', 'MU3R95', 'MMAD', 'SQDTC05', 'SQDTC10', 'SQDTC90', 'SQDTC95', LMU3Rnm, SQDnm,
           'Mean Low', 'Mean Upp', 'Median Low', 'Median Upp', 'F Mean Low', 'F Mean Upp', 'F Median Low',
           'F Median Upp', 'M Mean Low', 'M Mean Upp', 'M Median Low', 'M Median Upp', 'Mean Diff Low',
           'Mean Diff Upp', 'Med Diff Low', 'Med Diff Upp', 'd Low', 'd Upp', 'U3 Low', 'U3 Upp',
           'PS Low', 'PS Upp', 'SDR Low', 'SDR Upp', 'SDR_L Low', 'SDR_L Upp', 'SDR_R Low', 'SDR_R Upp',
           'MADR Low', 'MADR Upp', 'MADR_L Low', 'MADR_L Upp', 'MADR_R Low', 'MADR_R Upp', 'GMDR Low',
           'GMDR Upp', 'SQDTC05 Low', 'SQDTC05 Upp', 'SQDTC10 Low', 'SQDTC10 Upp', 'SQDTC90 Low',
           'SQDTC90 Upp', 'SQDTC95 Low', 'SQDTC95 Upp', 'U3R05 Low', 'U3R05 Upp', 'U3R10 Low', 'U3R10 Upp',
           'U3R90 Low', 'U3R90 Upp', 'U3R95 Low', 'U3R95 Upp', 'TPR05 Low', 'TPR05 Upp', 'TPR10 Low',
           'TPR10 Upp', 'TPR90 Low', 'TPR90 Upp', 'TPR95 Low', 'TPR95 Upp', 'MU3R05 Low', 'MU3R05 Upp',
           'MU3R10 Low', 'MU3R10 Upp', 'MU3R90 Low', 'MU3R90 Upp', 'MU3R95 Low', 'MU3R95 Upp', 'LSDR',
           'LSDR_L', 'LSDR_R', 'LMADR', 'LMADR_L', 'LMADR_R', 'LGMDR', 'LU3R05', 'LU3R10', 'LU3R90',
           'LU3R95', 'LTPR05', 'LTPR10', 'LTPR90', 'LTPR95', 'd TV', 'U3 TV', 'PS TV', 'LSDR TV',
           'LSDR_L TV', 'LSDR_R TV', 'LMADR TV', 'LMADR_L TV', 'LMADR_R TV', 'LGMDR TV', 'SQDTC05 TV',
           'SQDTC10 TV', 'SQDTC90 TV', 'SQDTC95 TV', 'LU3R05 TV', 'LU3R10 TV', 'LU3R90 TV', 'LU3R95 TV',
           'LTPR05 TV', 'LTPR10 TV', 'LTPR90 TV', 'LTPR95 TV', 'LMU3R05 TV', 'LMU3R10 TV', 'LMU3R90 TV',
           'LMU3R95 TV', 'Age U3', 'Age MADR', 'd A', 'U3 A', 'PS A', 'SDR A', 'SDR_L A', 'SDR_R A',
           'MADR A', 'MADR_L A', 'MADR_R A', 'GMDR A')

Variables <- c(Country, CNT, Size, FSize, MSize, WtRatio, Low, M, E[1:3], exp(c(E[4:10], R)), MMAD,
               Q, LMU3R, SQD, MCI, ECI[1:6], exp(ECI[7:20]), QCI, exp(RCI), E[4:10], R[1:8], TVE, TVQ, TVR,
               AgeU3, AgeMADR, EA[1:3], exp(EA[4:10]))

Output <- format(data.frame(Names, Variables), scientific = F) # put everything in this

# select file name manually, store as a csv
write.csv(x = Output, file = paste0('PIRLS output/2016/', test, '/Countries/', CNT, '.csv'))




















