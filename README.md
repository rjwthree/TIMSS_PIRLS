# Sex Differences in TIMSS 2015 Math
These are the most recently released data on 4th and 8th grade math exams from the quadrennial Trends in International Mathematics and Science Study (TIMSS).

The script proceeds as follows: (1) reading and formatting data; (2) weighted user-defined functions; (3) a variety of effect sizes expressing sex differences in central tendency, in variability, and in the tails; (4) effect sizes adjusted for age; (5) standard error and confidence interval calculations; and (6) writing a compilation of results as a csv file.

The data are complex for two reasons: (1) each student responds to only a subset of all questions, so each student is given five imputed scores, also known as plausible values (PVs); and (2) the sampling procedures employ clustering and stratification, so a jackknife method is needed to compute standard errors. Due to the first point standard errors incorporate not only sampling variance but also imputation variance. The sampling design also means that each student is assigned their own weight. The data files are designed for SPSS and SAS and are analyzed almost exclusively as such.

Ratios are log-transformed to place them on a linear scale.

## Read and Format Data
The data are read from SPSS format using the 'haven' package and converted to a dataframe. Some objects are created that will be needed throughout the script.

## Weighted Functions
All functions used for analysis are user-defined, for a few reasons:

(1) Some functions were not available, especially for weighted statistics and effect sizes I created because I wasn't aware of any that were suitable (see next section).

(2) Some existing functions were inefficient and runtime was far too long because of the large sample sizes and extensive repetition during jackknife resampling. All efforts were made to maximize the functions' efficiency.

(3) User-defined functions are entirely transparent and explicit, which is particularly good when different packages produce slightly different results due to differing methods.

## Effect Sizes
*Effect sizes with an asterisk do not exist elsewhere; they were adapted from other statistics.

Means and Medians - Means and medians of total group, females, and males, as well as male-female mean and median differences.

Tail Proportion Ratios (TPRs) - The proportion of males above a threshold divided by the proportion of females above that threshold. The thresholds are defined by the total group: the mean and every 5th percentile from 5 to 95.

LTPR Tail-Center Differences* - Differences between log-transformed TPRs (LTPRs) are used to express changes in sex differences throughout the distribution. In the right tail, the differences are between LTPRs at the 95th percentile and median, 90th percentile and median, 95th percentile and mean, 90th percentile and mean. In the left tail, the differences are between LTPRs at the median and 10th percentile, median and 5th percentile, mean and 10th percentile, mean and 5th percentile.

U3 Ratios (U3Rs)* - The proportion of males above a threshold defined by the female subgroup, divided by the natural proportion of females above that threshold (e.g., 75% above the 25th percentile). The thresholds are every 5th percentile from 5 to 95. Adapted from Cohen's U3.

LU3R Tail-Center Differences* - Differences between log-transformed U3Rs (LU3Rs) are used to express changes in sex differences throughout the distribution. In the right tail, the differences are between LU3Rs at the 95th percentile and median and at the 90th percentile and median. In the left tail, the differences are between LU3Rs at the median and 10th percentile and at the median and 5th percentile.

Cohen's d - This is a very common effect size and it is implemented as usual: the raw mean difference divided by the quadratic mean of standard deviations.

U3* - The proportion of males above the female median. Adapted from Cohen's U3, which is based on the mean rather than median.

Probability of Superiority (PS) - The probability that a random male will have a higher score than a random female.

Variance Ratio (VR) - The male-female variance ratio.

Left VR and Right VR* - The VR divided into the left and right tails. That is, the mean squared deviation from the mean among the subset of scores below or above the mean, respectively.

Mean Absolute Deviation Ratio (MADR) - The male-female ratio of mean absolute deviation from the median.

Left MADR and Right MADR* - The MADR divided into the left and right tails. That is, the mean absolute deviation from the median among the subset of scores below or above the median, respectively.

Gini's Mean Difference Ratio (GMDR) - Imagine randomly selecting two students and finding the absolute difference between their scores; GMD gives the expected value of this interval. The GMDR is the male-female ratio of this distance.
