# Sex Differences in TIMSS 2015 Math
These are the most recently released data on math exams from the quadrennial Trends in International Mathematics and Science Study (TIMSS) of 4th and 8th grade students.

The script proceeds as follows:

(1) Read and Format Data

(2) Weighted User-defined Functions

(3) A variety of effect sizes expressing sex differences in central tendency, in variability, and in the tails. Scores are then controlled for age and effect sizes are recalculated.

(4) Standard Errors and Confidence Intervals

(5) Summary and output

The data are complex for two reasons: (1) each student responds to only a subset of all questions, so each student is given five imputed scores, also known as plausible values (PVs); and (2) the sampling procedures employ clustering and stratification, so a jackknife method is needed to compute standard errors. Due to the first point standard errors incorporate not only sampling variance but also imputation variance. The sampling design also means that each student is assigned their own weight. The data files are designed for SPSS and SAS and are analyzed almost exclusively as such.

Effect sizes were selected (or created) to be robust, intuitive and useful. They examine sex differences throughout the distribution, not just at the center. Common effect sizes (e.g., Cohen's d) were also used. An attitude of 'multiverse analysis' was used: when multiple reasonable analyses exist, they are all conducted so that the sensitivity of substantive conclusions to specific analytic choices can be observed.

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

The following section adjusts the scores linearly for age and recalculates all the effect sizes from Cohen's d down.

## Standard Errors and Confidence Intervals
This section is the most computationally intensive. The jackknife resampling section, in combination with the five plausible values, means that each effect size must be reweighted and recalculated 750 times. The jackknife variance and imputation variance are then calculated according to procedures described by 'Methods and Procedures in TIMSS 2015'. The sum of jackknife and imputation variance is the total variance, and its square root is the standard error. The standard errors are then converted to 95% confidence intervals through typical methods, and the lower and upper bounds of these confidence intervals are computed.

The results are then summarized in a dataframe called 'Output' and written to the working directory.
