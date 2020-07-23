# Sex Differences in TIMSS 2015 Math
These are the most recently released data on 4th and 8th grade math exams from the Trends in International Mathematics and Science Study (TIMSS). It is administered every four years by the International Association for the Evaluation of Educational Achievement (IEA), which also administers the Progress in International Reading Literacy Study (PIRLS) every five years. https://timssandpirls.bc.edu/

The script proceeds as follows:

(1) Read and Format Data

(2) Weighted User-defined Functions

(3) Effect Sizes and effect sizes adjusted for age

(4) Standard Errors and Confidence Intervals

(5) Summary and Output

Effect sizes were selected (or created) to be robust, intuitive and useful. They examine sex differences throughout the distribution, not just at the center. Common effect sizes (e.g., Cohen's d) are also used. An attitude of 'multiverse analysis' is applied: when multiple reasonable analyses exist, they are all conducted so that the sensitivity of substantive conclusions to specific analytic choices can be observed.

The data are complex for two reasons: (1) each student responds to only a subset of all questions, so each student is given five imputed scores, also known as plausible values (PVs); and (2) the sampling procedures employ clustering and stratification, so a jackknife method is used to compute standard errors. Due to the first point, standard errors incorporate not only sampling variance but also imputation variance. The sampling design also means that each student is assigned their own weight. The data files are designed for SPSS and SAS and are analyzed almost exclusively as such.

Ratios are log-transformed to place them on a linear scale.

## Read and Format Data
The data are read from SPSS format using the 'haven' package and converted to a dataframe. Some objects are created that will be needed throughout the script.

## Weighted Functions
All functions used for analysis are user-defined, for a few reasons:

(1) Some functions were not available, especially for weighted statistics and effect sizes I created because I wasn't aware of any that were suitable (see next section).

(2) Some existing functions were inefficient and runtime was far too long because of the large sample sizes and extensive repetition during jackknife resampling. All efforts were made to maximize the functions' efficiency.

(3) User-defined functions are entirely transparent and explicit, which is particularly good when different packages produce slightly different results due to differing methods.

## Effect Sizes
*Effect sizes with an asterisk are novel; they were adapted from other statistics.

Means and Medians - Means and medians of total group, females, and males, as well as male-female mean and median differences.

Cohen's d - This is a very common effect size and it is implemented as usual: the raw mean difference divided by the quadratic mean of standard deviations.

U3* - The proportion of males above the female median. Adapted from the canonical Cohen's U3, which is incorrectly based on the mean rather than the median.

Probability of Superiority (PS) - The probability that a random male will have a higher score than a random female.

Standard Deviation Ratio (SDR) - The ratio of male SD to female SD.

Left SDR and Right SDR* - The SDR divided into the left and right tails. That is, the square root of the mean squared deviation from the mean among the subset of scores below or above the mean, respectively.

Median Absolute Deviation Ratio (MADR) - The male-female ratio of median absolute deviation (MAD) from the median.

Left MADR and Right MADR* - The MADR divided into the left and right tails. That is, the median absolute deviation from the median among the subset of scores below or above the median, respectively.

Gini's Mean Difference Ratio (GMDR) - Imagine randomly selecting two students and finding the distance between their scores. GMD gives the expected value of this interval: the mean absolute difference. The GMDR is the male-female ratio.

Tail Proportion Ratios (TPRs) - The proportion of males below (left tail) or above (right tail) a threshold, divided by the proportion of females below or above the threshold. The thresholds are defined by the total group. The left tail TPRs are below the 5th and 10th percentiles, and the right tail TPRs are above the 90th and 95th percentiles.

U3 Ratios (U3Rs)* - The proportion of males below (left tail) or above (right tail) a threshold defined by the female subgroup, divided by the natural proportion of females below or above that threshold (e.g., 10% above the 90th percentile). The left tail U3Rs are below the 5th and 10th percentiles, and the right tail U3Rs are above the 90th and 95th percentiles. Adapted from U3 (see above).

Median-aligned U3 Ratios (MU3Rs)* - U3Rs calculated after aligning the male and female medians, for each percentile. The MU3R below/above the median is 1 by definition. A U shape (values decrease to the left of the median and increase to the right) indicates higher male variability, and an upside-down U shape (values increase to the left of the median and decrease to the right) indicates higher female variability.

Standardized Quantile Differences (SQDs)* - Raw male-female differences at each percentile, as a percentage of the mean of male and female MADs. An increasing trend from left to right indicates higher male variability; a decreasing trend indicates higher female variability. The purpose of dividing by MAD is to standardize the quantile differences with a robust measure of scale.

The next section adjusts the scores linearly for age and recalculates the effect sizes from Cohen's d to GMDR.

## Standard Errors and Confidence Intervals
Computing standard errors is the most computationally intensive procedure. The use of 75 jackknife zones, each subdivided into two groups — in combination with the five sets of plausible values — means that each effect size must be recalculated 750 times. The sampling variance and imputation variance are then calculated according to formulas described by 'Methods and Procedures in TIMSS 2015'. The sum of sampling and imputation variance is the total variance, and its square root is the standard error. The standard errors are then converted to 95% confidence intervals through typical methods, and the lower and upper bounds of these confidence intervals are recorded.

The results are then summarized in a dataframe called 'Output' and written to the appropriate directory.
