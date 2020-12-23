Work on this repository paused until [TIMSS 2019 data](https://timss2019.org/reports/) are released, at an unspecified point in the near future.

# Sex Differences in TIMSS 2015 Math and PIRLS 2016 Reading
These scripts analyze the most recently released data on 4th and 8th grade math exams from the [Trends in International Mathematics and Science Study (TIMSS)](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/TIMSS%202015.R), and 4th grade reading exams from the [Progress in International Reading Literacy Study (PIRLS)](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/PIRLS%202016.R). TIMSS is administered every four years and PIRLS every five years by the International Association for the Evaluation of Educational Achievement (IEA). https://timssandpirls.bc.edu/

The first major section describes the results; the second section describes the code and statistics. See [here](https://github.com/rjwthree/PISA_2018) for similar assessments of 15-year-olds by the OECD.

## Results

The full output from all 50 countries from [PIRLS](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/Full%20results%20(reading).csv) and 57 countries from [TIMSS](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/Full%20results%20(math).csv) is available in csv files, but is difficult to read because there are so many analyses and countries. I've visualized some of the results and hyperlinked them in the discussion below. See further down for definitions of all my effect sizes.

#### Gender inequality

Although flawed, the World Economic Forum's Global Gender Gap Index (GGGI) is arguably the best quantification of women's disadvantages around the world. The reports published in [Oct 2016](https://www.weforum.org/reports/the-global-gender-gap-report-2016) and [Nov 2017](https://www.weforum.org/reports/the-global-gender-gap-report-2017), which correspond best with the assessment timelines of TIMSS 2015 and PIRLS 2016, respectively, [contain](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/GGGI%20(math).csv) [data](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/GGGI%20(reading).csv) on nearly all the countries tested by [TIMSS](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/Country%20codes%20(math).csv) and [PIRLS](https://github.com/rjwthree/TIMSS_PIRLS/blob/master/Country%20codes%20(reading).csv). (The exception is Taiwan: China's insistence that other countries officially deny its sovereignty continues to be effective.)

The index is an average of four subindices which can range from 0 to 1. For each subindex, a score of 1 would indicate equality or female advantage. Spearman's correlation coefficient (_r<sub>s</sub>_) will be used here because non-parametric statistics are most suitable.

## Code

The scripts proceed as follows:

(1) Read and Format Data

(2) Weighted User-defined Functions

(3) Effect Sizes and effect sizes adjusted for age

(4) Standard Errors and Confidence Intervals

(5) Summary and Output

Effect sizes were selected (or created) to be robust, intuitive and useful. They examine sex differences throughout the distribution, not just at the center. Common effect sizes which are not robust or intuitive, like Cohen's _d_, are also used for the sake of external comparison. An attitude of [multiverse analysis](https://journals.sagepub.com/doi/full/10.1177/1745691616658637) is applied: when multiple reasonable analyses exist, they are all conducted so that the sensitivity of substantive conclusions to specific analytic choices can be observed.

The data are complex for two reasons: (1) each student responds to only a subset of all questions, so each student is given five imputed scores, also known as [plausible values](https://timssandpirls.bc.edu/publications/timss/2015-methods/chapter-12.html) (PVs); and (2) the sampling procedures employ clustering and stratification, so a [jackknife method](https://timssandpirls.bc.edu/publications/timss/2015-methods/chapter-4.html) is used to compute standard errors. Due to the first point, standard errors incorporate not only sampling variance but also imputation variance. The sampling design also means that each student is assigned their own weight. The data files are designed for SPSS and SAS and are analyzed almost exclusively as such.

Ratios are log-transformed to place them on a linear scale.

## Read and Format Data
The data are read from SPSS format using the 'haven' package and converted to a dataframe. Some objects are created that will be needed throughout the script.

## Weighted Functions
All functions used for analysis are user-defined, for various reasons related to transparency, efficiency, and necessity - some effect sizes had to be created because I wasn't aware of any that were suitable (see next section).

## Effect Sizes
*Effect sizes with an asterisk are novel; they were adapted from other statistics.

Means and Medians — Means and medians of total group, females, and males, as well as male-female mean and median differences.

Cohen's _d_ — This is a very common effect size and it is implemented as usual: the raw mean difference divided by the quadratic mean of standard deviations.

U3* — The proportion of males above the female median. Adapted from the canonical Cohen's U3, which is incorrectly based on the mean rather than the median.

Probability of Superiority (PS) — The probability that a random male will have a higher score than a random female.

Median Absolute Deviation Ratio (MADR) — The male-female ratio of median absolute deviation (MAD) from the median.

Left MADR and Right MADR* — The MADR divided into the left and right tails. That is, the male-female ratio of the median absolute deviation from the median among the subset of scores below or above the median, respectively.

Standard Deviation Ratio (SDR) — The ratio of male SD to female SD.

Left SDR and Right SDR* — The SDR divided into the left and right tails. That is, the male-female ratio of the square root of the mean squared deviation from the mean among the subset of scores below or above the mean, respectively.

Gini's Mean Difference Ratio (GMDR) — Imagine randomly selecting two students and finding the distance between their scores. GMD gives the expected value of this interval: the mean absolute difference. The GMDR is the male-female ratio.

Tail Proportion Ratios (TPRs) — The proportion of males below (left tail) or above (right tail) a threshold, divided by the proportion of females below or above the threshold. The thresholds are defined by the total group. The left tail TPRs are below the 5th and 10th percentiles, and the right tail TPRs are above the 90th and 95th percentiles.

U3 Ratios (U3Rs)* — The proportion of males below (left tail) or above (right tail) a threshold defined by the female subgroup, divided by the natural proportion of females below or above that threshold (e.g., 10% above the 90th percentile). The left tail U3Rs are below the 5th and 10th percentiles, and the right tail U3Rs are above the 90th and 95th percentiles. Adapted from U3 (see above).

Median-aligned U3 Ratios (MU3Rs)* — U3Rs calculated after aligning the male and female medians, for each percentile. The MU3R below/above the median is 1 by definition. A U shape (values decrease to the left of the median and increase to the right) indicates higher male variability, and an upside-down U shape (values increase to the left of the median and decrease to the right) indicates higher female variability.

Standardized Quantile Differences (SQDs)* — Raw male-female differences at each percentile, as a percentage of the mean of male and female MADs. An increasing trend from left to right indicates higher male variability; a decreasing trend indicates higher female variability. The purpose of dividing by MAD is to standardize the quantile differences with a robust measure of scale.

The next section adjusts the scores linearly for age and recalculates the effect sizes from Cohen's _d_ to GMDR.

## Standard Errors and Confidence Intervals
Computing standard errors is the most computationally intensive procedure. The use of 75 jackknife zones, each subdivided into two groups — in combination with the five sets of plausible values — means that each effect size must be recalculated 750 times. The sampling variance and imputation variance are then calculated according to procedures described by '[Methods and Procedures in TIMSS 2015](https://timssandpirls.bc.edu/publications/timss/2015-methods.html)' and '[Methods and Procedures in PIRLS 2016](https://timssandpirls.bc.edu/publications/pirls/2016-methods.html)'. The sum of sampling and imputation variance is the total variance, and its square root is the standard error. The standard errors are then converted to 95% confidence intervals through typical methods, and the lower and upper bounds of these confidence intervals are recorded.

The results are then summarized in a dataframe called 'Output' and written to the appropriate directory.
