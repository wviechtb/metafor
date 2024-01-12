\name{conv.fivenum}
\alias{conv.fivenum}
\title{Estimate Means and Standard Deviations from Five-Number Summary Values}
\description{
   Function to estimate means and standard deviations from five-number summary values.
}
\usage{
conv.fivenum(min, q1, median, q3, max, n, data, include,
             method="default", dist="norm", transf=TRUE, test=TRUE,
             var.names=c("mean","sd"), append=TRUE, replace="ifna", \dots)
}
\arguments{
   \item{min}{vector with the minimum values.}
   \item{q1}{vector with the lower/first quartile values.}
   \item{median}{vector with the median values.}
   \item{q3}{vector with the upper/third quartile values.}
   \item{max}{vector with the maximum values.}
   \item{n}{vector with the sample sizes.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{include}{optional (logical or numeric) vector to specify the subset of studies for which means and standard deviations should be estimated.}
   \item{method}{character string indicating the method to use. Either \code{"default"} (same as \code{"luo/wan/shi"} which is the current default), \code{"qe"}, \code{"bc"}, \code{"mln"}, or \code{"blue"}. Can be abbreviated. See \sQuote{Details}.}
   \item{dist}{character string indicating the distribution assumed for the underlying data (either \code{"norm"} for a normal distribution or \code{"lnorm"} for a log-normal distribution). Can also be a string vector if different distributions are assumed for different studies. Only relevant when \code{method="default"}.}
   \item{transf}{logical to specify whether the estimated means and standard deviations of the log-transformed data should be back-transformed as described by Shi et al. (2020b) (the default is \code{TRUE}). Only relevant when \code{dist="lnorm"} and when \code{method="default"}.}
   \item{test}{logical to specify whether a study should be excluded from the estimation if the test for skewness is significant (the default is \code{TRUE}, but whether this is applicable depends on the method; see \sQuote{Details}).}
   \item{var.names}{character vector with two elements to specify the name of the variable for the estimated means and the name of the variable for the estimated standard deviations (the defaults are \code{"mean"} and \code{"sd"}).}
   \item{append}{logical to specify whether the data frame provided via the \code{data} argument should be returned together with the estimated values (the default is \code{TRUE}).}
   \item{replace}{character string or logical to specify how values in \code{var.names} should be replaced (only relevant when using the \code{data} argument and if variables in \code{var.names} already exist in the data frame). See the \sQuote{Value} section for more details.}
   \item{\dots}{other arguments.}
}
\details{
   Various effect size measures require means and standard deviations (SDs) as input (e.g., raw or standardized mean differences, ratios of means / response ratios; see \code{\link{escalc}} for further details). For some studies, authors may not report means and SDs, but other statistics, such as the so-called \sQuote{five-number summary}, consisting of the minimum, lower/first quartile, median, upper/third quartile, and the maximum of the sample values (plus the sample sizes). Occasionally, only a subset of these values are reported.

   The present function can be used to estimate means and standard deviations from five-number summary values based on various methods described in the literature (Bland, 2015; Cai et al. 2021; Hozo et al., 2005; Luo et al., 2016; McGrath et al., 2020; Shi et al., 2020a; Walter & Yao, 2007; Wan et al., 2014; Yang et al., 2022).

   When \code{method="default"} (which is the same as \code{"luo/wan/shi"}), the following methods are used:

   \subsection{Case 1: Min, Median, Max}{

      In case only the minimum, median, and maximum is available for a study (plus the sample size), then the function uses the method by Luo et al. (2016), equation (7), to estimate the mean and the method by Wan et al. (2014), equation (9), to estimate the SD.

   }

   \subsection{Case 2: Q1, Median, Q3}{

      In case only the lower/first quartile, median, and upper/third quartile is available for a study (plus the sample size), then the function uses the method by Luo et al. (2016), equation (11), to estimate the mean and the method by Wan et al. (2014), equation (16), to estimate the SD.

   }

   \subsection{Case 3: Min, Q1, Median, Q3, Max}{

      In case the full five-number summary is available for a study (plus the sample size), then the function uses the method by Luo et al. (2016), equation (15), to estimate the mean and the method by Shi et al. (2020a), equation (10), to estimate the SD.

   }

   ---------

   The median is not actually needed in the methods by Wan et al. (2014) and Shi et al. (2020a) and hence it is possible to estimate the SD even if the median is unavailable (this can be useful if a study reports the mean directly, but instead of the SD, it reports the minimum/maximum and/or first/third quartile values).

   Note that the sample size must be at least 5 to apply these methods. Studies where the sample size is smaller are not included in the estimation. The function also checks that \code{min <= q1 <= median <= q3 <= max} and throws an error if any studies are found where this is not the case.

   \subsection{Test for Skewness}{

      The methods described above were derived under the assumption that the data are normally distributed. Testing this assumption would require access to the raw data, but based on the three cases above, Shi et al. (2023) derived tests for skewness that only require the reported quantile values and the sample sizes. These tests are automatically carried out. When \code{test=TRUE} (which is the default), a study is automatically excluded from the estimation if the test is significant. If all studies should be included, set \code{test=FALSE}, but note that the accuracy of the methods will tend to be poorer when the data come from an apparently skewed (and hence non-normal) distribution.

   }

   \subsection{Log-Normal Distribution}{

      When setting \code{dist="lnorm"}, the raw data are assumed to follow a log-normal distribution. In this case, the methods as described by Shi et al. (2020b) are used to estimate the mean and SD of the log transformed data for the three cases above. When \code{transf=TRUE} (the default), the estimated mean and SD of the log transformed data are back-transformed to the estimated mean and SD of the raw data (using the bias-corrected back-transformation as described by Shi et al., 2020b). Note that the test for skewness is also carried out when \code{dist="lnorm"}, but now testing if the log transformed data exhibit skewness.

   }

   \subsection{Alternative Methods}{

      As an alternative to the methods above, one can make use of the methods implemented in the \href{https://cran.r-project.org/package=estmeansd}{estmeansd} package to estimate means and SDs based on the three cases above. Available are the quantile estimation method (\code{method="qe"}; using the \code{\link[estmeansd]{qe.mean.sd}} function; McGrath et al., 2020), the Box-Cox method (\code{method="bc"}; using the \code{\link[estmeansd]{bc.mean.sd}} function; McGrath et al., 2020), and the method for unknown non-normal distributions (\code{method="mln"}; using the \code{\link[estmeansd]{mln.mean.sd}} function; Cai et al. 2021). The advantage of these methods is that they do not assume that the data underlying the reported values are normally distributed (and hence the \code{test} argument is ignored), but they can only be used when the values are positive (except for the quantile estimation method, which can also be used when one or more of the values are negative, but in this case the method does assume that the data are normally distributed and hence the test for skewness is applied when \code{test=TRUE}). Note that all of these methods may struggle to provide sensible estimates when some of the values are equal to each other (which can happen when the data include a lot of ties and/or the reported values are rounded). Also, the Box-Cox method and the method for unknown non-normal distributions involve simulated data and hence results will slightly change on repeated runs. Setting the seed of the random number generator (with \code{\link{set.seed}}) ensures reproducibility.

      Finally, by setting \code{method="blue"}, one can make use of the \code{\link[metaBLUE]{BLUE_s}} function from the \href{https://cran.r-project.org/package=metaBLUE}{metaBLUE} package to estimate means and SDs based on the three cases above (Yang et al., 2022). The method assumes that the underlying data are normally distributed (and hence the test for skewness is applied when \code{test=TRUE}).

   }
}
\value{
   If the \code{data} argument was not specified or \code{append=FALSE}, a data frame with two variables called \code{var.names[1]} (by default \code{"mean"}) and \code{var.names[2]} (by default \code{"sd"}) with the estimated means and SDs.

   If \code{data} was specified and \code{append=TRUE}, then the original data frame is returned. If \code{var.names[1]} is a variable in \code{data} and \code{replace="ifna"} (or \code{replace=FALSE}), then only missing values in this variable are replaced with the estimated means (where possible) and otherwise a new variable called \code{var.names[1]} is added to the data frame. Similarly, if \code{var.names[2]} is a variable in \code{data} and \code{replace="ifna"} (or \code{replace=FALSE}), then only missing values in this variable are replaced with the estimated SDs (where possible) and otherwise a new variable called \code{var.names[2]} is added to the data frame.

   If \code{replace="all"} (or \code{replace=TRUE}), then all values in \code{var.names[1]} and \code{var.names[2]} where an estimated mean and SD can be computed are replaced, even for cases where the value in \code{var.names[1]} and \code{var.names[2]} is not missing.

   When missing values in \code{var.names[1]} are replaced, an attribute called \code{"est"} is added to the variable, which is a logical vector that is \code{TRUE} for values that were estimated. The same is done when missing values in \code{var.names[2]} are replaced.

   Attributes called \code{"tval"}, \code{"crit"}, \code{"sig"}, and \code{"dist"} are also added to \code{var.names[1]} corresponding to the test statistic and critical value for the test for skewness, whether the test was significant, and the assumed distribution (for the quantile estimation method, this is the distribution that provides the best fit to the given values).
}
\note{
   \bold{A word of caution:} Under the given distributional assumptions, the estimated means and SDs are approximately unbiased and hence so are any effect size measures computed based on them (assuming a measure is unbiased to begin with when computed with directly reported means and SDs). However, the estimated means and SDs are less precise (i.e., are more variable) than directly reported means and SDs (especially under case 1) and hence computing the sampling variance of a measure with equations that assume that directly reported means and SDs are available will tend to underestimate the actual sampling variance of the measure, giving too much weight to estimates computed based on estimated means and SDs (see also McGrath et al., 2023). It would therefore be prudent to treat effect size estimates computed from estimated means and SDs with caution (e.g., by examining in a moderator analysis whether there are systematic differences between studies directly reporting means and SDs and those where the means and SDs needed to be estimated and/or as part of a sensitivity analysis). McGrath et al. (2023) also suggest to use bootstrapping to estimate the sampling variance of effect size measures computed based on estimated means and SDs. See also the \href{https://cran.r-project.org/package=metamedian}{metamedian} package for this purpose.

   Also note that the development of methods for estimating means and SDs based on five-number summary values is an active area of research. Currently, when \code{method="default"}, then this is identical to \code{method="luo/wan/shi"}, but this might change in the future. For reproducibility, it is therefore recommended to explicitly set \code{method="luo/wan/shi"} (or one of the other methods) when running this function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Bland, M. (2015). Estimating mean and standard deviation from the sample size, three quartiles, minimum, and maximum. \emph{International Journal of Statistics in Medical Research}, \bold{4}(1), 57--64. \verb{https://doi.org/10.6000/1929-6029.2015.04.01.6}

   Cai, S., Zhou, J., & Pan, J. (2021). Estimating the sample mean and standard deviation from order statistics and sample size in meta-analysis. \emph{Statistical Methods in Medical Research}, \bold{30}(12), 2701--2719. \verb{https://doi.org/10.1177/09622802211047348}

   Hozo, S. P., Djulbegovic, B. & Hozo, I. (2005). Estimating the mean and variance from the median, range, and the size of a sample. \emph{BMC Medical Research Methodology}, \bold{5}, 13. \verb{https://doi.org/10.1186/1471-2288-5-13}

   Luo, D., Wan, X., Liu, J. & Tong, T. (2016). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical Methods in Medical Research}, \bold{27}(6), 1785--1805. \verb{https://doi.org/10.1177/0962280216669183}

   McGrath, S., Zhao, X., Steele, R., Thombs, B. D., Benedetti, A., & the DEPRESsion Screening Data (DEPRESSD) Collaboration (2020). Estimating the sample mean and standard deviation from commonly reported quantiles in meta-analysis. \emph{Statistical Methods in Medical Research}, \bold{29}(9), 2520--2537. \verb{https://doi.org/10.1177/0962280219889080}

   McGrath, S., Katzenschlager, S., Zimmer, A. J., Seitel, A., Steele, R., & Benedetti, A. (2023). Standard error estimation in meta-analysis of studies reporting medians. \emph{Statistical Methods in Medical Research}, \bold{32}(2), 373--388. \verb{https://doi.org/10.1177/09622802221139233}

   Shi, J., Luo, D., Weng, H., Zeng, X.-T., Lin, L., Chu, H. & Tong, T. (2020a). Optimally estimating the sample standard deviation from the five-number summary. \emph{Research Synthesis Methods}, \bold{11}(5), 641--654. \verb{https://doi.org/https://doi.org/10.1002/jrsm.1429}

   Shi, J., Tong, T., Wang, Y. & Genton, M. G. (2020b). Estimating the mean and variance from the five-number summary of a log-normal distribution. \emph{Statistics and Its Interface}, \bold{13}(4), 519--531. https://doi.org/10.4310/sii.2020.v13.n4.a9

   Shi, J., Luo, D., Wan, X., Liu, Y., Liu, J., Bian, Z. & Tong, T. (2023). Detecting the skewness of data from the five-number summary and its application in meta-analysis. \emph{Statistical Methods in Medical Research}, \bold{32}(7), 1338--1360. \verb{https://doi.org/10.1177/09622802231172043}

   Walter, S. D. & Yao, X. (2007). Effect sizes can be calculated for studies reporting ranges for outcome variables in systematic reviews. \emph{Journal of Clinical Epidemiology}, \bold{60}(8), 849-852. \verb{https://doi.org/10.1016/j.jclinepi.2006.11.003}

   Wan, X., Wang, W., Liu, J. & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC Medical Research Methodology}, \bold{14}, 135. \verb{https://doi.org/10.1186/1471-2288-14-135}

   Yang, X., Hutson, A. D., & Wang, D. (2022). A generalized BLUE approach for combining location and scale information in a meta-analysis. \emph{Journal of Applied Statistics}, \bold{49}(15), 3846--3867. \verb{https://doi.org/10.1080/02664763.2021.1967890}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to compute various effect size measures based on means and standard deviations.
}
\examples{
# example data frame
dat <- data.frame(case=c(1:3,NA), min=c(2,NA,2,NA), q1=c(NA,4,4,NA),
                  median=c(6,6,6,NA), q3=c(NA,10,10,NA), max=c(14,NA,14,NA),
                  mean=c(NA,NA,NA,7.0), sd=c(NA,NA,NA,4.2), n=c(20,20,20,20))
dat

# note that study 4 provides the mean and SD directly, while studies 1-3 provide five-number
# summary values or a subset thereof (corresponding to cases 1-3 above)

# estimate means/SDs (note: existing values in 'mean' and 'sd' are not touched)
dat <- conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat)
dat

# check attributes (none of the tests are significant, so means/SDs are estimated for studies 1-3)
dfround(data.frame(attributes(dat$mean)), digits=3)

# calculate the log transformed coefficient of variation and corresponding sampling variance
dat <- escalc(measure="CVLN", mi=mean, sdi=sd, ni=n, data=dat)
dat

# fit equal-effects model to the estimates
res <- rma(yi, vi, data=dat, method="EE")
res

# estimated coefficient of variation (with 95\% CI)
predict(res, transf=exp, digits=2)

############################################################################

# example data frame
dat <- data.frame(case=c(1:3,NA), min=c(2,NA,2,NA), q1=c(NA,4,4,NA),
                  median=c(6,6,6,NA), q3=c(NA,10,10,NA), max=c(14,NA,14,NA),
                  mean=c(NA,NA,NA,7.0), sd=c(NA,NA,NA,4.2), n=c(20,20,20,20))
dat

# try out different methods
conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat)
set.seed(1234)
conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat, method="qe")
conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat, method="bc")
conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat, method="mln")
conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat, method="blue")

############################################################################

# example data frame
dat <- data.frame(case=c(1:3,NA), min=c(2,NA,2,NA), q1=c(NA,4,4,NA),
                  median=c(6,6,6,NA), q3=c(NA,10,14,NA), max=c(14,NA,20,NA),
                  mean=c(NA,NA,NA,7.0), sd=c(NA,NA,NA,4.2), n=c(20,20,20,20))
dat

# for study 3, the third quartile and maximum value suggest that the data have
# a right skewed distribution (they are much further away from the median than
# the minimum and first quartile)

# estimate means/SDs
dat <- conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat)
dat

# note that the mean and SD are not estimated for study 3; this is because the
# test for skewness is significant for this study
dfround(data.frame(attributes(dat$mean)), digits=3)

# estimate means/SDs, but assume that the data for study 3 come from a log-normal distribution
# and back-transform the estimated mean/SD of the log-transformed data back to the raw data
dat <- conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat,
                    dist=c("norm","norm","lnorm","norm"), replace="all")
dat

# this works now because the test for skewness of the log-transformed data is not significant
dfround(data.frame(attributes(dat$mean)), digits=3)
}
\keyword{manip}
