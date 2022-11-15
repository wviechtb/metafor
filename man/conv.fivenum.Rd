\name{conv.fivenum}
\alias{conv.fivenum}
\title{Convert Five-Number Summary Values to Means and Standard Deviations}
\description{
   Function to convert five-number summary values (or a subset thereof) to means and standard deviations.
}
\usage{
conv.fivenum(min, q1, median, q3, max, n, data, include,
             var.names=c("mean","sd"), append=TRUE, method)
}
\arguments{
   \item{min}{vector with the sample minimum values.}
   \item{q1}{vector with the lower/first quartile values.}
   \item{median}{vector with the median values.}
   \item{q3}{vector with the upper/third quartile values.}
   \item{max}{vector with the sample maximum values.}
   \item{n}{vector with the sample sizes.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{include}{optional (logical or numeric) vector to specify the subset of studies for which the conversion should be carried out (if possible).}
   \item{var.names}{character vector with two elements to specify the name of the variable for the means and the name of the variable for the standard deviations (the defaults are \code{"mean"} and \code{"sd"}).}
   \item{append}{logical to specify whether the data frame provided via the \code{data} argument should be returned together with the estimated values (the default is \code{TRUE}).}
   \item{method}{optional character string (or vector with two values) specifying the method for estimating the means and standard deviations, respectively. By default, some recommended methods are used (see below for details).}
}
\details{
   Various effect size measures require means and standard deviations (SDs) as input (e.g., raw or standardized mean differences, ratios of means / response ratios; see \code{\link{escalc}} for further details). For some studies, authors may not report means and SDs, but other statistics, such as the so-called \sQuote{five-number summary}, consisting of the minimum, lower/first quartile, median, upper/third quartile, and the maximum of the sample values (plus the sample sizes). Occasionally, only a subset of these values are reported.

   The present function can be used to convert five-number summary values (or a subset thereof) to means and standard deviations. To be precise, the function provides estimates of the means and standard deviations based on methods described in the literature (Hozo et al., 2005; Wan et al., 2014; Bland, 2015; Luo et al., 2016; Shi et al., 2020).

   We can distinguish three cases:

   \subsection{Case 1: Min, Median, Max}{

      In case only the minimum, median, and maximum is available for a study (plus the sample size), then the function uses by default the method by Luo et al. (2016), equation (7), to estimate the mean and the method by Wan et al. (2014), equation (9), to estimate the SD.

   }

   \subsection{Case 2: Q1, Median, Q3}{

      In case only the lower/first quartile, median, and upper/third quartile is available for a study (plus the sample size), then the function uses by default the method by Luo et al. (2016), equation (11), to estimate the mean and the method by Wan et al. (2014), equation (16), to estimate the SD.

   }

   \subsection{Case 3: Min, Q1, Median, Q3, Max}{

      In case the full five-number summary is available for a study (plus the sample size), then the function uses by default the method by Luo et al. (2016), equation (15), to estimate the mean and the method by Shi et al. (2020), equation (10), to estimate the SD.

   }

   Note that the sample size must be at least 5 to apply these methods. Studies where the sample size is smaller are not included in the conversion. The function also checks that \code{min <= q1 <= median <= q3 <= max} and throws an error if any studies are found where this is not the case.

}
\value{
   If \code{data} was not specified or \code{append=FALSE}, a data frame with two variables called \code{var.names[1]} and \code{var.names[2]} with the estimated means and SDs.

   If \code{data} was specified and \code{append=TRUE}, then the original data frame is returned. If \code{var.names[1]} is a variable in \code{data}, then missing values in this variable are replaced with the estimated means (where possible) and otherwise a new variable called \code{var.names[1]} is added to the data frame. Similarly, if \code{var.names[2]} is a variable in \code{data}, then missing values in this variable are replaced with the estimated SDs (where possible) and otherwise a new variable called \code{var.names[2]} is added to the data frame.

   When missing values in \code{var.names[1]} are replaced, an attribute called \code{"est"} is added to the variable, which is a logical vector that is \code{TRUE} for values that were estimated. The same is done when missing values in \code{var.names[2]} are replaced.
}
\note{
   A few alternative methods are also implemented, but these are not recommended as they have been shown to be less accurate than the methods that are used by default.

   For case 1, one can set \code{method} to \code{"Hozo2005"} in which case the methods as described on page 8, right column, from Hozo et al. (2005) are used to estimate the mean and SD. One can set \code{method} to \code{"Wan2014"} in which case equation (3) from Wan et al. (2014) is used to estimate the mean (and equation (9) for the SD per default).

   For case 2, one can set \code{method} to \code{"Wan2014"} in which case equation (14) from Wan et al. (2014) is used to estimate the mean (and equation (16) for the SD per default).

   For case 3, one can set \code{method[1]} to \code{"Bland2015"} or \code{"Wan2014"}, in which case equation (10) from Wan et al. (2014) is used to estimate the mean. One can set \code{method[2]} to \code{"Bland2015"} in which case equation (11) from Wan et al. (2014) is used to estimate the SD or to \code{"Wan2014"} in which case equation (13) from Wan et al. (2014) is used to estimate the SD.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Bland, M. (2015). Estimating mean and standard deviation from the sample size, three quartiles, minimum, and maximum. \emph{International Journal of Statistics in Medical Research}, \bold{4}(1), 57--64. \verb{https://doi.org/10.6000/1929-6029.2015.04.01.6}

   Hozo, S. P., Djulbegovic, B. & Hozo, I. (2005). Estimating the mean and variance from the median, range, and the size of a sample. \emph{BMC Medical Research Methodology}, \bold{5}, 13. \verb{https://doi.org/10.1186/1471-2288-5-13}

   Luo, D., Wan, X., Liu, J. & Tong, T. (2016). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical Methods in Medical Research}, \bold{27}(6), 1785--1805. \verb{https://doi.org/10.1177/0962280216669183}

   Shi, J., Luo, D., Weng, H., Zeng, X.-T., Lin, L., Chu, H. & Tong, T. (2020). Optimally estimating the sample standard deviation from the five-number summary. \emph{Research Synthesis Methods}, \bold{11}(5), 641--654. \verb{https://doi.org/https://doi.org/10.1002/jrsm.1429}

   Wan, X., Wang, W., Liu, J. & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC Medical Research Methodology}, \bold{14}, 135. \verb{https://doi.org/10.1186/1471-2288-14-135}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to compute various effect size measures based on means and standard deviations.
}
\examples{
# example data frame
dat <- data.frame(case=c(1:3,NA), min=c(2,NA,2,NA), q1=c(NA,4,4,NA),
                  median=c(6,6,6,NA), q3=c(NA,10,10,NA), max=c(16,NA,16,NA),
                  mean=c(NA,NA,NA,7.0), sd=c(NA,NA,NA,4.2), n=c(20,20,20,20))
dat

# estimate means/SDs (note: existing values in 'mean' and 'sd' are not touched)
dat <- conv.fivenum(min=min, q1=q1, median=median, q3=q3, max=max, n=n, data=dat)
dat

# check for 'est' attribute
attributes(dat$mean)
attributes(dat$sd)

# calculate the log transformed coefficient of variation and corresponding sampling variance
dat <- escalc(measure="CVLN", mi=mean, sdi=sd, ni=n, data=dat)
dat

# fit equal-effects model to the estimates
res <- rma(yi, vi, data=dat, method="EE")
res

# estimated coefficient of variation (with 95% CI)
predict(res, transf=exp, digits=2)
}
\keyword{manip}
