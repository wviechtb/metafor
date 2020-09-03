\name{ranktest}
\alias{ranktest}
\alias{ranktest.rma}
\alias{ranktest.default}
\title{Rank Correlation Test for Funnel Plot Asymmetry}
\description{
   The function can be used to carry out the rank correlation test for funnel plot asymmetry.
}
\usage{
ranktest(x, \dots)

\method{ranktest}{rma}(x, digits, \dots)

\method{ranktest}{default}(x, vi, sei, subset, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"} or a vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances (needed if \code{x} is a vector with the observed effect sizes or outcomes).}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{subset}{optional (logical or numeric) vector indicating the subset of studies that should be included in the test. Only relevant when passing a vector via \code{x}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is 4).}
   \item{\dots}{other arguments.}
}
\details{
   The function carries out the rank correlation test as described by Begg and Mazumdar (1994). The test can be used to examine whether the observed effect sizes or outcomes and the corresponding sampling variances are correlated. A high correlation would indicate that the funnel plot is asymmetric, which may be a result of publication bias.

   One can either pass an object of class \code{"rma"} to the function or a vector with the observed effect sizes or outcomes (via \code{x}) and the corresponding sampling variances via \code{vi} (or the standard errors via \code{sei}).
}
\value{
   An object of class \code{"ranktest"}. The object is a list containing the following components:
   \item{tau}{the estimated value of Kendall's tau rank correlation coefficient}
   \item{pval}{the corresponding p-value for the test that the true tau is equal to zero}

   The results are formatted and printed with the \code{\link{print.ranktest}} function.
}
\note{
   The method does not depend on the model fitted. Therefore, regardless of the model passed to the function, the results of the rank test will always be the same. See \code{\link{regtest}} for tests of funnel plot asymmetry that are based on regression models and model dependent.

   The function makes use of the \code{\link{cor.test}} function with \code{method="kendall"}. If possible, an exact p-value is provided; otherwise, a large-sample approximation is used.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Begg, C. B., & Mazumdar, M. (1994). Operating characteristics of a rank correlation test for publication bias. \emph{Biometrics}, \bold{50}(4), 1088--1101. \verb{https://doi.org/10.2307/2533446}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{regtest}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

ranktest(res)
ranktest(dat$yi, dat$vi)
}
\keyword{htest}
