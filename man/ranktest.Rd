\name{ranktest}
\alias{ranktest}
\alias{ranktest.rma}
\title{Rank Correlation Test for Funnel Plot Asymmetry for 'rma' Objects}
\description{Rank correlation test for funnel plot asymmetry for objects of class \code{"rma"}.}
\usage{
ranktest(x, \dots)

\method{ranktest}{rma}(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{\dots}{other arguments.}
}
\details{
   The function carries out the rank correlation test as described by Begg and Mazumdar (1994). The test can be used to examine whether the observed outcomes and the corresponding sampling variances are correlated. A high correlation would indicate that the funnel plot is asymmetric, which may be a result of publication bias.
}
\value{
   An object of class \code{"ranktest.rma"}. The object is a list containing the following components:
   \item{tau}{the estimated value of Kendall's tau rank correlation coefficient}
   \item{pval}{the corresponding p-value for the test that the true tau is equal to zero}

   The results are formated and printed with the \code{\link{print.ranktest.rma}} function.
}
\note{
   The method does not depend on the model fitted. Therefore, regardless of the model, the results of the rank test will always be the same. See \code{\link{regtest}} for tests of funnel plot asymmetry that are based on regression models and model dependent.

   The function makes use of the \code{\link{cor.test}} function with \code{method="kendall"}. If possible, an exact p-value is provided; otherwise, a large-sample approximation is used.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Begg, C. B., & Mazumdar, M. (1994). Operating characteristics of a rank correlation test for publication bias. \emph{Biometrics}, \bold{50}, 1088--1101.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{regtest}}
}
\examples{
### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

ranktest(res)
}
\keyword{htest}
