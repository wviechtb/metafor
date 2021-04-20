\name{fitted.rma}
\alias{fitted}
\alias{fitted.rma}
\title{Fitted Values for 'rma' Objects}
\description{
   The function computes the fitted values for objects of class \code{"rma"}.
}
\usage{
\method{fitted}{rma}(object, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{\dots}{other arguments.}
}
\value{
   A vector with the fitted values.
}
\note{
   The \code{\link{predict.rma}} function also provides standard errors and confidence intervals for the fitted values. Best linear unbiased predictions (BLUPs) that combine the fitted values based on the fixed effects and the estimated contributions of the random effects can be obtained with \code{\link{blup.rma.uni}} (only for objects of class \code{"rma.uni"}).

   For objects not involving moderators, the fitted values are all identical to the estimated value of the model intercept.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{predict.rma}}, \code{\link{blup.rma.uni}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### compute the fitted values
fitted(res)
}
\keyword{models}
