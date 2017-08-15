\name{print.regtest.rma}
\alias{print.regtest.rma}
\title{Print Method for 'regtest.rma' Objects}
\description{
   Print method for objects of class \code{"regtest.rma"}.
}
\usage{
\method{print}{regtest.rma}(x, digits, ret.fit, \dots)
}
\arguments{
   \item{x}{an object of class \code{"regtest.rma"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{ret.fit}{logical indicating whether the full results from the fitted model should also be returned (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item the model used for the regression test
   \item the predictor used for the regression test
   \item the results from the fitted model (only when \code{ret.fit=TRUE})
   \item the value of the test statistic for the test that the predictor is unreleated to the outcomes
   \item the degrees of freedom of the test statistic (only if the test statistic follows a t-distribution)
   \item the p-value for the test statistic
   }
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{regtest.rma}}
}
\keyword{print}
