\name{print.regtest}
\alias{print.regtest}
\title{Print Method for 'regtest' Objects}
\description{
   Print method for objects of class \code{"regtest"}.
}
\usage{
\method{print}{regtest}(x, digits=x$digits, ret.fit=x$ret.fit, \dots)
}
\arguments{
   \item{x}{an object of class \code{"regtest"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{ret.fit}{logical indicating whether the full results from the fitted model should also be returned (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item the model used for the regression test
   \item the predictor used for the regression test
   \item the results from the fitted model (only when \code{ret.fit=TRUE})
   \item the test statistic of the test that the predictor is unreleated to the outcomes
   \item the degrees of freedom of the test statistic (only if the test statistic follows a t-distribution)
   \item the corresponding p-value
   \item the \sQuote{limit estimate} and its corresponding CI (only for predictors \code{"sei"} \code{"vi"}, \code{"ninv"}, or \code{"sqrtninv"} and when the model does not contain any additional moderators)
   }
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{regtest}}
}
\keyword{print}
