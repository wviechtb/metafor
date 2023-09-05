\name{print.regtest}
\alias{print.regtest}
\title{Print Method for 'regtest' Objects}
\description{
   Function to print objects of class \code{"regtest"}.
}
\usage{
\method{print}{regtest}(x, digits=x$digits, ret.fit=x$ret.fit, \dots)
}
\arguments{
   \item{x}{an object of class \code{"regtest"} obtained with \code{\link{regtest}}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{ret.fit}{logical to specify whether the full results from the fitted model should also be returned. If unspecified, the default is to take the value from the object.}
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
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{regtest}} for the function to create \code{regtest} objects.
}
\keyword{print}
