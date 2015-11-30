\name{print.confint.rma}
\alias{print.confint.rma}
\title{Print Method for 'confint.rma' Objects}
\description{
   Print method for objects of class \code{"confint.rma"}.
}
\usage{
\method{print}{confint.rma}(x, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"confint.rma"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item estimate of the model coefficient or variance/correlation parameter
   \item lower bound of the confidence interval
   \item upper bound of the confidence interval
   }
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{confint.rma.uni}}
}
\keyword{print}
