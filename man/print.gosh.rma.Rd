\name{print.gosh.rma}
\alias{print.gosh.rma}
\title{Print Method for 'gosh.rma' Objects}
\description{
   Print method for objects of class \code{"gosh.rma"}.
}
\usage{
\method{print}{gosh.rma}(x, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"gosh.rma"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output shows how many model fits were attempted, how many succeeded, and summary statistics (i.e., the mean and the five-number summary) for the various measures of (residual) heterogeneity and the model coefficient(s) computed across all of the subsets.
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
   \code{\link{gosh.rma}}
}
\keyword{print}
