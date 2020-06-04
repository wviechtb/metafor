\name{print.ranktest.rma}
\alias{print.ranktest.rma}
\title{Print Method for 'ranktest.rma' Objects}
\description{
   Print method for objects of class \code{"ranktest.rma"}.
}
\usage{
\method{print}{ranktest.rma}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"ranktest.rma"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item the estimated value of Kendall's tau rank correlation coefficient
   \item the corresponding p-value for the test that the true tau is equal to zero
   }
}
\value{
   The function does not return an object.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\seealso{
   \code{\link{ranktest.rma}}
}
\keyword{print}
