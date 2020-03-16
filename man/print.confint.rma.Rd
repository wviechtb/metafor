\name{print.confint.rma}
\alias{print.confint.rma}
\alias{print.list.confint.rma}
\title{Print Methods for 'confint.rma' and 'list.confint.rma' Objects}
\description{
   Print methods for objects of class \code{"confint.rma"} and \code{"list.confint.rma"}.
}
\usage{
\method{print}{confint.rma}(x, digits=x$digits, \dots)
\method{print}{list.confint.rma}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"confint.rma"} or \code{"list.confint.rma"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
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
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\seealso{
   \code{\link{confint.rma.uni}}, \code{\link{confint.rma.mv}}
}
\keyword{print}
