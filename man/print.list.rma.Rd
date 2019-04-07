\name{print.list.rma}
\alias{print.list.rma}
\title{Print method for 'list.rma' Objects}
\description{
   Print method for objects of class \code{"list.rma"}.
}
\usage{
\method{print}{list.rma}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"list.rma"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\value{
   See the documentation of the function that creates the \code{"list.rma"} object for details on what is printed. Regardless of what is printed, a data frame with the results is also returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\keyword{print}
