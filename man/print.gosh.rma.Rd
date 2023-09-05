\name{print.gosh.rma}
\alias{print.gosh.rma}
\title{Print Method for 'gosh.rma' Objects}
\description{
   Function to print objects of class \code{"gosh.rma"}.
}
\usage{
\method{print}{gosh.rma}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"gosh.rma"} obtained with \code{\link{gosh}}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output shows how many model fits were attempted, how many succeeded, and summary statistics (i.e., the mean, minimum, first quartile, median, third quartile, and maximum) for the various measures of (residual) heterogeneity and the model coefficient(s) computed across all of the subsets.
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
   \code{\link{gosh}} for the function to create \code{gosh.rma} objects.
}
\keyword{print}
