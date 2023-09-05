\name{print.fsn}
\alias{print.fsn}
\title{Print Method for 'fsn' Objects}
\description{
   Function to print objects of class \code{"fsn"}.
}
\usage{
\method{print}{fsn}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"fsn"} obtained with \code{\link{fsn}}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The output shows the results from the fail-safe N calculation.
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
   \code{\link{fsn}} for the function to create \code{fsn} objects.
}
\keyword{print}
