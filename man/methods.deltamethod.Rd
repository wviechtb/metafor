\name{coef.deltamethod}
\alias{coef.deltamethod}
\alias{vcov.deltamethod}
\title{Extract the Estimates and Variance-Covariance Matrix from 'deltamethod' Objects}
\description{
   Methods for objects of class \code{"deltamethod"}.
}
\usage{
\method{coef}{deltamethod}(object, \dots)
\method{vcov}{deltamethod}(object, \dots)
}
\arguments{
   \item{object}{an object of class \code{"deltamethod"}.}
   \item{\dots}{other arguments.}
}
\details{
   The \code{coef} function extracts the transformed estimates from objects of class \code{"deltamethod"}. The \code{vcov} function extracts the corresponding variance-covariance matrix.
}
\value{
   Either a vector with the transformed estimates or a variance-covariance matrix.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{deltamethod}} for the function to create \code{deltamethod} objects.
}
\keyword{models}
