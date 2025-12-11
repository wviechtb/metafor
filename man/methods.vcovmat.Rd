\name{methods.vcovmat}
\alias{methods.vcovmat}
\alias{print.vcovmat}
\alias{[.vcovmat}
\title{Methods for 'vcovmat' Objects}
\description{
   Methods for objects of class \code{"vcovmat"}.
}
\usage{
\method{print}{vcovmat}(x, digits=4, tol, zero=".", na="", \dots)

\method{[}{vcovmat}(x, i, j, \dots)
}
\arguments{
   \item{x}{an object of class \code{"vcovmat"}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (the default is 4).}
   \item{tol}{numeric value giving the tolerance for values that will be considered to be zeros.}
   \item{zero}{character string to represent zero values.}
   \item{na}{character string to represent missing values.}
   \item{i,j}{indices to select rows/columns.}
   \item{\dots}{other arguments.}
}
\details{
   When printing a \code{"vcovmat"} object, values equal to zero are printed by default as a period. This makes it easier to see the structure of the variance-covariance matrix, which often has a block-diagonal structure. Values sufficiently close to zero are also treated as zero. This is controlled by the \code{tol} argument, which, if unspecified, is set by default to \code{10 * .Machine$double.eps}.
}
\note{
   To turn a matrix \code{x} of class \code{"vcovmat"} into a regular matrix, just use \code{unclass(x)}.
}
\author{
   The \code{print} method is based on the code from \code{\link{print.table}} with some minor tweaks.
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{vcalc}} and \code{\link{rcalc}} for functions that create \code{vcovmat} objects.
}
\keyword{print}
