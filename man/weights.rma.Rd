\name{weights.rma.uni}
\alias{weights.rma.uni}
\alias{weights.rma.mh}
\alias{weights.rma.peto}
\alias{weights.rma.glmm}
\alias{weights.rma.mv}
\title{Compute Weights for 'rma' Objects}
\description{
   The function computes the weights given to the observed effects or outcomes during the model fitting for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, and \code{"rma.mv"}.
}
\usage{
\method{weights}{rma.uni}(object, type="diagonal", \dots)
\method{weights}{rma.mh}(object, type="diagonal", \dots)
\method{weights}{rma.peto}(object, type="diagonal", \dots)
\method{weights}{rma.glmm}(object, \dots)
\method{weights}{rma.mv}(object, type="diagonal", \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, or \code{"rma.mv"}. The method is not yet implemented for objects of class \code{"rma.glmm"}.}
   \item{type}{character string indicating whether to return only the diagonal of the weight matrix (\code{"diagonal"}) or the entire weight matrix (\code{"matrix"}).}
   \item{\dots}{other arguments.}
}
\value{
   Either a vector with the diagonal elements of the weight matrix or the entire weight matrix. When only the diagonal elements are returned, they are given in \% (and they add up to 100\%). When the entire weight matrix is requested, this is always a diagonal matrix for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}.
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
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.mv}}, \code{\link{influence.rma.uni}}
}
\examples{
### meta-analysis of the log risk ratios using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, 
           mods = ~ ablat + year, data=dat.bcg)
weights(res)
}
\keyword{models}
