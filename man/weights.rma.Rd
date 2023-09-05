\name{weights.rma}
\alias{weights}
\alias{weights.rma}
\alias{weights.rma.uni}
\alias{weights.rma.mh}
\alias{weights.rma.peto}
\alias{weights.rma.glmm}
\alias{weights.rma.mv}
\title{Compute Weights for 'rma' Objects}
\description{
   Functions to compute the weights given to the observed effect sizes or outcomes during the model fitting for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, and \code{"rma.mv"}.
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
   \item{type}{character string to specify whether to return only the diagonal of the weight matrix (\code{"diagonal"}) or the entire weight matrix (\code{"matrix"}). For \code{"rma.mv"}, this can also be \code{"rowsum"} for \sQuote{row-sum weights} (for intercept-only models).}
   \item{\dots}{other arguments.}
}
\value{
   Either a vector with the diagonal elements of the weight matrix or the entire weight matrix. When only the diagonal elements are returned, they are given in \% (and they add up to 100\%).

   When the entire weight matrix is requested, this is always a diagonal matrix for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}.

   For \code{"rma.mv"}, the structure of the weight matrix depends on the model fitted (i.e., the random effects included and the variance-covariance matrix of the sampling errors) but is often more complex and not just diagonal.

   For intercept-only \code{"rma.mv"} models, one can also take the sum over the rows in the weight matrix, which are actually the weights assigned to the observed effect sizes or outcomes when estimating the model intercept. These weights can be obtained with \code{type="rowsum"} (as with \code{type="diagonal"}, they are also given in \%). See \href{https://www.metafor-project.org/doku.php/tips:weights_in_rma.mv_models}{here} for a discussion of this.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W. (2021). Model checking in meta-analysis. In C. H. Schmid, T. Stijnen, & I. R. White (Eds.), \emph{Handbook of meta-analysis} (pp. 219--254). Boca Raton, FL: CRC Press. \verb{https://doi.org/10.1201/9781315119403}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.mv}} for functions to fit models for which model fitting weights can be extracted.

   \code{\link{influence.rma.uni}} and \code{\link{influence.rma.mv}} for other model diagnostics.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### extract the model fitting weights (in %)
weights(res)

### extract the weight matrix
round(weights(res, type="matrix"), 4)
}
\keyword{models}
