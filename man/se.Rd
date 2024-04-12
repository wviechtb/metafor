\name{se}
\alias{se}
\alias{se.default}
\alias{se.rma}
\title{Extract the Standard Errors from 'rma' Objects}
\description{
   Function to extract the standard errors from objects of class \code{"rma"}.
}
\usage{
se(object, \dots)
\method{se}{default}(object, \dots)
\method{se}{rma}(object, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{\dots}{other arguments.}
}
\value{
   A vector with the standard errors.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which standard errors can be extracted.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)
res

### extract model coefficients
coef(res)

### extract the standard errors
se(res)
}
\keyword{models}
