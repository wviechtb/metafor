\name{model.matrix.rma}
\alias{model.matrix}
\alias{model.matrix.rma}
\title{Extract the Model Matrix from 'rma' Objects}
\description{
   Function to extract the model matrix from objects of class \code{"rma"}.
}
\usage{
\method{model.matrix}{rma}(object, asdf, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{asdf}{logical to specify whether the model matrix should be turned into a data frame (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\value{
   The model matrix.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which a model matrix can be extracted.

   \code{\link[=fitted.rma]{fitted}} for a function to extract the fitted values.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### extract the model matrix
model.matrix(res)
}
\keyword{models}
