\name{update.rma}
\alias{update}
\alias{update.rma}
\title{Model Updating for 'rma' Objects}
\description{
   The function can be used to update and (by default) re-fit \code{"rma"} models. It does this by extracting the call stored in the object, updating the call and (by default) evaluating that call.
}
\usage{
\method{update}{rma}(object, formula., \dots, evaluate = TRUE)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{formula.}{changes to the formula. See \sQuote{Details}.}
   \item{\dots}{additional arguments to the call, or arguments with changed values.}
   \item{evaluate}{logical indicating whether to evaluate the new call or just return the call.}
}
\details{
   For objects of class \code{"rma.uni"}, \code{"rma.glmm"}, and \code{"rma.mv"}, the \code{formula.} argument can be used to update the set of moderators included in the model (see \sQuote{Examples}).
}
\value{
   If \code{evaluate=TRUE} the fitted object, otherwise the updated call.
}
\author{
   The present function is based on \code{\link{update.default}}, with changes made by Wolfgang Viechtbauer (\email{wvb@metafor-project.org}) so that the formula updating works with the (somewhat non-standard) interface of the \code{\link{rma.uni}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} functions.
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model (method="REML" is default)
res <- rma(yi, vi, data=dat, digits=3)
res

### mixed-effects model with two moderators (absolute latitude and publication year)
res <- update(res, ~ ablat + year)
res

### remove 'year' moderator
res <- update(res, ~ . - year)
res

### fit model with ML estimation
update(res, method="ML")

### example with rma.glmm()
res <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, digits=3)
res <- update(res, mods = ~ ablat)
res

### conditional model with approximate likelihood
update(res, model="CM.AL")
}
\keyword{models}
