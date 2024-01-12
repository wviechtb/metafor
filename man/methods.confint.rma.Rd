\name{methods.confint.rma}
\alias{methods.confint.rma}
\alias{as.data.frame.confint.rma}
\alias{as.data.frame.list.confint.rma}
\title{Methods for 'confint.rma' Objects}
\description{
   Methods for objects of class \code{"confint.rma"} and \code{"list.confint.rma"}.
}
\usage{
\method{as.data.frame}{confint.rma}(x, \dots)
\method{as.data.frame}{list.confint.rma}(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"confint.rma"} or \code{"list.confint.rma"}.}
   \item{\dots}{other arguments.}
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
### copy data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### fit random-effects model
res <- rma(yi, vi, data=dat)

### get 95\% CI for tau^2, tau, I^2, and H^2
sav <- confint(res)
sav

### turn object into a regular data frame
as.data.frame(sav)

############################################################################

### copy data into 'dat'
dat <- dat.berkey1998

### construct block diagonal var-cov matrix of the observed outcomes based on variables v1i and v2i
V <- vcalc(vi=1, cluster=author, rvars=c(v1i, v2i), data=dat)

### fit multivariate model
res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat)

### get 95\% CI for variance components and correlation
sav <- confint(res)
sav

### turn object into a regular data frame
as.data.frame(sav)
}
\keyword{internal}
