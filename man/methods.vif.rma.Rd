\name{methods.vif.rma}
\alias{methods.vif.rma}
\alias{as.data.frame.vif.rma}
\title{Methods for 'vif.rma' Objects}
\description{
   Methods for objects of class \code{"vif.rma"}.
}
\usage{
\method{as.data.frame}{vif.rma}(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"vif.rma"}.}
   \item{\dots}{other arguments.}
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
### copy data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### fit mixed-effects meta-regression model
res <- rma(yi, vi, mods = ~ ablat + year + alloc, data=dat)

### get variance inflation factors for all individual coefficients
sav <- vif(res)
sav

### turn object into a regular data frame
as.data.frame(sav)

### get VIFs for ablat and year and the generalized VIF for alloc
sav <- vif(res, btt=list("ablat","alloc","year"))
sav

### turn object into a regular data frame
as.data.frame(sav)
}
\keyword{internal}
