\name{methods.anova.rma}
\alias{methods.anova.rma}
\alias{as.data.frame.anova.rma}
\alias{as.data.frame.list.anova.rma}
\title{Methods for 'anova.rma' Objects}
\description{
   Methods for objects of class \code{"anova.rma"} and \code{"list.anova.rma"}.
}
\usage{
\method{as.data.frame}{anova.rma}(x, \dots)
\method{as.data.frame}{list.anova.rma}(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"anova.rma"} or \code{"list.anova.rma"}.}
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

### fit mixed-effects meta-regression model
res <- rma(yi, vi, mods = ~ alloc + ablat, data=dat)

### test the allocation factor
sav <- anova(res, btt="alloc")
sav

### turn object into a regular data frame
as.data.frame(sav)

### test the contrast between levels random and systematic
sav <- anova(res, X=c(0,1,-1,0))
sav

### turn object into a regular data frame
as.data.frame(sav)

### fit random-effects model
res0 <- rma(yi, vi, data=dat)

### LRT comparing the two models
sav <- anova(res, res0, refit=TRUE)
sav

### turn object into a regular data frame
as.data.frame(sav)
}
\keyword{internal}
