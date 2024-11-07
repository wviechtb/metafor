\name{plot.rma}
\alias{plot.rma}
\alias{plot.rma.uni}
\alias{plot.rma.mh}
\alias{plot.rma.mv}
\alias{plot.rma.peto}
\alias{plot.rma.glmm}
\title{Plot Method for 'rma' Objects}
\description{
   Functions to plot objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, and \code{"rma.glmm"}.
}
\usage{
\method{plot}{rma.uni}(x, qqplot=FALSE, \dots)

\method{plot}{rma.mh}(x, qqplot=FALSE, \dots)

\method{plot}{rma.peto}(x, qqplot=FALSE, \dots)

\method{plot}{rma.glmm}(x, qqplot=FALSE, \dots) # not currently implemented

\method{plot}{rma.mv}(x, qqplot=FALSE, \dots) # not currently implemented
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}. The method is not (yet) implemented for objects of class \code{"rma.glmm"} or \code{"rma.mv"}.}
   \item{qqplot}{logical to specify whether a normal QQ plot should be drawn (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   Four plots are produced. If the model does not contain any moderators, then a forest plot, funnel plot, radial plot, and a plot of the standardized residuals is provided. If \code{qqplot=TRUE}, the last plot is replaced by a normal QQ plot of the standardized residuals.

   If the model contains moderators, then a forest plot, funnel plot, plot of the standardized residuals against the fitted values, and a plot of the standardized residuals is provided. If \code{qqplot=TRUE}, the last plot is replaced by a normal QQ plot of the standardized residuals.
}
\note{
   If the number of studies is large, the forest plot may become difficult to read due to the small font size. Stretching the plotting device vertically should provide more space.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest}} for forest plots, \code{\link{funnel}} for funnel plots, \code{\link{radial}} for radial plots, and \code{\link[=qqnorm.rma]{qqnorm}} for normal QQ plots.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model
res <- rma(yi, vi, data=dat)

### plot results
plot(res, qqplot=TRUE)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### plot results
plot(res, qqplot=TRUE)
}
\keyword{hplot}
