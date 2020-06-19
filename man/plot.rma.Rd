\name{plot.rma}
\alias{plot.rma}
\alias{plot.rma.uni}
\alias{plot.rma.mh}
\alias{plot.rma.peto}
\alias{plot.rma.glmm}
\title{Plot Method for 'rma' Objects}
\description{
   Plot method for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, and \code{"rma.glmm"}.
}
\usage{
\method{plot}{rma.uni}(x, qqplot=FALSE, \dots)

\method{plot}{rma.mh}(x, qqplot=FALSE, \dots)

\method{plot}{rma.peto}(x, qqplot=FALSE, \dots)

\method{plot}{rma.glmm}(x, qqplot=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}. The method is not yet implemented for objects of class \code{"rma.glmm"}.}
   \item{qqplot}{logical indicating whether a normal QQ plot should be drawn (the default is \code{FALSE}).}
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
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest}}, \code{\link{funnel}}, \code{\link{radial}}, \code{\link{qqnorm.rma.uni}}
}
\examples{
### meta-analysis of the log risk ratios using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
plot(res, qqplot=TRUE)

### mixed-effects model with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year, data=dat.bcg)
plot(res, qqplot=TRUE)
}
\keyword{hplot}
