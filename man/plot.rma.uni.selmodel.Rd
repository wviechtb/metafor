\name{plot.rma.uni.selmodel}
\alias{plot.rma.uni.selmodel}
\title{Plot Method for 'plot.rma.uni.selmodel' Objects}
\description{
   Plot method for objects of class \code{"plot.rma.uni.selmodel"}. \loadmathjax
}
\usage{
\method{plot}{rma.uni.selmodel}(x, xlim, n=1001, prec="max",
      scale=FALSE, rug=TRUE, add=FALSE, lwd=2, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}.}
   \item{xlim}{x-axis limits. Essentially the range of p-values for which the selection function should be plotted. If unspecified, the function sets limits automatically.}
   \item{n}{numerical value indicating for how many p-values within the x-axis limits the function value should be computed.}
   \item{prec}{either a character string (with options \code{"max"}, \code{"min"}, \code{"mean"}, or \code{"median"}) or a numeric value.}
   \item{scale}{logical specifying whether the function values should be rescaled to a 0 to 1 range (the default is \code{FALSE}).}
   \item{rug}{logical specifying whether the observed p-values should be added as tick marks on the x-axis (the default is \code{TRUE}).}
   \item{add}{logical specifying whether the function should be added to an existing plot (the default is \code{FALSE}).}
   \item{lwd}{numerical value for the line width.}
   \item{\dots}{other arguments.}
}
\details{
   Blah.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{selmodel.rma.uni}}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.hackshaw1998

### random-effects model using the log odds ratios
res <- rma(yi, vi, data=dat)
res

### fit step selection model
sel <- selmodel(res, type="stepfun", steps=c(0.05, 0.10, 0.50, 1.00))

### plot selection function
plot(sel, scale=TRUE)

### fit negative exponential selection model
sel <- selmodel(res, type="negexp")

### add selection function to the existing plot
plot(sel, add=TRUE, col="blue")
}
\keyword{hplot}
