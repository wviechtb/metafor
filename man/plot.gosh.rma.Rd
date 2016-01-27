\name{plot.gosh.rma}
\alias{plot.gosh.rma}
\title{Plot Method for 'gosh.rma' Objects}
\description{
   Plot method for objects of class \code{"gosh.rma"}.
}
\usage{
\method{plot}{gosh.rma}(x, het="I2", pch=19, cex, out, col,
     alpha, breaks, adjust, lwd, labels, \ldots)
}
\arguments{
   \item{x}{an object of class \code{"gosh.rma"}.}
   \item{het}{character string indicating the heterogeneity measure to plot. Either \code{"I2"}, \code{"H2"}, \code{"QE"}, or \code{"tau2"} (the last only for random/mixed-effects models).}
   \item{pch}{plotting symbol to use. By default, a filled circle is used. See \code{\link{points}} for other options.}
   \item{cex}{symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{out}{optional integer indicating the number of a study that may be a potential outlier. If specified, subsets containing the specified study are drawn in a different color than those not containing the study.}
   \item{col}{optional character string specifying the name of a color to use for the points (if not provided, points are drawn in black). When \code{out} is used, two colors should be specified (if not provided, red is used for subsets containing the specified study and blue otherwise).}
   \item{alpha}{optional alpha transparency value for the points (0 means fully transparent and 1 means opaque). If unspecified, the function tries to set this to a sensible value. Ignored when \code{col} is specified.}
   \item{breaks}{optional argument passed on to \code{\link{hist}} for choosing the (number of) breakpoints of the histograms.}
   \item{adjust}{optional argument passed on to \code{\link{density}} for adjusting the bandwidth of the kernel density estimates.}
   \item{lwd}{optional numeric value to specify the line width of the kernel density estimates.}
   \item{labels}{optional argument passed on to \code{\link{pairs}} for specifying the names of the variables.}
   \item{\dots}{other arguments.}
}
\details{
   The function draws a scatterplot matrix with the \code{\link{pairs}} function. Histograms of the variables plotted are shown along the diagonal, with a kernel density estimate of the distribution superimposed.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Olkin, I., Dahabreh, I. J., & Trikalinos, T. A. (2012). GOSH - a graphical display of study heterogeneity. \emph{Research Synthesis Methods}, \bold{3}, 214--223.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{gosh.rma}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fixed-effects model
res <- rma(yi, vi, data=dat, method="FE")
res

### fit all possible subsets
sav <- gosh(res)

### plot results
plot(sav, breaks=80)

### plot results (red for subsets containing study 8; blue otherwise)
plot(sav, breaks=80, out=8)
}
\keyword{hplot}
