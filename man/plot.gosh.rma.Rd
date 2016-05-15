\name{plot.gosh.rma}
\alias{plot.gosh.rma}
\title{Plot Method for 'gosh.rma' Objects}
\description{
   Plot method for objects of class \code{"gosh.rma"}.
}
\usage{
\method{plot}{gosh.rma}(x, het="I2", pch=16, cex, out, col, alpha,
     xlim, ylim, xhist=TRUE, yhist=TRUE, hh=.3, hcol, breaks,
     adjust, lwd, labels, \ldots)
}
\arguments{
   \item{x}{an object of class \code{"gosh.rma"}.}
   \item{het}{character string indicating the heterogeneity measure to plot. Either \code{"I2"}, \code{"H2"}, \code{"QE"}, or \code{"tau2"} (the last only for random/mixed-effects models).}
   \item{pch}{plotting symbol to use. By default, a borderless filled circle is used. See \code{\link{points}} for other options.}
   \item{cex}{symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{out}{optional integer indicating the number of a study that may be a potential outlier. If specified, subsets containing the specified study are drawn in a different color than those not containing the study.}
   \item{col}{optional character string specifying the name of a color to use for the points (if not provided, points are drawn in black). When \code{out} is used, two colors should be specified (if not provided, red is used for subsets containing the specified study and blue otherwise).}
   \item{alpha}{optional alpha transparency value for the points (0 means fully transparent and 1 means opaque). If unspecified, the function tries to set this to a sensible value. Ignored when \code{col} is specified.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xhist}{logical specifying whether a histogram should be drawn for the x-axis (the default is \code{TRUE}).}
   \item{yhist}{logical specifying whether a histogram should be drawn for the y-axis (the default is \code{TRUE}).}
   \item{hh}{optional numerical value (or vector of two values) to adjust the height of the histogram(s). Must be between 0 and 1, but should not be too close to 0 or 1, as otherwise the plot cannot be drawn.}
   \item{hcol}{optional character string specifying the name of a color to use for the histogram(s).}
   \item{breaks}{optional argument passed on to \code{\link{hist}} for choosing the (number of) breakpoints of the histogram(s).}
   \item{adjust}{optional argument passed on to \code{\link{density}} for adjusting the bandwidth of the kernel density estimation.}
   \item{lwd}{optional numeric value to specify the line width of the estimated densities. Set to \code{0} to omit the line(s).}
   \item{labels}{optional argument to specify the x-axis and y-axis labels (or passed on to \code{\link{pairs}} for specifying the names of the variables in the scatterplot matrix).}
   \item{\dots}{other arguments.}
}
\details{
   For models without moderators, the function draws a scatterplot of the model estimates on the x-axis against the chosen measure of heterogeneity on the y-axis. Histograms of the respective distributions (with kernel density estimates superimposed) are shown in the margins (when \code{xhist=TRUE} and \code{yhist=TRUE}).

   For models with moderators, the function draws a scatterplot matrix (with the \code{\link{pairs}} function) of the chosen measure of heterogeneity and each of the model coefficients. Histograms of the variables plotted are shown along the diagonal, with kernel density estimates of the distributions superimposed. Arguments \code{xlim}, \code{ylim}, \code{xhist}, and \code{yhist} are then ignored (argument \code{hh} can then be used to compress/stretch the height of the distributions shown along the diagonal).
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
### meta-analysis of all trials including ISIS-4 using a FE model
res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, method="FE")

\dontrun{
### fit FE model to all possible subsets (65535 models)
sav <- gosh(res)

### create GOSH plot
### red points for subsets that include and blue points
### for subsets that exclude study 16 (the ISIS-4 trial)
plot(sav, out=16, breaks=100, adjust=.3)}
}
\keyword{hplot}
