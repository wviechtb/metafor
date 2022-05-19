\name{plot.cumul.rma}
\alias{plot.cumul.rma}
\title{Plot Method for 'cumul.rma' Objects}
\description{
   Plot method for objects of class \code{"cumul.rma"}. \loadmathjax
}
\usage{
\method{plot}{cumul.rma}(x, yaxis, xlim, ylim, xlab, ylab,
     at, transf, atransf, targs, digits, cols=c("gray80","gray10"),
     grid=TRUE, pch=19, cex=1, lwd=2, \dots)
}
\arguments{
   \item{x}{an object of class \code{"cumul.rma"} obtained with \code{\link{cumul}}.}
   \item{yaxis}{either \code{"tau2"}, \code{"I2"}, or \code{"H2"} to indicate what values should be placed on the y-axis. See \sQuote{Details}.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{transf}{optional argument to specify a function to transform the summary estimates (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument to specify a function to transform the x-axis labels (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{digits}{optional integer to specify the number of decimal places to which the tick mark labels of the x- and y-axis should be rounded. Can also be a vector of two integers, the first to specify the number of decimal places for the x-axis, the second for the y-axis labels (e.g., \code{digits=c(2,3)}). If unspecified, the function tries to set the argument to some sensible values.}
   \item{cols}{vector with two or more colors to use for visualizing the order of the cumulative results.}
   \item{grid}{logical to specify whether a grid should be added to the plot. Can also be a color name.}
   \item{pch}{plotting symbol to use. By default, a filled circle is used. See \code{\link{points}} for other options.}
   \item{cex}{symbol expansion factor.}
   \item{lwd}{line width.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to visualize the results from a cumulative meta-analysis as obtained with the \code{\link{cumul}} function.

   The plot shows the model estimate (i.e., the estimated overall/average outcome) on the x-axis and some measure of heterogeneity on the y-axis in the cumulative order of the results in the \code{"cumul.rma"} object. By default, \mjseqn{\tau^2} is shown on the y-axis for a random-effects model and \mjseqn{I^2} otherwise, but one can also use argument \code{yaxis} to specify the measure of heterogeneity to place on the y-axis.

   The color gradient of the points/lines indicates the order of the cumulative results (by default, light gray at the beginning, dark gray at the end). A different set of colors can be chosen via the \code{cols} argument. See \sQuote{Examples}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{cumul.rma.uni}} for the function to conduct a cumulative meta-analysis.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

### cumulative meta-analysis (in the order of publication year)
sav <- cumul(res, transf=exp, order=year)

### plot of model estimate and tau^2 over time
plot(sav)

### illustrate some other plot options
plot(sav, yaxis="I2", ylim=c(0,100), atransf=exp, at=log(seq(1.3, 1.6, by=.1)),
     lwd=5, cex=1.5, cols=c("green","blue","red"))
}
\keyword{hplot}
