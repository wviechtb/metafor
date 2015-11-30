\name{plot.cumul.rma}
\alias{plot.cumul.rma}
\title{Plot Method for 'cumul.rma' Objects}
\description{
   Plot method for objects of class \code{"cumul.rma"}.
}
\usage{
\method{plot}{cumul.rma}(x, yaxis="tau2", xlim, ylim, xlab, ylab,
     at, transf, atransf, targs, digits, cols=c("gray90","gray10"),
     addgrid=TRUE, pch=19, cex=1, lwd=2, \dots)
}
\arguments{
   \item{x}{an object of class \code{"cumul.rma"}.}
   \item{yaxis}{either \code{"tau2"}, \code{"I2"}, or \code{"H2"} to indicate what values should be placed on the y-axis. See \sQuote{Details}.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the summary estimates (e.g., \code{transf=exp}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument specifying the name of a function that should be used to transform the x-axis labels (e.g., \code{transf=exp}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{digits}{integer specifying the number of decimal places to which the tick mark labels of the x- and y-axis should be rounded. Can also be a vector of two integers, the first specifying the number of decimal places for the x-axis, the second for the y-axis labels (e.g., c(2,3)). If unspecified, the function tries to set the argument to some sensible values.}
   \item{cols}{vector with two or more colors to be used for indicating the order of the cumulative results.}
   \item{addgrid}{logical specifying whether a grid should be added to the plot (default is \code{TRUE}).}
   \item{pch}{plotting symbol to use. By default, a solid circle is used. See \code{\link{points}} for other options.}
   \item{cex}{symbol expansion factor.}
   \item{lwd}{line width.}
   \item{\dots}{other arguments.}
}
\details{
   The function is meant to be used for cumulative random-effects models results as obtained with the \code{\link{cumul.rma.uni}} function.

   By default, the plot shows the model estimate on the x-axis and the corresponding estimate of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} on the y-axis in the cumulative order of the results in the \code{"cumul.rma"} object. Alternatively, one can also place the value of \if{latex}{\eqn{I^2}}\if{html}{\eqn{I²}}\if{text}{\eqn{I^2}} (\code{yaxis="I2"}) or \if{latex}{\eqn{H^2}}\if{html}{\eqn{H²}}\if{text}{\eqn{H^2}} (\code{yaxis="H2"}) on the y-axis.

   The color gradient of the points/lines indicates the order of the cumulative results (by default, light gray at the beginning, dark gray at the end). A different set of colors to be used can be chosen via the \code{cols} argument. See \sQuote{Examples}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{cumul.rma.uni}}
}
\examples{
### load BCG vaccine data
data(dat.bcg)

### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

### cumulative meta-analysis (in the order of publication year)
sav <- cumul(res, transf=exp, order=order(dat$year))

### plot of model estimate and tau^2 over time
plot(sav)

### illustrate some other plot options
plot(sav, yaxis="I2", ylim=c(0,100), atransf=exp, at=log(seq(1.3, 1.6, by=.1)),
     lwd=5, cex=0, cols=c("green","blue","red"))
}
\keyword{hplot}
