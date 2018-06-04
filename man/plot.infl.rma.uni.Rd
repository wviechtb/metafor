\name{plot.infl.rma.uni}
\alias{plot.infl.rma.uni}
\title{Plot Method for 'infl.rma.uni' Objects}
\description{
   Plot method for objects of class \code{"infl.rma.uni"}.
}
\usage{
\method{plot}{infl.rma.uni}(x, plotinf=TRUE, plotdfbs=FALSE, dfbsnew=FALSE, logcov=TRUE,
     layout, slab.style=1, las=0, pch=21, bg="black",
     bg.infl="red", col.na="lightgray", \dots)
}
\arguments{
   \item{x}{an object of class \code{"infl.rma.uni"}.}
   \item{plotinf}{logical indicating whether the various case diagnostics should be plotted (the default is \code{TRUE}). Can also be a vector of up to 8 integers indicating which plots to draw (see \sQuote{Details} for the numbers corresponding to the various plots).}
   \item{plotdfbs}{logical indicating whether the DFBETAS values should be plotted (the default is \code{FALSE}). Can also be a vector of integers to indicate for which coefficient(s) to plot the DFBETAS values.}
   \item{dfbsnew}{logical indicating whether a new device should be opened for plotting the DFBETAS values (the default is \code{FALSE}).}
   \item{logcov}{logical indicating whether the covariance ratios should be plotted on a log scale (the default is \code{TRUE}).}
   \item{layout}{optional vector of two numbers, indicating the number of rows and columns for the layout of the figure.}
   \item{slab.style}{integer to indicate the style of the x-axis labels: 1 = study number, 2 = study label, 3 = abbreviated study label. Note that study labels, even when abbreviated, may be too long to fit in the margin.)}
   \item{las}{integer between 0 and 3 to specify the alignment of the axis labels (see \code{\link{par}}). The most useful alternative to 0 is 3, so that the x-axis labels are drawn vertical to the axis.}
   \item{pch}{plotting symbol to use. By default, a filled circle is used. See \code{\link{points}} for other options.}
   \item{bg}{color to use for filling the plotting symbol (the default is \code{"black"}).}
   \item{bg.infl}{color to use for filling the plotting symbol when the point is considered influential (the default is \code{"red"}).}
   \item{col.na}{color to use for lines connecting two points with \code{NA} values in between (the default is \code{"lightgray"}).}
   \item{\dots}{other arguments.}
}
\details{
   When \code{plotinf=TRUE}, the function plots the (1) externally standardized residuals, (2) DFFITS values, (3) Cook's distances, (4) covariance ratios, (5) leave-one-out \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} estimates, (6) leave-one-out (residual) heterogeneity test statistics, (7) hat values, and (8) weights. If \code{plotdfbs=TRUE}, the DFBETAS values are also plotted either after confirming the page change (if \code{dfbsnew=FALSE}) or on a separate device (if \code{dfbsnew=TRUE}).

   A study may be considered to be \sQuote{influential} if at least one of the following is true:
   \itemize{
   \item The absolute DFFITS value is larger than \eqn{3\sqrt{p/(k-p)}}{3\sqrt(p/(k-p))}, where \eqn{p} is the number of model coefficients and \eqn{k} the number of studies.
   \item The lower tail area of a chi-square distribution with \eqn{p} degrees of freedom cut off by the Cook's distance is larger than 50\%.
   \item The hat value is larger than \eqn{3(p/k)}.
   \item Any DFBETAS value is larger than \eqn{1}.
   } Studies which are considered influential with respect to any of these measures are indicated by the color specified for the \code{bg.infl} argument (the default is \code{"red"}).

   The cut-offs described above are indicated in the plot with horizontal reference lines. In addition, on the plot of the externally standardized residuals, horizontal reference lines are drawn at -1.96, 0, and 1.96. On the plot of the hat values, a horizontal reference line is drawn at \eqn{p/k}. Since the sum of the hat values is equal to \eqn{p}, the value \eqn{p/k} indicates equal hat values for all \eqn{k} studies. Finally, on the plot of weights, a horizontal reference line is drawn at \eqn{100/k}, corresponding to the value for equal weights (in \%) for all \eqn{k} studies. Note that all weights will automatically be equal to each other when using unweighted model fitting. Also, the hat values will be equal to the weights values (except for their scaling) in models without moderators.

   The chosen cut-offs are (somewhat) arbitrary. Substantively informed judgment should always be used when examining the influence of each study on the results.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}, 112--125.
}
\seealso{
   \code{\link{influence.rma.uni}}
}
\examples{
### meta-analysis of the log risk ratios using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year, data=dat.bcg)
plot(influence(res))
plot(influence(res), plotinf=1:4)
plot(influence(res), plotinf=1:4, layout=c(4,1))
plot(influence(res), plotdfbs=TRUE, dfbsnew=TRUE)
}
\keyword{hplot}
