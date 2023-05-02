\name{plot.infl.rma.uni}
\alias{plot.infl.rma.uni}
\title{Plot Method for 'infl.rma.uni' Objects}
\description{
   Plot method for objects of class \code{"infl.rma.uni"}. \loadmathjax
}
\usage{
\method{plot}{infl.rma.uni}(x, plotinf=TRUE, plotdfbs=FALSE, dfbsnew=FALSE, logcov=TRUE,
     layout, slab.style=1, las=0, pch=21, bg, bg.infl, col.na, \dots)
}
\arguments{
   \item{x}{an object of class \code{"infl.rma.uni"} obtained with \code{\link[=influence.rma.uni]{influence}}.}
   \item{plotinf}{logical to specify whether the various case diagnostics should be plotted (the default is \code{TRUE}). Can also be a vector of up to 8 integers to specify which plots to draw. See \sQuote{Details} for the numbers corresponding to the various plots.}
   \item{plotdfbs}{logical to specify whether the DFBETAS values should be plotted (the default is \code{FALSE}). Can also be a vector of integers to specify for which coefficient(s) to plot the DFBETAS values.}
   \item{dfbsnew}{logical to specify whether a new device should be opened for plotting the DFBETAS values (the default is \code{FALSE}).}
   \item{logcov}{logical to specify whether the covariance ratios should be plotted on a log scale (the default is \code{TRUE}).}
   \item{layout}{optional vector of two numbers to specify the number of rows and columns for the layout of the figure.}
   \item{slab.style}{integer to indicate the style of the x-axis labels: 1 = study number, 2 = study label, 3 = abbreviated study label. Note that study labels, even when abbreviated, may be too long to fit in the margins.)}
   \item{las}{integer between 0 and 3 to specify the alignment of the axis labels (see \code{\link{par}}). The most useful alternative to 0 is 3, so that the x-axis labels are drawn vertical to the axis.}
   \item{pch}{plotting symbol to use. By default, an open circle is used. See \code{\link{points}} for other options.}
   \item{bg}{optional character string to specify the background color of open plotting symbols. If unspecified, gray is used by default.}
   \item{bg.infl}{optional character string to specify the background color when the point is considered influential. If unspecified, red is used by default.}
   \item{col.na}{optional character string to specify the color for lines connecting two points with \code{NA} values in between. If unspecified, a light shade of gray is used by default.}
   \item{\dots}{other arguments.}
}
\details{
   When \code{plotinf=TRUE}, the function plots the (1) externally standardized residuals, (2) DFFITS values, (3) Cook's distances, (4) covariance ratios, (5) leave-one-out \mjseqn{\tau^2} estimates, (6) leave-one-out (residual) heterogeneity test statistics, (7) hat values, and (8) weights. If \code{plotdfbs=TRUE}, the DFBETAS values are also plotted either after confirming the page change (if \code{dfbsnew=FALSE}) or on a separate device (if \code{dfbsnew=TRUE}).

   A case (which is typically synonymous with study) may be considered to be \sQuote{influential} if at least one of the following is true:
   \itemize{
   \item The absolute DFFITS value is larger than \mjeqn{3 \times \sqrt{p/(k-p)}}{3*\sqrt(p/(k-p))}, where \mjseqn{p} is the number of model coefficients and \mjseqn{k} the number of cases.
   \item The lower tail area of a chi-square distribution with \mjseqn{p} degrees of freedom cut off by the Cook's distance is larger than 50\%.
   \item The hat value is larger than \mjeqn{3 \times (p/k)}{3*(p/k)}.
   \item Any DFBETAS value is larger than \mjseqn{1}.
   } Cases which are considered influential with respect to any of these measures are indicated by the color specified for the \code{bg.infl} argument (the default is \code{"red"}).

   The cut-offs described above are indicated in the plot with horizontal reference lines. In addition, on the plot of the externally standardized residuals, horizontal reference lines are drawn at -1.96, 0, and 1.96. On the plot of the hat values, a horizontal reference line is drawn at \mjseqn{p/k}. Since the sum of the hat values is equal to \mjseqn{p}, the value \mjseqn{p/k} indicates equal hat values for all \mjseqn{k} cases. Finally, on the plot of weights, a horizontal reference line is drawn at \mjseqn{100/k}, corresponding to the value for equal weights (in \%) for all \mjseqn{k} cases. Note that all weights will automatically be equal to each other when using unweighted model fitting. Also, the hat values will be equal to the weights values (except for their scaling) in models without moderators.

   The chosen cut-offs are (somewhat) arbitrary. Substantively informed judgment should always be used when examining the influence of each case on the results.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}(2), 112--125. \verb{https://doi.org/10.1002/jrsm.11}
}
\seealso{
   \code{\link{influence.rma.uni}} for the function to compute the various model diagnostics.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### compute the diagnostics
inf <- influence(res)

### plot the values
plot(inf)

### select which plots to show
plot(inf, plotinf=1:4)
plot(inf, plotinf=1:4, layout=c(4,1))

### plot the DFBETAS values
plot(inf, plotinf=FALSE, plotdfbs=TRUE)
}
\keyword{hplot}
