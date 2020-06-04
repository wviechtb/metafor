\name{forest.cumul.rma}
\alias{forest.cumul.rma}
\title{Forest Plots (Method for 'cumul.rma' Objects)}
\description{
   Function to create forest plots for objects of class \code{"cumul.rma"}.
}
\usage{
\method{forest}{cumul.rma}(x, annotate=TRUE, header=FALSE,
       xlim, alim, clim, ylim, top=3, at, steps=5,
       level=x$level, refline=0, digits=2L, width,
       xlab, ilab, ilab.xpos, ilab.pos,
       transf, atransf, targs, rows,
       efac=1, pch=15, psize=1, col, lty, fonts,
       cex, cex.lab, cex.axis, annosym, \dots)
}
\arguments{
   \item{x}{an object of class \code{"cumul.rma"}.}
   \item{annotate}{logical specifying whether annotations should be added to the plot (the default is \code{TRUE}).}
   \item{header}{logical specifying whether column headings should be added to the plot (the default is \code{FALSE}). Can also be a character vector to specify the left and right headings.}
   \item{xlim}{horizontal limits of the plot region. If unspecified, the function tries to set the horizontal plot limits to some sensible values.}
   \item{alim}{the actual x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{clim}{limits for the confidence intervals. If unspecified, no limits are used.}
   \item{ylim}{the y limits of the plot. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{top}{the amount of space to leave empty at the top of the plot (e.g., for adding headers) (the default is 3 rows).}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{steps}{the number of tick marks for the x-axis (the default is 5). Ignored when the user specifies the positions via the \code{at} argument.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is to take the value from the object).}
   \item{refline}{value at which a vertical \sQuote{reference} line should be drawn (the default is 0). The line can be suppressed by setting this argument to \code{NA}.}
   \item{digits}{integer specifying the number of decimal places to which the tick mark labels of the x-axis and the annotations should be rounded (the default is \code{2L}). Can also be a vector of two integers, the first specifying the number of decimal places for the annotations, the second for the x-axis labels. When specifying an integer (e.g., \code{2L}), trailing zeros after the decimal mark are dropped for the x-axis labels. When specifying a numerical value (e.g., \code{2}), trailing zeros are retained.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ilab}{optional vector, matrix, or data frame providing additional information about the studies that should be added to the plot.}
   \item{ilab.xpos}{vector of numerical value(s) specifying the x-axis position(s) of the variable(s) given via \code{ilab} (must be specified if \code{ilab} is specified).}
   \item{ilab.pos}{integer(s) (either 1, 2, 3, or 4) specifying the alignment of the vector(s) given via \code{ilab} (2 means right, 4 mean left aligned). If unspecified, the default is to center the labels.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the observed effect sizes, summary estimates, fitted values, and confidence interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument specifying the name of a function that should be used to transform the x-axis labels and annotations (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{rows}{optional vector specifying the rows (or more generally, the horizontal positions) for plotting the outcomes. If unspecified, the function sets this value automatically. Can also be a single value specifying the row (horizontal position) of the first outcome (the remaining outcomes are then plotted below this starting row).}
   \item{efac}{vertical expansion factor for confidence interval limits and arrows. The default value of 1 should usually work okay. Can also be a vector of two numbers, the first for CI limits, the second for arrows.}
   \item{pch}{plotting symbol to use for the observed effect sizes or outcomes. By default, a filled square is used. See \code{\link{points}} for other options. Can also be a vector of values.}
   \item{psize}{optional vector with point sizes for the observed effect sizes or outcomes. Default is 1.}
   \item{col}{optional character string specifying the name of a color to use for plotting ("black" is used by default if not specified). Can also be a vector of color names.}
   \item{lty}{optional character string specifying the line type for the confidence intervals (if unspecified, the function sets this to \code{"solid"} by default).}
   \item{fonts}{optional character string specifying the font to use for the study labels, annotations, and the extra information (if specified via \code{ilab}). If unspecified, the default font is used.}
   \item{cex}{optional character and symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.lab}{optional expansion factor for the x-axis title. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.axis}{optional expansion factor for the x-axis labels. If unspecified, the function tries to set this to a sensible value.}
   \item{annosym}{optional vector of length 3 to change the left bracket, separation, and right bracket symbols for the annotations.}
   \item{\dots}{other arguments.}
}
\details{
   The plot shows the estimated (average) outcome with corresponding confidence interval as one study at a time is added to the analysis.
}
\note{
   The function tries to set some sensible values for the optional arguments, but it may be necessary to tweak these in certain circumstances.

   The function actually returns some information about the chosen defaults invisibly. Printing this information is useful as a starting point to make adjustments to the plot.

   If the number of studies is quite large, the labels, annotations, and symbols may become quite small and impossible to read. Stretching the plot window vertically may then provide a more readable figure (one should call the function again after adjusting the window size, so that the label/symbol sizes can be properly adjusted). Also, the \code{cex}, \code{cex.lab}, and \code{cex.axis} arguments are then useful to adjust the symbol and text sizes.

   If the horizontal plot and/or x-axis limits are set by the user, then the horizontal plot limits (\code{xlim}) must be at least as wide as the x-axis limits (\code{alim}). This restriction is enforced inside the function.

   If the outcome measure used for creating the plot is bounded (e.g., correlations are bounded between \code{-1} and \code{1}, proportions are bounded between \code{0} and \code{1}), one can use the \code{clim} argument to enforce those limits (confidence intervals cannot exceed those bounds then).

   The \code{lty} argument can also be a vector of two elements, the first for specifying the line type of the individual CIs (\code{"solid"} by default), the second for the line type of the horizontal line that is automatically added to the plot (\code{"solid"} by default; set to \code{"blank"} to remove it).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Chalmers, T. C., & Lau, J. (1993). Meta-analytic stimulus for changes in clinical trials. \emph{Statistical Methods in Medical Research}, \bold{2}, 161--172.

   Lau, J., Schmid, C. H., & Chalmers, T. C. (1995). Cumulative meta-analysis of clinical trials builds evidence for exemplary medical care. \emph{Journal of Clinical Epidemiology}, \bold{48}, 45--57.

   Lewis, S., & Clarke, M. (2001). Forest plots: Trying to see the wood and the trees. \emph{British Medical Journal}, \bold{322}, 1479--1480.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\seealso{
   \code{\link{forest}}, \code{\link{cumul}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat, slab=paste(author, year, sep=", "))

x <- cumul(res, order=order(dat$year))
forest(x, cex=.8)
forest(x, alim=c(-2,1), cex=.8)

### meta-analysis of the (log) risk ratios using the Mantel-Haenszel method
res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg,
              slab=paste(author, year, sep=", "))
x <- cumul(res, order=order(dat$year))
forest(x, alim=c(-2,1), cex=.8)
}
\keyword{hplot}
