\name{forest.default}
\alias{forest.default}
\title{Forest Plots (Default Method)}
\description{Function to create forest plots for a given set of data.}
\usage{
\method{forest}{default}(x, vi, sei, ci.lb, ci.ub, annotate=TRUE, showweights=FALSE,
       xlim, alim, clim, ylim, at, steps=5,
       level=95, refline=0, digits=2, width,
       xlab, slab, ilab, ilab.xpos, ilab.pos,
       subset, transf, atransf, targs, rows,
       efac=1, pch=15, psize, col, lty,
       cex, cex.lab, cex.axis, \dots)
}
\arguments{
   \item{x}{vector of length \eqn{k} with the observed effect sizes or outcomes.}
   \item{vi}{vector of length \eqn{k} with the corresponding sampling variances.}
   \item{sei}{vector of length \eqn{k} with the corresponding standard errors.\cr(note: only one of the two, \code{vi} or \code{sei}, needs to be specified)}
   \item{ci.lb}{vector of length \eqn{k} with the corresponding lower confidence intervals bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{ci.ub}{vector of length \eqn{k} with the corresponding upper confidence intervals bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{annotate}{logical specifying whether annotations should be added to the plot (the default is \code{TRUE}).}
   \item{showweights}{logical specifying whether the annotations should also include inverse variance weights (the default is \code{FALSE}).}
   \item{xlim}{horizontal limits of the plot region. If unspecified, the function tries to set the horizontal plot limits to some sensible values.}
   \item{alim}{the actual x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{clim}{limits for the confidence intervals. If unspecified, no limits are used.}
   \item{ylim}{the y limits of the plot. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{steps}{the number of tick marks for the x-axis (the default is 5). Ignored when the user specifies the positions via the \code{at} argument.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is 95).}
   \item{refline}{value at which a vertical \sQuote{reference} line should be drawn (the default is 0). The line can be suppressed by setting this argument to \code{NA}.}
   \item{digits}{integer specifying the number of decimal places to which the tick mark labels of the x-axis and the annotations should be rounded (the default is 2). Can also be a vector of two integers, the first specifying the number of decimal places for the annotations, the second for the x-axis labels.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{slab}{optional vector with labels for the \eqn{k} studies. If unspecified, simple labels are created within the function. To suppress labels, set this argument to \code{NA}.}
   \item{ilab}{optional vector or matrix of character strings providing additional information about the studies.}
   \item{ilab.xpos}{numerical value(s) specifying the x-axis position(s) of the character vector(s) given via \code{ilab} (must be specified if \code{ilab} is specified).}
   \item{ilab.pos}{integer(s) (either 1, 2, 3, or 4) specifying the alignment of the character vector(s) given via \code{ilab} (2 means right, 4 mean left aligned). If unspecified, the default is to center the labels.}
   \item{subset}{optional vector indicating the subset of studies that should be used for the plot. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the observations to include. Note that this argument can also be used for reordering the studies.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the observed effect sizes or outcomes and corresponding confidence interval bounds (e.g., \code{transf=exp}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument specifying the name of a function that should be used to transform the x-axis labels and annotations (e.g., \code{transf=exp}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{rows}{optional vector specifying the rows (or more generally, the horizontal positions) for plotting the outcomes. If unspecified, the function sets this value automatically. Can also be a single value specifying the row (horizontal position) of the first outcome (the remaining outcomes are then plotted below this starting row).}
   \item{efac}{vertical expansion factor for arrows and confidence interval limits. The default value of 1 should usually work okay.}
   \item{pch}{plotting symbol to use for the observed effect sizes or outcomes. By default, a filled square is used. See \code{\link{points}} for other options. Can also be a vector of values.}
   \item{psize}{optional vector with point sizes for the observed effect sizes or outcomes. If unspecified, the point sizes are drawn proportional to the precision of the estimates.}
   \item{col}{optional character string specifying the name of a color to use for plotting ("black" is used by default if not specified). Can also be a vector of color names.}
   \item{lty}{optional character string specifying the line type for the confidence intervals (if unspecified, the function sets this to \code{"solid"} by default).}
   \item{cex}{optional character and symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.lab}{optional expansion factor for the x-axis titel. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.axis}{optional expansion factor for the x-axis labels. If unspecified, the function tries to set this to a sensible value.}
   \item{\dots}{other arguments.}
}
\details{
   The plot shows the individual observed effect sizes or outcomes with corresponding confidence intervals. To use the function, one should specify the observed effect sizes or outcomes (via the \code{x} argument) together with the corresponding sampling variances (via the \code{vi} argument) or with the corresponding standard errors (via the \code{sei} argument). Alternatively, one can specify the observed effect sizes or outcomes together with the corresponding confidence interval bounds (via the \code{ci.lb} and \code{ci.ub} arguments).

   With the \code{transf} argument, the observed effect sizes or outcomes and corresponding confidence interval bounds can be transformed with some suitable function. For example, when plotting log odds ratios, then one could use \code{transf=exp} to obtain a forest plot showing the odds ratios. Alternatively, one can use the \code{atransf} argument to transform the x-axis labels and annotations (e.g., \code{atransf=exp}). The examples below illustrate the use of these arguments.

   By default, the studies are ordered from top to bottom (i.e., the first study in the dataset will be placed in row \eqn{k}, the second study in row \eqn{k-1}, and so on, until the last study, which is placed in the first row). The studies can be reordered with the \code{subset} argument (by specifying a vector with indices with the desired order).

   Summary estimates can also be added to the plot with the \code{\link{addpoly}} function. See the documentation for that function for examples.
}
\note{
   The function tries to set some sensible values for the optional arguments, but it may be necessary to tweak these in certain circumstances.

   The function actually returns some information about the chosen defaults invisibly. Printing this information is useful as a starting point to make adjustments to the plot.

   If the number of studies is quite large, the labels, annotations, and symbols may become quite small and impossible to read. Stretching the plot window vertically may then provide a more readable figure (one should call the function again after adjusting the window size, so that the label/symbol sizes can be properly adjusted). Also, the \code{cex}, \code{cex.lab}, and \code{cex.axis} arguments are then useful to adjust the symbol and text sizes.

   If the horizontal plot and/or x-axis limits are set by the user, then the horizontal plot limits (\code{xlim}) must be at least as wide as the x-axis limits (\code{alim}). This restriction is enforced inside the function.

   If the outcome measure used for creating the plot has absolute limits (e.g., correlations are bounded between \code{-1} and \code{1}, proportions are bounded between \code{0} and \code{1}), one can use the \code{clim} argument to enforce those limits (confidence intervals cannot exceed those bounds then).

   The \code{lty} argument can also be a vector of two elements, the first for specifying the line type of the individual CIs (\code{"solid"} by default), the second for the line type of the horizontal line that is automatically added to the plot (\code{"solid"} by default; set to \code{"blank"} to remove it).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Lewis, S., & Clarke, M. (2001). Forest plots: Trying to see the wood and the trees. \emph{British Medical Journal}, \bold{322}, 1479--1480.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{forest}}, \code{\link{forest.rma}}, \code{\link{addpoly}}
}
\examples{
### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### default forest plot of the observed log relative risks
forest(dat$yi, dat$vi)

### forest plot of the observed relative risks
forest(dat$yi, dat$vi, slab=paste(dat$author, dat$year, sep=", "), transf=exp,
       alim=c(0,2), steps=5, xlim=c(-2.5,4), refline=1, cex=.9)

### forest plot of the observed relative risks
forest(dat$yi, dat$vi, slab=paste(dat$author, dat$year, sep=", "), atransf=exp,
       at=log(c(.05,.25,1,4,20)), xlim=c(-10,8), cex=.9)

### see also examples for the forest.rma function
}
\keyword{hplot}
