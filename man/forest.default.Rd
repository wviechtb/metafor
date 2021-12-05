\name{forest.default}
\alias{forest.default}
\title{Forest Plots (Default Method)}
\description{
   Function to create forest plots for a given set of data. \loadmathjax
}
\usage{
\method{forest}{default}(x, vi, sei, ci.lb, ci.ub,
       annotate=TRUE, showweights=FALSE, header=FALSE,
       xlim, alim, olim, ylim, top=3, at, steps=5,
       level=95, refline=0, digits=2L, width,
       xlab, slab, ilab, ilab.xpos, ilab.pos,
       order, subset, transf, atransf, targs, rows,
       efac=1, pch=15, psize, plim=c(0.5,1.5), col,
       lty, fonts, cex, cex.lab, cex.axis, annosym, \dots)
}
\arguments{
   \item{x}{vector of length \mjseqn{k} with the observed effect sizes or outcomes.}
   \item{vi}{vector of length \mjseqn{k} with the corresponding sampling variances.}
   \item{sei}{vector of length \mjseqn{k} with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{ci.lb}{vector of length \mjseqn{k} with the corresponding lower confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{ci.ub}{vector of length \mjseqn{k} with the corresponding upper confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{annotate}{logical to specify whether annotations should be added to the plot (the default is \code{TRUE}).}
   \item{showweights}{logical to specify whether the annotations should also include the inverse variance weights (the default is \code{FALSE}).}
   \item{header}{logical to specify whether column headings should be added to the plot (the default is \code{FALSE}). Can also be a character vector to specify the left and right headings  (or only the left one).}
   \item{xlim}{horizontal limits of the plot region. If unspecified, the function tries to set the horizontal plot limits to some sensible values.}
   \item{alim}{the x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{olim}{optional argument to specify observation/outcome limits. If unspecified, no limits are used.}
   \item{ylim}{the y-axis limits of the plot. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{top}{the amount of space to leave empty at the top of the plot (e.g., for adding headers) (the default is 3 rows).}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{steps}{the number of tick marks for the x-axis (the default is 5). Ignored when the positions are specified via the \code{at} argument.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95).}
   \item{refline}{numeric value to specify the location of the vertical \sQuote{reference} line (the default is 0). The line can be suppressed by setting this argument to \code{NA}.}
   \item{digits}{integer to specify the number of decimal places to which the tick mark labels of the x-axis and the annotations should be rounded (the default is \code{2L}). Can also be a vector of two integers, the first to specify the number of decimal places for the annotations, the second for the x-axis labels. When specifying an integer (e.g., \code{2L}), trailing zeros after the decimal mark are dropped for the x-axis labels. When specifying a numeric value (e.g., \code{2}), trailing zeros are retained.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{slab}{optional vector with labels for the \mjseqn{k} studies. If unspecified, simple labels are created within the function. To suppress labels, set this argument to \code{NA}.}
   \item{ilab}{optional vector, matrix, or data frame providing additional information about the studies that should be added to the plot.}
   \item{ilab.xpos}{numeric vector to specify the x-axis position(s) of the variable(s) given via \code{ilab} (must be specified if \code{ilab} is specified).}
   \item{ilab.pos}{integer(s) (either 1, 2, 3, or 4) to specify the alignment of the vector(s) given via \code{ilab} (2 means right, 4 mean left aligned). If unspecified, the default is to center the labels.}
   \item{order}{optional character string to specify how the studies should be ordered. Can also be a variable based on which the studies will be ordered. See \sQuote{Details}.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be included in the plot.}
   \item{transf}{optional argument to specify a function that should be used to transform the observed outcomes and corresponding confidence interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument to specify a function that should be used to transform the x-axis labels and annotations (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{rows}{optional vector to specify the rows (or more generally, the horizontal positions) for plotting the outcomes. Can also be a single value to specify the row (horizontal position) of the first outcome (the remaining outcomes are then plotted below this starting row). If unspecified, the function sets this value automatically.}
   \item{efac}{vertical expansion factor for confidence interval limits and arrows. The default value of 1 should usually work okay. Can also be a vector of two numbers, the first for CI limits, the second for arrows.}
   \item{pch}{plotting symbol to use for the observed outcomes. By default, a filled square is used. See \code{\link{points}} for other options. Can also be a vector of values.}
   \item{psize}{optional numeric value to specify the point sizes for the observed outcomes. If unspecified, the point sizes are a function of the precision of the estimates. Can also be a vector of values.}
   \item{plim}{numeric vector of length 2 to scale the point sizes (ignored when \code{psize} is specified). See \sQuote{Details}.}
   \item{col}{optional character string to specify the name of a color to use for plotting the observed outcomes (\code{"black"} is used by default if not specified). Can also be a vector of color names.}
   \item{lty}{optional character string to specify the line type for the confidence intervals. If unspecified, the function sets this to \code{"solid"} by default.}
   \item{fonts}{optional character string to specify the font to use for the study labels, annotations, and the extra information (if specified via \code{ilab}). If unspecified, the default font is used.}
   \item{cex}{optional character and symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.lab}{optional expansion factor for the x-axis title. If unspecified, the function tries to set this to a sensible value.}
   \item{cex.axis}{optional expansion factor for the x-axis labels. If unspecified, the function tries to set this to a sensible value.}
   \item{annosym}{optional vector of length 3 to change the left bracket, separation, and right bracket symbols for the annotations.}
   \item{\dots}{other arguments.}
}
\details{
   The plot shows the observed effect sizes or outcomes with corresponding confidence intervals. To use the function, one should specify the observed outcomes (via the \code{x} argument) together with the corresponding sampling variances (via the \code{vi} argument) or with the corresponding standard errors (via the \code{sei} argument). Alternatively, one can specify the observed outcomes together with the corresponding confidence interval bounds (via the \code{ci.lb} and \code{ci.ub} arguments).

   With the \code{transf} argument, the observed outcomes and corresponding confidence interval bounds can be transformed with some suitable function. For example, when plotting log odds ratios, then one could use \code{transf=exp} to obtain a forest plot showing the odds ratios. Alternatively, one can use the \code{atransf} argument to transform the x-axis labels and annotations (e.g., \code{atransf=exp}). See also \link{transf} for some other useful transformation functions in the context of a meta-analysis. The examples below illustrate the use of these arguments.

   By default, the studies are ordered from top to bottom (i.e., the first study in the dataset will be placed in row \mjseqn{k}, the second study in row \mjseqn{k-1}, and so on, until the last study, which is placed in the first row). The studies can be reordered with the \code{order} argument:
   \itemize{
   \item \code{order="obs"}: the studies are ordered by the observed outcomes,
   \item \code{order="prec"}: the studies are ordered by their sampling variances.
   } Alternatively, it is also possible to set \code{order} equal to a variable based on which the studies will be ordered (see \sQuote{Examples}).

   Additional columns with information about the studies can be added to the plot via the \code{ilab} argument. This can either be a single variable or an entire matrix / data frame (with as many rows as there are studies in the forest plot). The \code{ilab.xpos} argument must then also be specified to indicate the x-axis position of the variables specified via \code{ilab}.

   Summary estimates can be added to the plot with the \code{\link{addpoly}} function. See the documentation for that function for examples.

   By default (i.e., when \code{psize} is not specified), the point sizes are a function of the precision (i.e., inverse standard errors) of the outcomes. This way, more precise estimates are visually more prominent in the plot. By making the point sizes a function of the inverse standard errors of the estimates, their areas are proportional to the inverse sampling variances, which corresponds to the weights they would receive in an equal-effects model. However, the point sizes are rescaled so that the smallest point size is \code{plim[1]} and the largest point size is \code{plim[2]}. As a result, their relative sizes (i.e., areas) no longer exactly correspond to their relative weights in such a model. If exactly relative point sizes are desired, one can set \code{plim[2]} to \code{NA}, in which case the points are rescaled so that the smallest point size corresponds to \code{plim[1]} and all other points are scaled accordingly. As a result, the largest point may be very large. Alternatively, one can set \code{plim[1]} to \code{NA}, in which case the points are rescaled so that the largest point size corresponds to \code{plim[2]} and all other points are scaled accordingly. As a result, the smallest point may be very small and essentially indistinguishable from the confidence interval line. To avoid the latter, one can also set \code{plim[3]}, which enforces a minimal point size.
}
\note{
   The function tries to set some sensible values for the optional arguments, but it may be necessary to adjust these in certain circumstances.

   The function actually returns some information about the chosen defaults invisibly. Printing this information is useful as a starting point to make adjustments to the plot.

   If the number of studies is quite large, the labels, annotations, and symbols may become quite small and impossible to read. Stretching the plot window vertically may then provide a more readable figure (one should call the function again after adjusting the window size, so that the label/symbol sizes can be properly adjusted). Also, the \code{cex}, \code{cex.lab}, and \code{cex.axis} arguments are then useful to adjust the symbol and text sizes.

   If the horizontal plot and/or x-axis limits are set manually, then the horizontal plot limits (\code{xlim}) must be at least as wide as the x-axis limits (\code{alim}). This restriction is enforced inside the function.

   If the outcome measure used for creating the plot is bounded (e.g., correlations are bounded between -1 and +1, proportions are bounded between 0 and 1), one can use the \code{olim} argument to enforce those limits (the observed outcomes and confidence intervals cannot exceed those bounds then).

   The \code{lty} argument can also be a vector of two elements, the first for specifying the line type of the individual CIs (\code{"solid"} by default), the second for the line type of the horizontal line that is automatically added to the plot (\code{"solid"} by default; set to \code{"blank"} to remove it).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Lewis, S., & Clarke, M. (2001). Forest plots: Trying to see the wood and the trees. \emph{British Medical Journal}, \bold{322}(7300), 1479--1480. \verb{https://doi.org/10.1136/bmj.322.7300.1479}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest}} for an overview of the various \code{forest} functions and especially \code{\link{forest.rma}} for the function draw forest plots including a summary polygon.

   \code{\link{addpoly}} for a function to add polygons to forest plots.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### default forest plot of the observed log risk ratios
forest(dat$yi, dat$vi)

### forest plot of the observed risk ratios (transform outcomes)
forest(dat$yi, dat$vi, slab=paste(dat$author, dat$year, sep=", "), transf=exp,
       alim=c(0,2), steps=5, xlim=c(-2.5,4), refline=1, cex=.9, header=TRUE)

### forest plot of the observed risk ratios (transformed x-axis)
forest(dat$yi, dat$vi, slab=paste(dat$author, dat$year, sep=", "), atransf=exp,
       at=log(c(.05,.25,1,4,20)), xlim=c(-10,8), cex=.9, header=TRUE)

### forest plot of the observed risk ratios with studies ordered by the RRs
forest(dat$yi, dat$vi, slab=paste(dat$author, dat$year, sep=", "), atransf=exp,
       at=log(c(.05,.25,1,4,20)), xlim=c(-10,8), cex=.9, header=TRUE, order="obs")

### forest plot of the observed risk ratios with studies ordered by absolute latitude
forest(dat$yi, dat$vi, slab=paste(dat$author, dat$year, sep=", "), atransf=exp,
       at=log(c(.05,.25,1,4,20)), xlim=c(-10,8), cex=.9, header=TRUE, order=dat$ablat)

### see also examples for the forest.rma function
}
\keyword{hplot}
