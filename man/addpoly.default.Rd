\name{addpoly.default}
\alias{addpoly.default}
\title{Add Polygons to Forest Plots (Default Method)}
\description{
   Function to add one or more polygons to a forest plot.
}
\usage{
\method{addpoly}{default}(x, vi, sei, ci.lb, ci.ub, pi.lb, pi.ub,
        rows=-1, level, annotate, predstyle, predlim, digits, width, mlab,
        transf, atransf, targs, efac, col, border, lty, fonts, cex,
        constarea=FALSE, \dots)
}
\arguments{
   \item{x}{vector with the values at which the polygons should be drawn.}
   \item{vi}{vector with the corresponding variances.}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{ci.lb}{vector with the corresponding lower confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{ci.ub}{vector with the corresponding upper confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{pi.lb}{optional vector with the corresponding lower prediction interval bounds.}
   \item{pi.ub}{optional vector with the corresponding upper prediction interval bounds.}
   \item{rows}{vector to specify the rows (or more generally, the positions) for plotting the polygons (defaults is \code{-1}). Can also be a single value to specify the row of the first polygon (the remaining polygons are then plotted below this starting row). When \code{predstyle} is not \code{"line"}, can also be a vector of two numbers, the first for the position of the polygon, the second for the position of the prediction interval/distribution.}
   \item{level}{optional numeric value between 0 and 100 to specify the confidence interval level (see \link[=misc-options]{here} for details).}
   \item{annotate}{optional logical to specify whether annotations should be added to the plot for the polygons that are drawn.}
   \item{predstyle}{character string to specify the style of the prediction interval (either \code{"line"} (the default), \code{"polygon"}, \code{"bar"}, \code{"shade"}, or \code{"dist"}; the last three only when adding a single polygon). Can be abbreviated.}
   \item{predlim}{optional argument to specify the limits of the predictive distribution when \code{predstyle="dist"}.}
   \item{digits}{optional integer to specify the number of decimal places to which the annotations should be rounded.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{mlab}{optional character vector of the same length as \code{x} giving labels for the polygons that are drawn.}
   \item{transf}{optional argument to specify a function to transform the \code{x} values and confidence interval bounds (e.g., \code{transf=exp}; see also \link{transf}).}
   \item{atransf}{optional argument to specify a function to transform the annotations (e.g., \code{atransf=exp}; see also \link{transf}).}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{efac}{optional vertical expansion factor for the polygons.}
   \item{col}{optional character string to specify the color of the polygons.}
   \item{border}{optional character string to specify the border color of the polygons.}
   \item{lty}{optional argument to specify the line type for the prediction interval.}
   \item{fonts}{optional character string to specify the font for the labels and annotations.}
   \item{cex}{optional symbol expansion factor.}
   \item{constarea}{logical to specify whether the height of the polygons (when adding multiple) should be adjusted so that the area of the polygons is constant (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to add one or more polygons to an existing forest plot created with the \code{\link{forest}} function. For example, pooled estimates based on a model involving moderators can be added to the plot this way (see \sQuote{Examples}).

   To use the function, one should specify the values at which the polygons should be drawn (via the \code{x} argument) together with the corresponding variances (via the \code{vi} argument) or with the corresponding standard errors (via the \code{sei} argument). Alternatively, one can specify the values at which the polygons should be drawn together with the corresponding confidence interval bounds (via the \code{ci.lb} and \code{ci.ub} arguments). Optionally, one can also specify the bounds of the corresponding prediction interval bounds via the \code{pi.lb} and \code{pi.ub} arguments. If the latter are specified, then they are added by default as lines around the summary polygons. When adding a single polygon to the plot, one can also use the \code{predstyle} argument to change the way the prediction interval is visualized (see \code{\link{forest.rma}} for details).

   If unspecified, arguments \code{level}, \code{annotate}, \code{digits}, \code{width}, \code{transf}, \code{atransf}, \code{targs}, \code{efac} (only if the forest plot was created with \code{\link{forest.rma}}), \code{fonts}, \code{cex}, \code{annosym}, and \code{textpos} are automatically set equal to the same values that were used when creating the forest plot.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest}} for functions to draw forest plots to which polygons can be added.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude as a moderator
res <- rma(yi, vi, mods = ~ ablat, slab=paste(author, year, sep=", "), data=dat)

### forest plot of the observed risk ratios
forest(res, addfit=FALSE, atransf=exp, xlim=c(-9,5), ylim=c(-5,16), cex=0.9,
       order=ablat, ilab=ablat, ilab.lab="Lattitude", ilab.xpos=-4.5,
       header="Author(s) and Year")

### predicted average log risk ratios for 10, 30, and 50 degrees absolute latitude
x <- predict(res, newmods=c(10, 30, 50))

### add predicted average risk ratios to the forest plot
addpoly(x$pred, sei=x$se, rows=-2, mlab=c("- at 10 Degrees", "- at 30 Degrees", "- at 50 Degrees"))
abline(h=0)
text(-9, -1, "Model-Based Estimates:", pos=4, cex=0.9, font=2)
}
\keyword{aplot}
