\name{addpoly.predict.rma}
\alias{addpoly.predict.rma}
\title{Add Polygons to Forest Plots (Method for 'predict.rma' Objects)}
\description{
   Function to add one or more polygons to a forest plot based on an object of class \code{"predict.rma"}.
}
\usage{
\method{addpoly}{predict.rma}(x, rows=-2, annotate,
        addpred=FALSE, digits, width, mlab, transf, atransf, targs,
        efac, col, border, lty, fonts, cex, constarea=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"predict.rma"}.}
   \item{rows}{vector to specify the rows (or more generally, the horizontal positions) for plotting the polygons (defaults is \code{-2}). Can also be a single value to specify the row (horizontal position) of the first polygon (the remaining polygons are then plotted below this starting row).}
   \item{annotate}{optional logical to specify whether annotations should be added to the plot for the polygons that are drawn.}
   \item{addpred}{logical to specify whether the bounds of the prediction interval should be added to the plot (the default is \code{FALSE}).}
   \item{digits}{optional integer to specify the number of decimal places to which the annotations should be rounded.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{mlab}{optional character vector with the same length as \code{x} giving labels for the polygons that are drawn.}
   \item{transf}{optional argument to specify a function to transform the \code{x} values and confidence interval bounds (e.g., \code{transf=exp}; see also \link{transf}).}
   \item{atransf}{optional argument to specify a function to transform the annotations (e.g., \code{atransf=exp}; see also \link{transf}).}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{efac}{optional vertical expansion factor for the polygons.}
   \item{col}{optional character string to specify the color of the polygons.}
   \item{border}{optional character string to specify the border color of the polygons.}
   \item{lty}{optional character string to specify the line type for the prediction interval.}
   \item{fonts}{optional character string to specify the font for the labels and annotations.}
   \item{cex}{optional symbol expansion factor.}
   \item{constarea}{logical to specify whether the height of the polygons (when adding multiple) should be adjusted so that the area of the polygons is constant (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to add one or more polygons to an existing forest plot created with the \code{\link{forest}} function. For example, summary estimates based on a model involving moderators can be added to the plot this way (see \sQuote{Examples}).

   To use the function, one should specify the values at which the polygons should be drawn (via the \code{x} argument) together with the corresponding variances (via the \code{vi} argument) or with the corresponding standard errors (via the \code{sei} argument). Alternatively, one can specify the values at which the polygons should be drawn together with the corresponding confidence interval bounds (via the \code{ci.lb} and \code{ci.ub} arguments). Optionally, one can also specify the bounds of the corresponding prediction interval bounds via the \code{pi.lb} and \code{pi.ub} arguments.

   If unspecified, arguments \code{annotate}, \code{digits}, \code{width}, \code{transf}, \code{atransf}, \code{targs}, \code{efac} (only if the forest plot was created with \code{\link{forest.rma}}), \code{fonts}, \code{cex}, \code{annosym}, and \code{textpos} are automatically set equal to the same values that were used when creating the forest plot.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest}} for functions to draw forest plots to which polygons can be added.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
              data=dat.bcg, slab=paste(author, year, sep=", "))

### forest plot of the observed risk ratios
with(dat, forest(yi, vi, atransf=exp, xlim=c(-9,5), ylim=c(-4.5,15),
                 at=log(c(0.05, 0.25, 1, 4)), cex=0.9, order=alloc,
                 ilab=alloc, ilab.xpos=-4.5,
                 header="Author(s) and Year", top=2))

### fit mixed-effects model with allocation method as a moderator
res <- rma(yi, vi, mods = ~ 0 + alloc, data=dat)

### predicted log risk ratios for the different allocation methods
x <- predict(res, newmods=diag(3))

### add predicted average risk ratios to the forest plot
addpoly(x, efac=1.3, col="gray", addpred=TRUE,
        mlab=c("Alternate Allocation", "Random Allocation", "Systematic Allocation"))
abline(h=0)
text(-9, -1, "Model-Based Estimates:", pos=4, cex=0.9, font=2)
text(-4.5, res$k+2, "Allocation", cex=0.9, font=2)
}
\keyword{aplot}
