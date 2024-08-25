\name{addpoly.rma}
\alias{addpoly.rma}
\title{Add Polygons to Forest Plots (Method for 'rma' Objects)}
\description{
   Function to add a polygon to a forest plot showing the pooled estimate with corresponding confidence interval based on an object of class \code{"rma"}.
}
\usage{
\method{addpoly}{rma}(x, row=-2, level=x$level, annotate,
        addpred=FALSE, digits, width, mlab, transf, atransf, targs,
        efac, col, border, lty, fonts, cex, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{row}{numeric value to specify the row (or more generally, the position) for plotting the polygon (the default is \code{-2}).}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (see \link[=misc-options]{here} for details). The default is to take the value from the object.}
   \item{annotate}{optional logical to specify whether annotations for the pooled estimate should be added to the plot.}
   \item{addpred}{logical to specify whether the bounds of the prediction interval should be added to the plot (the default is \code{FALSE}).}
   \item{digits}{optional integer to specify the number of decimal places to which the annotations should be rounded.}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{mlab}{optional character string giving a label for the pooled estimate. If unspecified, the function sets a default label.}
   \item{transf}{optional argument to specify a function to transform the pooled estimate and confidence interval bounds (e.g., \code{transf=exp}; see also \link{transf}).}
   \item{atransf}{optional argument to specify a function to transform the annotations (e.g., \code{atransf=exp}; see also \link{transf}).}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{efac}{optional vertical expansion factor for the polygon.}
   \item{col}{optional character string to specify the color of the polygon.}
   \item{border}{optional character string to specify the border color of the polygon.}
   \item{lty}{optional argument to specify the line type for the prediction interval.}
   \item{fonts}{optional character string to specify the font for the label and annotations.}
   \item{cex}{optional symbol expansion factor.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to add a four-sided polygon, sometimes called a summary \sQuote{diamond}, to an existing forest plot created with the \code{\link{forest}} function. The polygon shows the pooled estimate (with its confidence interval bounds) based on an equal- or a random-effects model. Using this function, pooled estimates based on different types of models can be shown in the same plot. Also, pooled estimates based on a subgrouping of the studies can be added to the plot this way. See \sQuote{Examples}.

   If unspecified, arguments \code{annotate}, \code{digits}, \code{width}, \code{transf}, \code{atransf}, \code{targs}, \code{efac} (only if the forest plot was created with \code{\link{forest.rma}}), \code{fonts}, \code{cex}, \code{annosym}, and \code{textpos} are automatically set equal to the same values that were used when creating the forest plot.
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
### meta-analysis of the log risk ratios using the Mantel-Haenszel method
res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg,
              slab=paste(author, year, sep=", "))

### forest plot of the observed risk ratios with the pooled estimate
forest(res, atransf=exp, xlim=c(-8,6), ylim=c(-3,16))

### meta-analysis of the log risk ratios using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### add the pooled estimate from the random-effects model to the forest plot
addpoly(res)

### forest plot with subgrouping of studies and summaries per subgroup
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg,
              slab=paste(author, year, sep=", "))
res <- rma(yi, vi, data=dat)
tmp <- forest(res, xlim=c(-16, 4.6), at=log(c(0.05, 0.25, 1, 4)), atransf=exp,
              ilab=cbind(tpos, tneg, cpos, cneg), ilab.lab=c("TB+","TB-","TB+","TB-"),
              ilab.xpos=c(-9.5,-8,-6,-4.5), cex=0.75, ylim=c(-2, 27), order=alloc,
              rows=c(3:4,9:15,20:23), mlab="RE Model for All Studies",
              header="Author(s) and Year")
op <- par(cex=tmp$cex)
text(c(-8.75,-5.25), tmp$ylim[2]-0.2, c("Vaccinated", "Control"), font=2)
text(-16, c(24,16,5), c("Systematic Allocation", "Random Allocation",
                        "Alternate Allocation"), font=4, pos=4)
par(op)
res <- rma(yi, vi, data=dat, subset=(alloc=="systematic"))
addpoly(res, row=18.5, mlab="RE Model for Subgroup")
res <- rma(yi, vi, data=dat, subset=(alloc=="random"))
addpoly(res, row=7.5, mlab="RE Model for Subgroup")
res <- rma(yi, vi, data=dat, subset=(alloc=="alternate"))
addpoly(res, row=1.5, mlab="RE Model for Subgroup")
}
\keyword{aplot}
