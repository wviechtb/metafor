\name{addpoly.default}
\alias{addpoly.default}
\title{Add Polygons to Forest Plots (Default Method)}
\description{
   Function to add one or more polygons to a forest plot.
}
\usage{
\method{addpoly}{default}(x, vi, sei, ci.lb, ci.ub, pi.lb, pi.ub,
        rows=-1, level=95, annotate=TRUE, digits=2, width, mlab,
        transf, atransf, targs, efac=1, col, border, fonts, cex, \dots)
}
\arguments{
   \item{x}{vector with the values at which the polygons should be drawn.}
   \item{vi}{vector with the corresponding variances.}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{ci.lb}{vector with the corresponding lower confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{ci.ub}{vector with the corresponding upper confidence interval bounds. Not needed if \code{vi} or \code{sei} is specified. See \sQuote{Details}.}
   \item{pi.lb}{optional vector with the corresponding lower prediction interval bounds.}
   \item{pi.ub}{optional vector with the corresponding upper prediction interval bounds.}
   \item{rows}{vector specifying the rows (or more generally, the horizontal positions) for plotting the polygons (defaults is \code{-1}). Can also be a single value specifying the row (horizontal position) of the first polygon (the remaining polygons are then plotted below this starting row).}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is 95).}
   \item{annotate}{logical specifying whether annotations should be added to the plot for the polygons that are drawn (the default is \code{TRUE}).}
   \item{digits}{integer specifying the number of decimal places to which the annotations should be rounded (the default is 2).}
   \item{width}{optional integer to manually adjust the width of the columns for the annotations.}
   \item{mlab}{optional character vector with the same length as \code{x} giving labels for the polygons that are drawn.}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the \code{x} values and confidence interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{atransf}{optional argument specifying the name of a function that should be used to transform the annotations (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf} or \code{atransf}.}
   \item{efac}{vertical expansion factor for the polygons. The default value of 1 should usually work okay.}
   \item{col}{optional character string specifying the name of a color to use for the polygons. If unspecified, the function sets a default color.}
   \item{border}{optional character string specifying the name of a color to use for the border of the polygons. If unspecified, the function sets a default color.}
   \item{fonts}{optional character string specifying the font to use for the labels and annotations. If unspecified, the default font is used.}
   \item{cex}{optional symbol expansion factor. If unspecified, the function tries to set this to a sensible value.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to add one or more polygons to an existing forest plot created with the \code{\link{forest}} function. For example, summary estimates based on a model involving moderators can be added to the plot this way (see \sQuote{Examples}).

   To use the function, one should specify the values at which the polygons should be drawn (via the \code{x} argument) together with the corresponding variances (via the \code{vi} argument) or with the corresponding standard errors (via the \code{sei} argument). Alternatively, one can specify the values at which the polygons should be drawn together with the corresponding confidence interval bounds (via the \code{ci.lb} and \code{ci.ub} arguments). Optionally, one can also specify the bounds of the corresponding prediction interval bounds via the \code{pi.lb} and \code{pi.ub} arguments.

   The arguments \code{transf}, \code{atransf}, \code{efac}, and \code{cex} should always be set equal to the same values used to create the forest plot.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest.rma}}, \code{\link{forest.default}}
}
\examples{
### meta-analysis of the log risk ratios using a mixed-
### effects model with absolute latitude as a moderator
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, mods = ~ ablat,
           slab=paste(author, year, sep=", "), data=dat.bcg)

### forest plot of the observed risk ratios
forest(res, addfit=FALSE, atransf=exp, xlim=c(-8,5), ylim=c(-4.5,16), cex=.8,
       order=order(dat.bcg$ablat), ilab=dat.bcg$ablat, ilab.xpos=-4,
       header="Author(s) and Year")

### predicted average log risk ratios for 10, 30, and 50 degrees absolute latitude
x <- predict(res, newmods=c(10, 30, 50))

### add predicted average risk ratios to forest plot
addpoly(x$pred, sei=x$se, atransf=exp, col="white", rows=-2,
        mlab=c("- at 10 Degrees", "- at 30 Degrees", "- at 50 Degrees"), cex=.8)
abline(h=0)
text(-8, -1, "Model-Based Estimates of RR:", pos=4, cex=.8)
text(-4, 15, "Latitude", cex=.8, font=2)
}
\keyword{aplot}
