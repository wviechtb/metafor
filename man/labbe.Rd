\name{labbe}
\alias{labbe}
\alias{labbe.rma}
\title{L'Abbe Plots for 'rma' Objects}
\description{
   Function to create \enc{L'Abbé}{L'Abbe} plots for objects of class \code{"rma"}. \loadmathjax
}
\usage{
labbe(x, \dots)

\method{labbe}{rma}(x, xlim, ylim, lim, xlab, ylab, ci=FALSE,
      add=x$add, to=x$to, transf, targs,
      pch=21, psize, plim=c(0.5,3.5),
      col, bg, grid=FALSE, lty, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{xlim}{x-axis limits. If unspecified, the function sets the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function sets the y-axis limits to some sensible values.}
   \item{lim}{axis limits. If specified, this is used for both \code{xlim} and \code{ylim}.}
   \item{xlab}{title for the x-axis. If unspecified, the function sets an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function sets an appropriate axis title.}
   \item{ci}{logical to indicate whether the confidence interval region should be shown in the plot (the default is \code{FALSE}). Can also be a string to specify the color for shading the confidence region.}
   \item{add}{See the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{See the documentation of the \code{\link{escalc}} function for more details.}
   \item{transf}{optional argument to specify a function to transform the outcomes (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{pch}{plotting symbol to use for the outcomes. By default, an open circle is used. Can also be a vector of values. See \code{\link{points}} for other options.}
   \item{psize}{optional numeric vector to specify the point sizes for the outcomes. If unspecified, the point sizes are a function of the precision of the outcomes. Can also be a vector of values.}
   \item{plim}{numeric vector of length 2 to scale the point sizes (ignored when \code{psize} is specified). See \sQuote{Details}.}
   \item{col}{optional character string to specify the (border) color of the points. Can also be a vector.}
   \item{bg}{optional character string to specify the background color of open plot symbols. Can also be a vector. Set to \code{NA} to make the plotting symbols transparent.}
   \item{grid}{logical to specify whether a grid should be added to the plot (the default is \code{FALSE}). Can also be a color name.}
   \item{lty}{optional character vector to specify the line type for the diagonal reference line of no effect, the line that indicates the estimated effect based on the fitted model, and the corresponding confidence interval region. If unspecified, the function sets this to \code{c("solid","dashed","dotted")} by default (use \code{"blank"} to suppress a line).}
   \item{\dots}{other arguments.}
}
\details{
   The model specified via \code{x} must be a model without moderators (i.e., either an equal- or a random-effects model) fitted with either the \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, or \code{\link{rma.glmm}} functions. Moreover, the model must have been fitted with \code{measure} set equal to \code{"RD"} (for risk differences), \code{"RR"} (for risk ratios), \code{"OR"} (for odds ratios), \code{"AS"} (for arcsine square root transformed risk differences), \code{"IRR"} (for incidence rate ratios), \code{"IRD"} (for incidence rate differences), or \code{"IRSD"} (for square root transformed incidence rate differences).

   The function calculates the arm-level outcomes for the two groups (e.g., treatment and control) and plots them against each other. In particular, the function plots the raw proportions of the two groups against each other when analyzing risk differences, the log of the proportions when analyzing (log) risk ratios, the log odds when analyzing (log) odds ratios, the arcsine square root transformed proportions when analyzing arcsine square root transformed risk differences, the raw incidence rates when analyzing incidence rate differences, the log of the incidence rates when analyzing (log) incidence rate ratios, and the square root transformed incidence rates when analyzing square root transformed incidence rate differences. The \code{transf} argument can be used to transform these values (e.g., \code{transf=exp} to transform the log of the proportions back to raw proportions; see also \link{transf}).

   As described under the documentation for the \code{\link{escalc}} function, zero cells can lead to problems when calculating particular outcomes. Adding a small constant to the cells of the \mjeqn{2 \times 2}{2x2} tables is a common solution to this problem. By default, the functions adopts the same method for handling zero cells as was used when fitting the model.

   By default (i.e., when \code{psize} is not specified), the point sizes are a function of the precision (i.e., inverse standard errors) of the outcomes. This way, more precise estimates are visually more prominent in the plot. By making the point sizes a function of the inverse standard errors of the estimates, their areas are proportional to the inverse sampling variances, which corresponds to the weights they would receive in an equal-effects model. However, the point sizes are rescaled so that the smallest point size is \code{plim[1]} and the largest point size is \code{plim[2]}. As a result, their relative sizes (i.e., areas) no longer exactly correspond to their relative weights in such a model. If exactly relative point sizes are desired, one can set \code{plim[2]} to \code{NA}, in which case the points are rescaled so that the smallest point size corresponds to \code{plim[1]} and all other points are scaled accordingly. As a result, the largest point may be very large. Alternatively, one can set \code{plim[1]} to \code{NA}, in which case the points are rescaled so that the largest point size corresponds to \code{plim[2]} and all other points are scaled accordingly. As a result, the smallest point may be very small. To avoid the latter, one can also set \code{plim[3]}, which enforces a minimal point size.

   The solid line corresponds to identical outcomes in the two groups (i.e., the absence of a difference between the two groups). The dashed line indicates the estimated effect based on the fitted model. If \code{ci=TRUE}, then the shaded region indicates the corresponding confidence interval.
}
\value{
   A data frame with components:
   \item{x}{the x-axis coordinates of the points that were plotted.}
   \item{y}{the y-axis coordinates of the points that were plotted.}
   \item{cex}{the point sizes.}
   \item{pch}{the plotting symbols.}
   \item{col}{the point colors.}
   \item{bg}{the background colors.}
   \item{ids}{the study id numbers.}
   \item{slab}{the study labels.}
   Note that the data frame is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   \enc{L'Abbé}{L'Abbe}, K. A., Detsky, A. S., & O'Rourke, K. (1987). Meta-analysis in clinical research. \emph{Annals of Internal Medicine}, \bold{107}(2), 224--233. \verb{https://doi.org/10.7326/0003-4819-107-2-224}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.glmm}} for functions to fit models for which \enc{L'Abbé}{L'Abbe} plots can be drawn.
}
\examples{
### meta-analysis of the log risk ratios using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### default plot
labbe(res)

### plot with risk values on the x- and y-axis and some further customization
labbe(res, transf=exp, grid=TRUE, ci=TRUE, bty="l", lty=c(1,2,0))
legend("topleft", inset=0.02, bg="white", pch=c(NA,NA,22), pt.cex=2.5, pt.lwd=0,
       pt.bg="lightgray", lty=c(1,2,0), legend=c("Reference Line of No Effect",
       "Line for the Estimated Effect", "95\% Confidence Interval"))
}
\keyword{hplot}
