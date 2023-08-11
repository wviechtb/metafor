\name{radial}
\alias{radial}
\alias{galbraith}
\alias{radial.rma}
\title{Radial (Galbraith) Plots for 'rma' Objects}
\description{
   Function to create radial (also called Galbraith) plots for objects of class \code{"rma"}. \loadmathjax
}
\usage{
radial(x, \dots)
galbraith(x, \dots)

\method{radial}{rma}(x, center=FALSE, xlim, zlim, xlab, zlab,
       atz, aty, steps=7, level=x$level, digits=2,
       transf, targs, pch=21, col, bg, back, arc.res=100,
       cex, cex.lab, cex.axis, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{center}{logical to indicate whether the plot should be centered horizontally at the model estimate (the default is \code{FALSE}).}
   \item{xlim}{x-axis limits. If unspecified, the function sets the x-axis limits to some sensible values.}
   \item{zlim}{z-axis limits. If unspecified, the function sets the z-axis limits to some sensible values (note that the z-axis limits are the actual vertical limit of the plotting region).}
   \item{xlab}{title for the x-axis. If unspecified, the function sets an appropriate axis title.}
   \item{zlab}{title for the z-axis. If unspecified, the function sets an appropriate axis title.}
   \item{atz}{position for the z-axis tick marks and labels. If unspecified, these values are set by the function.}
   \item{aty}{position for the y-axis tick marks and labels. If unspecified, these values are set by the function.}
   \item{steps}{the number of tick marks for the y-axis (the default is 7). Ignored when argument \code{aty} is used.}
   \item{level}{numeric value between 0 and 100 to specify the level of the z-axis error region. The default is to take the value from the object.}
   \item{digits}{integer to specify the number of decimal places to which the tick mark labels of the y-axis should be rounded (the default is 2).}
   \item{transf}{argument to specify a function to transform the y-axis labels (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{transf}.}
   \item{pch}{plotting symbol. By default, an open circle is used. See \code{\link{points}} for other options.}
   \item{col}{character string to specify the (border) color of the points.}
   \item{bg}{character string to specify the background color of open plot symbols.}
   \item{back}{character string to specify the background color of the z-axis error region. If unspecified, a shade of gray is used. Set to \code{NA} to suppress shading of the region.}
   \item{arc.res}{integer to specify the number of line segments (i.e., the resolution) when drawing the y-axis and confidence interval arcs (the default is 100).}
   \item{cex}{symbol expansion factor.}
   \item{cex.lab}{character expansion factor for axis labels.}
   \item{cex.axis}{character expansion factor for axis annotations.}
   \item{\dots}{other arguments.}
}
\details{
   For an equal-effects model, the plot shows the inverse of the standard errors on the horizontal axis (i.e., \mjeqn{1/\sqrt{v_i}}{1/\sqrt(v_i)}, where \mjseqn{v_i} is the sampling variance of the observed effect size or outcome) against the observed effect sizes or outcomes standardized by their corresponding standard errors on the vertical axis (i.e., \mjeqn{y_i/\sqrt{v_i}}{y_i/\sqrt(v_i)}). Since the vertical axis corresponds to standardized values, it is referred to as the z-axis within this function. On the right hand side of the plot, an arc is drawn (referred to as the y-axis within this function) corresponding to the observed effect sizes or outcomes. A line projected from (0,0) through a particular point within the plot onto this arc indicates the value of the observed effect size or outcome for that point.

   For a random-effects model, the function uses \mjeqn{1/\sqrt{v_i + \tau^2}}{1/\sqrt(v_i + \tau^2)} for the horizontal axis, where \mjseqn{\tau^2} is the amount of heterogeneity as estimated based on the model. For the z-axis, \mjeqn{y_i/\sqrt{v_i + \tau^2}}{y_i/\sqrt(v_i + \tau^2)} is used to compute standardized values of the observed effect sizes or outcomes.

   The second (inner/smaller) arc that is drawn on the right hand side indicates the model estimate (in the middle of the arc) and the corresponding confidence interval (at the ends of the arc).

   The shaded region in the plot is the z-axis error region. For \code{level=95} (or if this was the \code{level} value when the model was fitted), this corresponds to z-axis values equal to \mjeqn{\pm 1.96}{+-1.96}. Under the assumptions of the equal/random-effects models, approximately 95\% of the points should fall within this region.

   When \code{center=TRUE}, the values on the y-axis are centered around the model estimate. As a result, the plot is centered horizontally at the model estimate.

   If the z-axis label on the left is too close to the actual z-axis and/or the arc on the right is clipped, then this can be solved by increasing the margins on the right and/or left (see \code{\link{par}} and in particular the \code{mar} argument).

   Note that radial plots cannot be drawn for models that contain moderators.
}
\value{
   A data frame with components:
   \item{x}{the x-axis coordinates of the points that were plotted.}
   \item{y}{the y-axis coordinates of the points that were plotted.}
   \item{ids}{the study id numbers.}
   \item{slab}{the study labels.}
   Note that the data frame is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Galbraith, R. F. (1988). Graphical display of estimates having differing standard errors. \emph{Technometrics}, \bold{30}(3), 271--281. \verb{https://doi.org/10.1080/00401706.1988.10488400}

   Galbraith, R. F. (1988). A note on graphical presentation of estimated odds ratios from several clinical trials. \emph{Statistics in Medicine}, \bold{7}(8), 889--894. \verb{https://doi.org/10.1002/sim.4780070807}

   Galbraith, R. F (1994). Some applications of radial plots. \emph{Journal of the American Statistical Association}, \bold{89}(428), 1232--1242. \verb{https://doi.org/10.1080/01621459.1994.10476864}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which radial plots can be drawn.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat

### fit equal-effects model
res <- rma(yi, vi, data=dat, method="EE")

### draw radial plot
radial(res)

### the line from (0,0) with a slope equal to the log risk ratio from the 4th study points
### to the corresponding effect size value on the arc (i.e., -1.44)
abline(a=0, b=dat$yi[4], lty="dotted")
dat$yi[4]

### meta-analysis of the log risk ratios using a random-effects model
res <- rma(yi, vi, data=dat)

### draw radial plot
radial(res)

### center the values around the model estimate
radial(res, center=TRUE)

### show risk ratio values on the y-axis arc
radial(res, transf=exp)
}
\keyword{hplot}
