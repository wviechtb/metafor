\name{qqnorm.rma}
\alias{qqnorm}
\alias{qqnorm.rma}
\alias{qqnorm.rma.uni}
\alias{qqnorm.rma.mh}
\alias{qqnorm.rma.peto}
\alias{qqnorm.rma.glmm}
\alias{qqnorm.rma.mv}
\title{Normal QQ Plots for 'rma' Objects}
\description{
   Function to create normal QQ plots for objects of class \code{"rma.uni"}, \code{"rma.mh"}, and \code{"rma.peto"}. \loadmathjax
}
\usage{
\method{qqnorm}{rma.uni}(y, type="rstandard", pch=21, col, bg,
       envelope=TRUE, level=y$level, bonferroni=FALSE, reps=1000, smooth=TRUE, bass=0,
       label=FALSE, offset=0.3, pos=13, lty, \dots)
\method{qqnorm}{rma.mh}(y, type="rstandard", pch=21, col, bg, label=FALSE, offset=0.3, pos=13, \dots)
\method{qqnorm}{rma.peto}(y, type="rstandard", pch=21, col, bg, label=FALSE, offset=0.3, pos=13, \dots)
\method{qqnorm}{rma.glmm}(y, \dots)
\method{qqnorm}{rma.mv}(y, \dots)
}
\arguments{
   \item{y}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}. The method is not yet implemented for objects of class \code{"rma.glmm"} or \code{"rma.mv"}.}
   \item{type}{character string (either \code{"rstandard"} (default) or \code{"rstudent"}) to specify whether standardized residuals or studentized deleted residuals should be used in creating the plot. See \sQuote{Details}.}
   \item{pch}{plotting symbol to use for the observed outcomes. By default, an open circle is used. See \code{\link{points}} for other options.}
   \item{col}{optional character string to specify the (border) color of the points. If unspecified, black is used by default.}
   \item{bg}{optional character string to specify the background color of open plot symbols. If unspecified, gray is used by default.}
   \item{envelope}{logical to specify whether a pseudo confidence envelope should be simulated and added to the plot (the default is \code{TRUE})). Only for objects of class \code{"rma.uni"}. See \sQuote{Details}.}
   \item{level}{numeric value between 0 and 100 to specify the level of the pseudo confidence envelope (see \link[=misc-options]{here} for details). The default is to take the value from the object.}
   \item{bonferroni}{logical to specify whether the bounds of the envelope should be Bonferroni corrected.}
   \item{reps}{numeric value to specify the number of iterations for simulating the pseudo confidence envelope (the default is 1000).}
   \item{smooth}{logical to specify whether the results from the simulation should be smoothed (the default is \code{TRUE}).}
   \item{bass}{numeric value that controls the degree of smoothing (the default is 0).}
   \item{label}{argument to control the labeling of the points (the default is \code{FALSE}). See \sQuote{Details}.}
   \item{offset}{argument to control the distance between the points and the corresponding labels.}
   \item{pos}{argument to control the position of the labels.}
   \item{lty}{optional character string to specify the line type for the diagonal line and the pseudo confidence envelope. If unspecified, the function sets this to \code{c("solid", "dotted")} by default.}
   \item{\dots}{other arguments.}
}
\details{
   The plot shows the theoretical quantiles of a normal distribution on the horizontal axis against the observed quantiles for either the standardized residuals (\code{type="rstandard"}, the default) or the externally standardized residuals (\code{type="rstudent"}) on the vertical axis (see \code{\link[=residuals.rma]{residuals}} for details on the definition of these residual types).

   For reference, a line is added to the plot with a slope of 1, going through the (0,0) point.

   For objects of class \code{"rma.uni"}, it is also possible to add a pseudo confidence envelope to the plot. The envelope is created based on the quantiles of sets of pseudo residuals simulated from the given model (for details, see Cook & Weisberg, 1982). The number of sets simulated can be controlled with the \code{reps} argument. When \code{smooth=TRUE}, the simulated bounds are smoothed with Friedman's SuperSmoother (see \code{\link{supsmu}}). The \code{bass} argument can be set to a number between 0 and 10, with higher numbers indicating increasing smoothness. If \code{bonferroni=TRUE}, the envelope bounds are Bonferroni corrected, so that the envelope can be regarded as a confidence region for all \mjseqn{k} residuals simultaneously. The default however is \code{bonferroni=FALSE}, which makes the plot more sensitive to deviations from normality.

   With the \code{label} argument, one can control whether points in the plot will be labeled (e.g., to identify outliers). If \code{label="all"} (or \code{label=TRUE}), all points in the plot will be labeled. If \code{label="out"}, points falling outside of the confidence envelope will be labeled (only available for objects of class \code{"rma.uni"}). Finally, one can also set this argument to a numeric value (between 1 and \mjseqn{k}), indicating how many of the most extreme points should be labeled (for example, with \code{label=1} only the most extreme point is labeled, while with \code{label=3}, the most extreme, and the second and third most extreme points is labeled). With the \code{offset} argument, one can adjust the distance between the labels and the corresponding points. The \code{pos} argument is the position specifier for the labels (\code{1}, \code{2}, \code{3}, and \code{4}, respectively indicate positions below, to the left of, above, and to the right of the points; \code{13} places the labels below the points for points that fall below the reference line and above otherwise; \code{24} places the labels to the left of the points for points that fall above the reference line and to the right otherwise).
}
\value{
   A list with components:
   \item{x}{the x-axis coordinates of the points that were plotted.}
   \item{y}{the y-axis coordinates of the points that were plotted.}
   Note that the list is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Cook, R. D., & Weisberg, S. (1982). \emph{Residuals and influence in regression}. London: Chapman and Hall.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W. (2021). Model checking in meta-analysis. In C. H. Schmid, T. Stijnen, & I. R. White (Eds.), \emph{Handbook of meta-analysis} (pp. 219--254). Boca Raton, FL: CRC Press. \verb{https://doi.org/10.1201/9781315119403}

   Wang, M. C., & Bushman, B. J. (1998). Using the normal quantile plot to explore meta-analytic data sets. \emph{Psychological Methods}, \bold{3}(1), 46--54. \verb{https://doi.org/10.1037/1082-989X.3.1.46}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, and \code{\link{rma.peto}} for functions to fit models for which normal QQ plots can be drawn.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model
res <- rma(yi, vi, data=dat)

### draw QQ plot
qqnorm(res)

### fit mixed-effects model with absolute latitude as moderator
res <- rma(yi, vi, mods = ~ ablat, data=dat)

### draw QQ plot
qqnorm(res)
}
\keyword{hplot}
