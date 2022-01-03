\name{funnel}
\alias{funnel}
\alias{funnel.rma}
\alias{funnel.default}
\title{Funnel Plots}
\description{
   Function to create funnel plots. \loadmathjax
}
\usage{
funnel(x, \dots)

\method{funnel}{rma}(x, yaxis="sei",
       xlim, ylim, xlab, ylab,
       steps=5, at, atransf, targs, digits, level=x$level,
       addtau2=FALSE, type="rstandard",
       back="lightgray", shade="white", hlines="white",
       refline, lty=3, pch=19, pch.fill=21, col, bg,
       label=FALSE, offset=0.4, legend=FALSE, \dots)

\method{funnel}{default}(x, vi, sei, ni, subset, yaxis="sei",
       xlim, ylim, xlab, ylab,
       steps=5, at, atransf, targs, digits, level=95,
       back="lightgray", shade="white", hlines="white",
       refline=0, lty=3, pch=19, col, bg,
       label=FALSE, offset=0.4, legend=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"} or a vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances (needed if \code{x} is a vector with the observed effect sizes or outcomes).}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{ni}{vector with the corresponding sample sizes. Only relevant when passing a vector via \code{x}.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be included in the plot. Only relevant when passing a vector via \code{x}.}
   \item{yaxis}{either \code{"sei"}, \code{"vi"}, \code{"seinv"}, \code{"vinv"}, \code{"ni"}, \code{"ninv"}, \code{"sqrtni"}, \code{"sqrtninv"}, \code{"lni"}, or \code{"wi"} to indicate what values should be placed on the y-axis. See \sQuote{Details}.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{steps}{the number of tick marks for the y-axis (the default is 5).}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{atransf}{optional argument to specify a function that should be used to transform the x-axis labels (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{atransf}.}
   \item{digits}{optional integer to specify the number of decimal places to which the tick mark labels of the x- and y-axis should be rounded. Can also be a vector of two integers, the first to specify the number of decimal places for the x-axis, the second for the y-axis labels (e.g., \code{digits=c(2,3)}). If unspecified, the function tries to set the argument to some sensible values.}
   \item{level}{numeric value between 0 and 100 to specify the level of the pseudo confidence interval region (for \code{"rma"} objects, the default is to take the value from the object). May also be a vector of values to obtain multiple regions. See \sQuote{Examples}.}
   \item{addtau2}{logical to indicate whether the amount of heterogeneity should be accounted for when drawing the pseudo confidence interval region (the default is \code{FALSE}). Ignored when \code{x} is a meta-regression model and residuals are plotted. See \sQuote{Details}.}
   \item{type}{either \code{"rstandard"} (default) or \code{"rstudent"} to specify whether the usual or deleted residuals should be used in creating the funnel plot when \code{x} is a meta-regression model. See \sQuote{Details}.}
   \item{back}{color to use for the background of the plotting region (default is \code{"lightgray"}).}
   \item{shade}{color to use for shading the pseudo confidence interval region (default is \code{"white"}). When \code{level} is a vector of values, different shading colors can be specified for each region.}
   \item{hlines}{color of the horizontal reference lines (default is \code{"white"}).}
   \item{refline}{numeric value to specify the location of the vertical \sQuote{reference} line and where the pseudo confidence interval should be centered. If unspecified, the reference line is drawn at the equal- or random-effects model estimate and at zero for meta-regression models (in which case the residuals are plotted) or when directly plotting observed outcomes.}
   \item{lty}{line type for the pseudo confidence interval region and the reference line. The default is to draw dotted lines (see \code{\link{par}} for other options). Can also be a vector to specify the two line types separately.}
   \item{pch}{plotting symbol to use for the observed outcomes. By default, a filled circle is used. Can also be a vector of values. See \code{\link{points}} for other options.}
   \item{pch.fill}{plotting symbol to use for the outcomes filled in by the trim and fill method. By default, a circle is used. Only relevant when plotting an object created by the \code{\link{trimfill}} function.}
   \item{col}{optional character string to specify the color to use for the points (\code{"black"} is used by default if not specified). Can also be a vector.}
   \item{bg}{optional character string to specify the background color for open plot symbols (\code{"white"} is used by default if not specified). Can also be a vector.}
   \item{label}{argument to control the labeling of the points (the default is \code{FALSE}). See \sQuote{Details}.}
   \item{offset}{argument to control the distance between the points and the corresponding labels.}
   \item{legend}{logical to indicate whether a legend should be added to the plot (the default is \code{FALSE}). Can also be a keyword to indicate the position of the legend (see \code{\link{legend}}).}
   \item{\dots}{other arguments.}
}
\details{
   For equal- and random-effects models (i.e., models not involving moderators), the plot shows the observed effect sizes or outcomes on the x-axis against the corresponding standard errors (i.e., the square root of the sampling variances) on the y-axis. A vertical line indicates the estimate based on the model. A pseudo confidence interval region is drawn around this value with bounds equal to \mjeqn{\pm 1.96 \mbox{SE}}{±1.96*SE}, where \mjeqn{\mbox{SE}}{SE} is the standard error value from the y-axis (assuming \code{level=95}). If \code{addtau2=TRUE} (only for models of class \code{"rma.uni"}), then the bounds of the pseudo confidence interval region are equal to \mjeqn{\pm 1.96 \sqrt{\mbox{SE}^2 + \hat{\tau}^2}}{±1.96*\sqrt(SE^2 + \tau^2)}, where \mjeqn{\hat{\tau}^2}{\tau^2} is the amount of heterogeneity as estimated by the model.

   For (mixed-effects) meta-regression models (i.e., models involving moderators), the plot shows the residuals on the x-axis against their corresponding standard errors. Either the usual or deleted residuals can be used for that purpose (set via the \code{type} argument). See \code{\link[=residuals.rma]{residuals}} for more details on the different types of residuals.

   With the \code{atransf} argument, the labels on the x-axis can be transformed with some suitable function. For example, when plotting log odds ratios, one could use \code{transf=exp} to obtain a funnel plot with the values on the x-axis corresponding to the odds ratios. See also \link{transf} for some other useful transformation functions in the context of a meta-analysis.

   Instead of placing the standard errors on the y-axis, several other options are available by setting the \code{yaxis} argument to:
   \itemize{
   \item \code{yaxis="vi"} for the sampling variances,
   \item \code{yaxis="seinv"} for the inverse of the standard errors,
   \item \code{yaxis="vinv"} for the inverse of the sampling variances,
   \item \code{yaxis="ni"} for the sample sizes,
   \item \code{yaxis="ninv"} for the inverse of the sample sizes,
   \item \code{yaxis="sqrtni"} for the square root of the sample sizes,
   \item \code{yaxis="sqrtninv"} for the inverse square root of the sample sizes,
   \item \code{yaxis="lni"} for the log of the sample sizes,
   \item \code{yaxis="wi"} for the weights.
   } However, only when \code{yaxis="sei"} (the default) will the pseudo confidence region have the expected (upside-down) funnel shape with straight lines. Also, when placing (a function of) the sample sizes on the y-axis or the weights, then the pseudo confidence region cannot be drawn. See Sterne and Egger (2001) for more details on the choice of the y-axis.

   If the object passed to the function comes from the \code{\link{trimfill}} function, the outcomes that are filled in by the trim and fill method are also added to the funnel plot. The symbol to use for plotting the filled in values can be specified via the \code{pch.fill} argument.

   One can also directly pass a vector with the observed effect sizes or outcomes (via \code{x}) and the corresponding sampling variances (via \code{vi}), standard errors (via \code{sei}), and/or sample sizes (via \code{ni}) to the function. By default, the vertical reference line is then drawn at zero.

   The arguments \code{back}, \code{shade}, and \code{hlines} can be set to \code{NULL} to suppress the shading and the horizontal reference line.

   With the \code{label} argument, one can control whether points in the plot will be labeled. If \code{label="all"} (or \code{label=TRUE}), all points in the plot will be labeled. If \code{label="out"}, points falling outside of the pseudo confidence region will be labeled. Finally, one can also set this argument to a numeric value (between 1 and \mjseqn{k}) to specify how many of the most extreme points should be labeled (e.g., with \code{label=1} only the most extreme point would be labeled, while with \code{label=3}, the most extreme, and the second and third most extreme points would be labeled). With the \code{offset} argument, one can adjust the distance between the labels and the corresponding points.
}
\note{
   Placing (a function of) the sample sizes on the y-axis (i.e., using \code{yaxis="ni"}, \code{yaxis="ninv"}, \code{yaxis="sqrtni"}, \code{yaxis="sqrtninv"}, or \code{yaxis="lni"}) is only possible when information about the sample sizes is actually stored within the object passed to the \code{funnel} function. That should automatically be the case when the observed effect sizes or outcomes were computed with the \code{\link{escalc}} function or when the observed effect sizes or outcomes were computed within the model fitting function. On the other hand, this will not be the case when \code{\link{rma.uni}} was used together with the \code{yi} and \code{vi} arguments and the \code{yi} and \code{vi} values were \emph{not} computed with \code{\link{escalc}}. In that case, it is still possible to pass information about the sample sizes to the \code{\link{rma.uni}} function (e.g., use \code{rma.uni(yi, vi, ni=ni, data=dat)}, where data frame \code{dat} includes a variable called \code{ni} with the sample sizes).

   When using unweighted estimation, using \code{yaxis="wi"} will place all points on a horizontal line. When directly passing a vector with the observed effect sizes or outcomes to the function, \code{yaxis="wi"} is equivalent to \code{yaxis="vinv"}, except that the weights are expressed in percent.

   When specifying vectors for \code{pch}, \code{col}, and/or \code{bg}, the variables specified are assumed to be of the same length as the data passed to the funnel function or the model fitting function (when using \code{funnel} on a model object). Any subsetting and removal of studies with missing values is automatically applied to the variables specified via these arguments.
}
\value{
   A data frame with components:
   \item{x}{the x-axis coordinates of the points that were plotted.}
   \item{y}{the y-axis coordinates of the points that were plotted.}
   \item{slab}{the study labels.}
   Note that the data frame is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Light, R. J., & Pillemer, D. B. (1984). \emph{Summing up: The science of reviewing research}. Cambridge, MA: Harvard University Press.

   Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2008). Contour-enhanced meta-analysis funnel plots help distinguish publication bias from other causes of asymmetry. \emph{Journal of Clinical Epidemiology}, \bold{61}(10), 991--996. \verb{https://doi.org/10.1016/j.jclinepi.2007.11.010}

   Sterne, J. A. C., & Egger, M. (2001). Funnel plots for detecting bias in meta-analysis: Guidelines on choice of axis. \emph{Journal of Clinical Epidemiology}, \bold{54}(10), 1046--1055. \verb{https://doi.org/10.1016/s0895-4356(01)00377-8}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which funnel plots can be drawn.

   \code{\link{trimfill}} for the trim and fill method and \code{\link{regtest}} for the regression test.
}
\examples{
### copy BCG vaccine data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### fit random-effects model
res <- rma(yi, vi, data=dat, slab=paste(dat$author, dat$year, sep=", "))

### draw a standard funnel plot
funnel(res)

### show risk ratio values on x-axis (log scale)
funnel(res, atransf=exp)

### label points outside of the pseudo confidence interval region
funnel(res, atransf=exp, label="out")

### passing log risk ratios and sampling variances directly to the function
### note: same plot, except that reference line is centered at zero
funnel(dat$yi, dat$vi)

### can accomplish the same thing by setting refline=0
funnel(res, refline=0)

### adjust the position of the x-axis labels, number of digits, and y-axis limits
funnel(res, atransf=exp, at=log(c(.125, .25, .5, 1, 2)), digits=3L, ylim=c(0,.8))

### contour-enhanced funnel plot centered at 0 (see Peters et al., 2008)
funnel(res, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=TRUE)

### same, but show risk ratio values on the x-axis and some further adjustments
funnel(res, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), digits=3L, ylim=c(0,.8),
       refline=0, legend=TRUE, atransf=exp, at=log(c(.125, .25, .5, 1, 2, 4, 8)))

### illustrate the use of vectors for 'pch' and 'col'
res <- rma(yi, vi, data=dat, subset=2:10)
funnel(res, pch=ifelse(dat$yi > -1, 19, 21), col=ifelse(sqrt(dat$vi) > .3, "red", "blue"))

### can add a second funnel via (undocumented) argument refline2
funnel(res, atransf=exp, at=log(c(.125, .25, .5, 1, 2, 4)), digits=3L, ylim=c(0,.8), refline2=0)

### mixed-effects model with absolute latitude in the model
res <- rma(yi, vi, mods = ~ ablat, data=dat)

### funnel plot of the residuals
funnel(res)

### simulate a large meta-analytic dataset (correlations with rho = 0.2)
### with no heterogeneity or publication bias; then try out different
### versions of the funnel plot

gencor <- function(rhoi, ni) {
   x1 <- rnorm(ni, mean=0, sd=1)
   x2 <- rnorm(ni, mean=0, sd=1)
   x3 <- rhoi*x1 + sqrt(1-rhoi^2)*x2
   cor(x1, x3)
}

set.seed(1234)
k  <- 200                               ### number of studies to simulate
ni <- round(rchisq(k, df=2) * 20 + 20)  ### simulate sample sizes (skewed distribution)
ri <- mapply(gencor, rep(0.2,k), ni)    ### simulate correlations

res <- rma(measure="ZCOR", ri=ri, ni=ni, method="EE") ### use r-to-z transformed correlations

funnel(res, yaxis="sei")
funnel(res, yaxis="vi")
funnel(res, yaxis="seinv")
funnel(res, yaxis="vinv")
funnel(res, yaxis="ni")
funnel(res, yaxis="ninv")
funnel(res, yaxis="sqrtni")
funnel(res, yaxis="sqrtninv")
funnel(res, yaxis="lni")
funnel(res, yaxis="wi")
}
\keyword{hplot}
