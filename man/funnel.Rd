\name{funnel}
\alias{funnel}
\alias{funnel.rma}
\alias{funnel.default}
\title{Funnel Plots}
\description{Function to create funnel plots.}
\usage{
funnel(x, \dots)

\method{funnel}{rma}(x, yaxis="sei",
       xlim, ylim, xlab, ylab,
       steps=5, at, atransf, targs, digits, level=x$level,
       addtau2=FALSE, type="rstandard",
       back="lightgray", shade="white", hlines="white",
       refline, pch=19, pch.fill=21, col, bg, ci.res=1000, \dots)

\method{funnel}{default}(x, vi, sei, ni, subset, yaxis="sei",
       xlim, ylim, xlab, ylab,
       steps=5, at, atransf, targs, digits, level=95,
       back="lightgray", shade="white", hlines="white",
       refline=0, pch=19, col, bg, ci.res=1000, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"} or a vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{sei}{vector with the corresponding standard errors.}
   \item{ni}{vector with the corresponding sample sizes.}
   \item{subset}{optional vector indicating the subset of studies that should be included in the plot. This can be a logical vector of the same length as \code{x} or a numeric vector indicating the indices of the observations to include. Only relevant when passing a vector via \code{x}.}
   \item{yaxis}{either \code{"sei"}, \code{"vi"}, \code{"seinv"}, \code{"vinv"}, \code{"ni"}, \code{"ninv"}, \code{"sqrtni"}, \code{"sqrtninv"}, \code{"lni"}, or \code{"wi"} to indicate what values should be placed on the y-axis. See \sQuote{Details}.}
   \item{xlim}{x-axis limits. If unspecified, the function tries to set the x-axis limits to some sensible values.}
   \item{ylim}{y-axis limits. If unspecified, the function tries to set the y-axis limits to some sensible values.}
   \item{xlab}{title for the x-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{ylab}{title for the y-axis. If unspecified, the function tries to set an appropriate axis title.}
   \item{steps}{the number of tick marks for the y-axis (the default is 5).}
   \item{at}{position of the x-axis tick marks and corresponding labels. If unspecified, the function tries to set the tick mark positions/labels to some sensible values.}
   \item{atransf}{optional argument specifying the name of a function that should be used to transform the x-axis labels (e.g., \code{atransf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified via \code{atransf}.}
   \item{digits}{integer specifying the number of decimal places to which the tick mark labels of the x- and y-axis should be rounded. Can also be a vector of two integers, the first specifying the number of decimal places for the x-axis, the second for the y-axis labels (e.g., \code{c(2,3)}). If unspecified, the function tries to set the argument to some sensible values.}
   \item{level}{numerical value between 0 and 100 specifying the level of the pseudo confidence interval region (for \code{"rma"} objects, the default is to take the value from the object). May also be a vector of values to obtain multiple regions. See \sQuote{Examples}.}
   \item{addtau2}{logical to indicate whether the amount of heterogeneity should be accounted for when drawing the pseudo confidence interval region (the default is \code{FALSE}). Ignored when the model includes moderators and residuals are plotted. See \sQuote{Details}.}
   \item{type}{either \code{"rstandard"} (default) or \code{"rstudent"} indicating whether the usual or deleted residuals should be used in creating the funnel plot when the model involves moderators. See \sQuote{Details}.}
   \item{back}{color to use for the background of the plotting region.}
   \item{shade}{color to use for shading the pseudo confidence interval region. When \code{level} is a vector of values, different shading colors can be specified for each region.}
   \item{hlines}{color of the horizontal reference lines.}
   \item{refline}{value at which to draw the vertical reference line and, if drawn, where the pseudo confidence interval should be centered. If unspecified, the reference line is drawn at the fixed- or random-effects model estimate when the model does not include moderators and at zero when moderators are included (and therefore residuals are plotted) or when directly plotting observed effect sizes or outcomes.}
   \item{pch}{plotting symbol to use for the observed effect sizes or outcomes. By default, a solid circle is used. Can also be a vector of values. See \code{\link{points}} for other options.}
   \item{pch.fill}{plotting symbol to use for the effect sizes or outcomes filled in by the trim and fill method. By default, a circle is used. Only relevant when plotting an object created by the \code{\link{trimfill}} function.}
   \item{col}{optional character string specifying the name of a color to use for the points ("black" is used by default if not specified). Can also be a vector of color names.}
   \item{bg}{optional character string specifying the name of a background color for open plot symbols ("white" is used by default if not specified). Can also be a vector of color names.}
   \item{ci.res}{integer specifying the number of y-axis values at which to calculate the bounds of the pseudo confidence interval. The default is \code{1000}, which usually provides a sufficient resolution for the plotting.}
   \item{\dots}{other arguments.}
}
\details{
   For fixed- and random-effects models (i.e., models not involving moderators), the plot shows the individual observed effect sizes or outcomes on the x-axis against the corresponding standard errors (i.e., the square root of the sampling variances) on the y-axis. A vertical line indicates the estimate based on the model. A pseudo confidence interval region is drawn around this value with bounds equal to \if{latex}{\eqn{\pm 1.96 SE}}\if{html}{\eqn{± 1.96 SE}}\if{text}{\eqn{+- 1.96 SE}}, where \eqn{SE} is the standard error value from the y-axis (assuming \code{level=95}). If \code{addtau2=TRUE} (only for models of class \code{"rma.uni"}), then the bounds of the pseudo confidence interval region are equal to \if{latex}{\eqn{\pm 1.96 \sqrt{SE^2 + \tau^2}}}\if{html}{\eqn{± 1.96 \sqrt(SE² + \tau²)}}\if{text}{\eqn{+- 1.96 \sqrt(SE^2 + \tau^2)}}, where \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} is the amount of heterogeneity as estimated by the model.

   For models involving moderators, the plot shows the residuals on the x-axis against their corresponding standard errors. Either the usual or deleted residuals can be used for that purpose (set via the \code{type} argument). See \code{\link{residuals.rma}} for more details on the different types of residuals.

   With the \code{atransf} argument, the labels of the observed effect sizes or outcomes on the x-axis can be transformed with some suitable function. For example, when plotting log odds ratios, one could use \code{transf=exp} to obtain a funnel plot with the values on the x-axis corresponding to the odds ratios. See also \link{transf} for some transformation functions useful for meta-analyses.

   Instead of placing the standard error value on the y-axis, several other options are available by setting the \code{yaxis} argument to:
   \itemize{
   \item \code{yaxis="vi"} for the sampling variance,
   \item \code{yaxis="seinv"} for the inverse of the standard error,
   \item \code{yaxis="vinv"} for the inverse of the sampling variance,
   \item \code{yaxis="ni"} for the sample size,
   \item \code{yaxis="ninv"} for the inverse of the sample size,
   \item \code{yaxis="sqrtni"} for the square root sample size,
   \item \code{yaxis="sqrtninv"} for the inverse of the square root sample size,
   \item \code{yaxis="lni"} for the log of the sample size,
   \item \code{yaxis="wi"} for the weights.
   } However, only when \code{yaxis="sei"} (the default) will the pseudo confidence region have the expected (upside-down) funnel shape with straight lines. Also, when placing (a function of) the sample size on the y-axis or the weights, then the pseudo confidence region cannot be drawn. See Sterne and Egger (2001) for more details on the choice of the y-axis.

   If the object passed to the function comes from the \code{\link{trimfill}} function, the effect sizes or outcomes that are filled in by the trim and fill method are also added to the funnel plot. The symbol to use for plotting the filled in values can then be specified via the \code{pch.fill} argument.

   One can also directly pass a vector of observed effect sizes or outcomes (via \code{x}) and the corresponding sampling variances (via \code{vi}), standard errors (via \code{sei}), and/or sample sizes (via \code{ni}) to the function. By default, the vertical reference line is then drawn at zero.

   The arguments \code{back}, \code{shade}, and \code{hlines} can be set to \code{NULL} to suppress the shading and the horizontal reference lines.
}
\note{
   Placing (a function of) the sample size on the y-axis (i.e., using \code{yaxis="ni"}, \code{yaxis="ninv"}, \code{yaxis="sqrtni"}, \code{yaxis="sqrtninv"}, or \code{yaxis="lni"}) is only possible when information about the sample sizes is actually stored within the object passed to the \code{funnel} function. That should automatically be the case when the observed outcomes were computed with the \code{escalc} function or when the observed outcomes were computed within the model fitting function. On the other hand, this will not automatically be the case when \code{rma.uni} was used together with the \code{yi} and \code{vi} arguments and the \code{yi} and \code{vi} values were \emph{not} computed with \code{escalc}. In that case, it is still possible to pass information about the sample sizes to the \code{rma.uni} function (i.e., use \code{rma.uni(yi, vi, ni)}).

   When using unweighted estimation, using \code{yaxis="wi"} will place all points on a horizontal line. When directly passing a vector of observed effect sizes or outcomes to the function, \code{yaxis="wi"} is equivalent to \code{yaxis="vinv"}, except that the weights are expressed in percent.

   When specifying vectors for \code{pch}, \code{col}, and/or \code{bg}, the variables specified are assumed to be of the same length as the data passed to the funnel function or the model fitting function (when using \code{funnel} on a model object). Any subsetting and removal of studies with missing values is automatically applied to the variables specified via \code{pch}, \code{col}, and \code{bg}.
}
\value{
   A data frame with components:
   \item{x}{the x coordinates of the points that were plotted.}
   \item{y}{the y coordinates of the points that were plotted.}
   \item{slab}{study labels of the points that were plotted.}
   Note that the data frame is returned invisibly.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Light, R. J., & Pillemer, D. B. (1984). \emph{Summing up: The science of reviewing research}. Cambridge, MA: Harvard University Press.

   Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2008). Contour-enhanced meta-analysis funnel plots help distinguish publication bias from other causes of asymmetry. \emph{Journal of Clinical Epidemiology}, \bold{61}, 991--996.

   Sterne, J. A. C., & Egger, M. (2001). Funnel plots for detecting bias in meta-analysis: Guidelines on choice of axis. \emph{Journal of Clinical Epidemiology}, \bold{54}, 1046--1055.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}, \code{\link{influence.rma.uni}}, \code{\link{trimfill}}
}
\examples{
### copy BCG vaccine data into 'dat'
dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### random-effects model
res <- rma(yi, vi, data=dat)

### standard funnel plot
funnel(res)

### show risk ratio values on x-axis (log scale)
funnel(res, atransf=exp)

### passing log risk ratios and sampling variances directly to the function
### note: essentially the same plot, except that reference line is centered at zero
funnel(dat$yi, dat$vi)
funnel(res, refline=0)

### funnel plot with risk ratio values on the x-axis (log scale)
funnel(res, atransf=exp, at=log(c(.12, .25, .5, 1, 2)))

### contour-enhanced funnel plot centered at 0 (see Peters et al., 2008)
funnel(res, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0)

### same, but show risk ratio values on the x-axis
funnel(res, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"),
       refline=0, atransf=exp, at=log(c(.10, .25, .5, 1, 2, 4, 10)))

### illustrate the use of vectors for 'pch' and 'col'
res <- rma(yi, vi, data=dat, subset=2:10)
funnel(res, pch=ifelse(dat$yi > -1, 19, 21), col=ifelse(sqrt(dat$vi) > .3, "red", "blue"))

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

res <- rma(measure="ZCOR", ri=ri, ni=ni, method="FE") ### use r-to-z transformed correlations

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
