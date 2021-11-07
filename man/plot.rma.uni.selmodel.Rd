\name{plot.rma.uni.selmodel}
\alias{plot.rma.uni.selmodel}
\title{Plot Method for 'plot.rma.uni.selmodel' Objects}
\description{
   Plot method for objects of class \code{"plot.rma.uni.selmodel"}. \loadmathjax
}
\usage{
\method{plot}{rma.uni.selmodel}(x, xlim, ylim, n=1000, prec="max", scale=FALSE,
      ci=FALSE, reps=1000, rug=TRUE, add=FALSE,
      lty=c("solid","dotted"), lwd=c(2,1), \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni.selmodel"} obtained with \code{\link{selmodel}}.}
   \item{xlim}{x-axis limits. Essentially the range of p-values for which the selection function should be drawn. If unspecified, the function sets the limits automatically.}
   \item{ylim}{y-axis limits. If unspecified, the function sets the limits automatically.}
   \item{n}{numeric value to specify for how many p-values within the x-axis limits the function value should be computed (the default is 1000).}
   \item{prec}{either a character string (with options \code{"max"}, \code{"min"}, \code{"mean"}, or \code{"median"}) or a numeric value. See \sQuote{Details}.}
   \item{scale}{logical to specify whether the function values should be rescaled to a 0 to 1 range (the default is \code{FALSE}).}
   \item{ci}{logical to specify whether a confidence interval should be drawn around the selection function (the default is \code{FALSE}). Can also be a string (with options \code{"boot"} or \code{"wald"}). See \sQuote{Details}.}
   \item{reps}{numeric value to specify the number of bootstrap samples to draw for generating the confidence interval bounds (the default is 1000).}
   \item{rug}{logical to specify whether the observed p-values should be added as tick marks on the x-axis (the default is \code{TRUE}).}
   \item{add}{logical to specify whether the function should be added to an existing plot (the default is \code{FALSE}).}
   \item{lty}{the line types for the selection function and the confidence interval bounds.}
   \item{lwd}{the line widths for the selection function and the confidence interval bounds.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to draw the estimated selection function based on objects of class \code{"plot.rma.uni.selmodel"}.

   When the selection function incorporates a measure of precision (which, strictly speaking, is really a measure of imprecision), one can specify for which level of precision the selection function should be drawn. When \code{prec="max"}, then the function is drawn for the \emph{least} precise study (maximum imprecision), when \code{prec="min"}, then the function is drawn for the \emph{most} precise study (minimum imprecision), while \code{prec="mean"} and \code{prec="median"} will show the function for the mean and median level of imprecision, respectively. Alternatively, one can specify a numeric value for argument \code{prec} to specify the precision value (where \code{prec="max"} corresponds to \code{prec=1} and higher levels of precision to \code{prec} values below 1).

   When \code{ci=TRUE} (or equivalently, \code{ci="boot"}), a confidence interval is drawn around the selection function. The bounds of this interval are generated using parametric bootstrapping, with argument \code{reps} controlling the number of bootstrap samples to draw for generating the confidence interval bounds. When both \code{n} and \code{reps} are large, constructing the confidence interval can take a few moments to complete.

   For models where the selection function involves a single \mjseqn{\delta} parameter, one can also set \code{ci="wald"}, in which case the confidence interval will be constructed based on the Wald-type CI of the \mjseqn{\delta} parameter (doing so is much quicker than using parametric bootstrapping). This option is also available for step function models (even if they involve multiple \mjseqn{\delta} parameters).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{selmodel.rma.uni}}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.hackshaw1998

### fit random-effects model using the log odds ratios
res <- rma(yi, vi, data=dat, method="ML")
res

### fit step selection model
sel1 <- selmodel(res, type="stepfun", steps=c(0.05, 0.10, 0.50, 1.00))

### plot selection function
plot(sel1, scale=TRUE)

### fit negative exponential selection model
sel2 <- selmodel(res, type="negexp")

### add selection function to the existing plot
plot(sel2, add=TRUE, col="blue")

### plot selection function with CI
plot(sel1, ci="wald")
}
\keyword{hplot}
