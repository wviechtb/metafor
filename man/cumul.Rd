\name{cumul}
\alias{cumul}
\alias{cumul.rma.uni}
\alias{cumul.rma.mh}
\alias{cumul.rma.peto}
\title{Cumulative Meta-Analysis for 'rma' Objects}
\description{
   The functions repeatedly fit the specified model, adding one observation/study at a time to the model. \loadmathjax
}
\usage{
cumul(x, \dots)

\method{cumul}{rma.uni}(x, order, digits, transf, targs, progbar=FALSE, \dots)
\method{cumul}{rma.mh}(x, order, digits, transf, targs, progbar=FALSE, \dots)
\method{cumul}{rma.peto}(x, order, digits, transf, targs, progbar=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.mh"}, \code{"rma.peto"}, \code{"rma.uni"}.}
   \item{order}{optional vector with indices giving the desired order for the cumulative meta-analysis.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the model coefficients and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{progbar}{logical indicating whether a progress bar should be shown (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   For \code{"rma.uni"} objects, the model specified by \code{x} must be a model without moderators (i.e., either a fixed- or a random-effects model).
}
\value{
   An object of class \code{c("list.rma","cumul.rma")}. The object is a list containing the following components:
   \item{estimate}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{p-values for the test statistics.}
   \item{ci.lb}{lower bounds of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bounds of the confidence intervals for the coefficients.}
   \item{QE}{test statistics for the tests of heterogeneity.}
   \item{QEp}{p-values for the tests of heterogeneity.}
   \item{tau2}{estimated amounts of (residual) heterogeneity (only for random-effects models).}
   \item{I2}{values of \mjseqn{I^2} (only for random-effects models).}
   \item{H2}{values of \mjseqn{H^2} (only for random-effects models).}

   The object is formatted and printed with \code{\link{print.list.rma}}. A forest plot showing the results from the cumulative meta-analysis can be obtained with \code{\link{forest.cumul.rma}}. For random-effects models, \code{\link{plot.cumul.rma}} can also be used to visualize the results.
}
\note{
   When using the \code{transf} option, the transformation is applied to the estimated coefficients and the corresponding interval bounds. The standard errors are then set equal to \code{NA} and are omitted from the printed output.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Chalmers, T. C., & Lau, J. (1993). Meta-analytic stimulus for changes in clinical trials. \emph{Statistical Methods in Medical Research}, \bold{2}(2), 161--172. \url{https://doi.org/10.1177/096228029300200204}

   Lau, J., Schmid, C. H., & Chalmers, T. C. (1995). Cumulative meta-analysis of clinical trials builds evidence for exemplary medical care. \emph{Journal of Clinical Epidemiology}, \bold{48}(1), 45--57. \url{https://doi.org/10.1016/0895-4356(94)00106-z}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{forest.cumul.rma}}, \code{\link{plot.cumul.rma}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

### cumulative meta-analysis (in the order of publication year)
cumul(res, transf=exp, order=order(dat$year))

### meta-analysis of the (log) risk ratios using the Mantel-Haenszel method
res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### cumulative meta-analysis
cumul(res, order=order(dat.bcg$year))
cumul(res, order=order(dat.bcg$year), transf=TRUE)

### meta-analysis of the (log) odds ratios using Peto's method
res <- rma.mh(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### cumulative meta-analysis
cumul(res, order=order(dat.bcg$year))
cumul(res, order=order(dat.bcg$year), transf=TRUE)
}
\keyword{methods}
