\name{leave1out}
\alias{leave1out}
\alias{leave1out.rma.uni}
\alias{leave1out.rma.mh}
\alias{leave1out.rma.peto}
\title{Leave-One-Out Diagnostics for 'rma' Objects}
\description{
   The functions repeatedly fit the specified model, leaving out one observation/study at a time. \loadmathjax
}
\usage{
leave1out(x, \dots)

\method{leave1out}{rma.uni}(x, digits, transf, targs, progbar=FALSE, \dots)
\method{leave1out}{rma.mh}(x, digits, transf, targs, progbar=FALSE, \dots)
\method{leave1out}{rma.peto}(x, digits, transf, targs, progbar=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.mh"}, \code{"rma.peto"}, or \code{"rma.uni"}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{transf}{an optional argument to specify a function that should be used to transform the model coefficients and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{progbar}{logical to specify whether a progress bar should be shown (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   For \code{"rma.uni"} objects, the model specified via \code{x} must be a model without moderators (i.e., either an equal- or a random-effects model).
}
\value{
   An object of class \code{"list.rma"}. The object is a list containing the following components:
   \item{estimate}{estimated (average) outcomes.}
   \item{se}{corresponding standard errors.}
   \item{zval}{corresponding test statistics.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{lower bounds of the confidence intervals.}
   \item{ci.ub}{upper bounds of the confidence intervals.}
   \item{Q}{test statistics for the test of heterogeneity.}
   \item{Qp}{corresponding p-values.}
   \item{tau2}{estimated amount of heterogeneity (only for random-effects models).}
   \item{I2}{values of \mjseqn{I^2}.}
   \item{H2}{values of \mjseqn{H^2}.}

   When the model was fitted with \code{test="t"} or \code{test="knha"}, then \code{zval} is called \code{tval} in the object that is returned by the function.

   The object is formatted and printed with the \code{\link[=print.list.rma]{print}} function.
}
\note{
   When using the \code{transf} option, the transformation is applied to the estimated coefficients and the corresponding interval bounds. The standard errors are then set equal to \code{NA} and are omitted from the printed output.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}(2), 112--125. \verb{https://doi.org/10.1002/jrsm.11}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

### leave-one-out analysis
leave1out(res)
leave1out(res, transf=exp)

### meta-analysis of the (log) risk ratios using the Mantel-Haenszel method
res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### leave-one-out analysis
leave1out(res)
leave1out(res, transf=exp)

### meta-analysis of the (log) odds ratios using Peto's method
res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### leave-one-out analysis
leave1out(res)
leave1out(res, transf=exp)
}
\keyword{methods}
