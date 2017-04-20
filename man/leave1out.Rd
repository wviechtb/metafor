\name{leave1out}
\alias{leave1out}
\alias{leave1out.rma.uni}
\alias{leave1out.rma.mh}
\alias{leave1out.rma.peto}
\title{Leave-One-Out Diagnostics for 'rma' Objects}
\description{
   The functions repeatedly fit the specified model, leaving out one observation/study at a time.
}
\usage{
leave1out(x, \dots)

\method{leave1out}{rma.uni}(x, digits, transf, targs, progbar=FALSE, \dots)
\method{leave1out}{rma.mh}(x, digits, transf, targs, progbar=FALSE, \dots)
\method{leave1out}{rma.peto}(x, digits, transf, targs, progbar=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.mh"}, \code{"rma.peto"}, or \code{"rma.uni"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{transf}{an optional argument specifying the name of a function that should be used to transform the model coefficients and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{progbar}{logical indicating whether a progress bar should be shown (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   The model specified via \code{x} must be a model without moderators (i.e., either a fixed- or a random-effects model and not a fixed-effects with moderators or mixed-effects model).
}
\value{
   An object of class \code{"list.rma"}. The object is a list containing the following components:
   \item{estimate}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{p-values for the test statistics.}
   \item{ci.lb}{lower bounds of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bounds of the confidence intervals for the coefficients.}
   \item{Q}{test statistics for the tests of heterogeneity.}
   \item{Qp}{p-values for the tests of heterogeneity.}
   \item{tau2}{estimated amounts of (residual) heterogeneity (only for random-effects models).}
   \item{I2}{values of \if{latex}{\eqn{I^2}}\if{html}{\eqn{I²}}\if{text}{\eqn{I^2}} (only for random-effects models).}
   \item{H2}{values of \if{latex}{\eqn{H^2}}\if{html}{\eqn{H²}}\if{text}{\eqn{H^2}} (only for random-effects models).}

   The \code{"list.rma"} object is formated and printed with \code{\link{print.list.rma}}.
}
\note{
   When using the \code{transf} option, the transformation is applied to the estimated coefficients and the corresponding interval bounds. The standard errors are then set equal to \code{NA} and are omitted from the printed output.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}, 112--125.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

### cumulative meta-analysis
leave1out(res)
leave1out(res, transf=exp)

### meta-analysis of the (log) risk ratios using the Mantel-Haenszel method
res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### cumulative meta-analysis
leave1out(res)
leave1out(res, transf=exp)

### meta-analysis of the (log) odds ratios using Peto's method
res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### cumulative meta-analysis
leave1out(res)
leave1out(res, transf=exp)
}
\keyword{methods}
