\name{leave1out}
\alias{leave1out}
\alias{leave1out.rma.uni}
\alias{leave1out.rma.mh}
\alias{leave1out.rma.peto}
\title{Leave-One-Out Diagnostics for 'rma' Objects}
\description{
   Functions to carry out a \sQuote{leave-one-out analysis}, by repeatedly fitting the specified model leaving out one study (or cluster level) at a time. \loadmathjax
}
\usage{
leave1out(x, \dots)

\method{leave1out}{rma.uni}(x, cluster, digits, transf, targs, progbar=FALSE, \dots)
\method{leave1out}{rma.mh}(x, cluster, digits, transf, targs, progbar=FALSE, \dots)
\method{leave1out}{rma.peto}(x, cluster, digits, transf, targs, progbar=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}.}
   \item{cluster}{optional vector to specify a clustering variable.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{transf}{optional argument to specify a function to transform the model coefficients and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{progbar}{logical to specify whether a progress bar should be shown (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   In a leave-one-out analysis, the same model is repeatedly fitted, leaving out one study at a time. By doing so, we can assess how much the results are influenced by each individual study. It is also possible to specify a \code{cluster} variable, in which case each cluster level is left out in turn.

   Note that for \code{"rma.uni"} objects, the model specified via \code{x} must be a model without moderators (i.e., either an equal- or a random-effects model).
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

   When the model was fitted with \code{test="t"}, \code{test="knha"}, \code{test="hksj"}, or \code{test="adhoc"}, then \code{zval} is called \code{tval} in the object that is returned by the function.

   The object is formatted and printed with the \code{\link[=print.list.rma]{print}} function. To format the results as a data frame, one can use the \code{\link[=as.data.frame.list.rma]{as.data.frame}} function.
}
\note{
   When using the \code{transf} option, the transformation is applied to the estimated coefficients and the corresponding interval bounds. The standard errors are then set equal to \code{NA} and are omitted from the printed output.

   The variable specified via \code{cluster} is assumed to be of the same length as the data originally passed to the model fitting function (and if the \code{data} argument was used in the original model fit, then the variable will be searched for within this data frame first). Any subsetting and removal of studies with missing values that was applied during the model fitting is also automatically applied to the variable specified via the \code{cluster} argument.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W. (2021). Model checking in meta-analysis. In C. H. Schmid, T. Stijnen, & I. R. White (Eds.), \emph{Handbook of meta-analysis} (pp. 219--254). Boca Raton, FL: CRC Press. \verb{https://doi.org/10.1201/9781315119403}

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}(2), 112--125. \verb{https://doi.org/10.1002/jrsm.11}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, and \code{\link{rma.peto}} for functions to fit models for which leave-one-out diagnostics can be computed.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res <- rma(yi, vi, data=dat)

### leave-one-out analysis
leave1out(res)
leave1out(res, transf=exp)

### leave-one-out analysis with a cluster variable
leave1out(res, cluster=alloc)

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
