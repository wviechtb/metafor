\name{cumul}
\alias{cumul}
\alias{cumul.rma.uni}
\alias{cumul.rma.mh}
\alias{cumul.rma.peto}
\title{Cumulative Meta-Analysis for 'rma' Objects}
\description{
   Function to carry out a \sQuote{cumulative meta-analysis}, by repeatedly fitting the specified model adding one study at a time. \loadmathjax
}
\usage{
cumul(x, \dots)

\method{cumul}{rma.uni}(x, order, digits, transf, targs, collapse=FALSE, progbar=FALSE, \dots)
\method{cumul}{rma.mh}(x, order, digits, transf, targs, collapse=FALSE, progbar=FALSE, \dots)
\method{cumul}{rma.peto}(x, order, digits, transf, targs, collapse=FALSE, progbar=FALSE, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, or \code{"rma.peto"}.}
   \item{order}{optional argument to specify a variable based on which the studies will be ordered for the cumulative meta-analysis.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{transf}{optional argument to specify a function to transform the model coefficients and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{collapse}{logical to specify whether studies with the same value of the \code{order} variable should be added simultaneously (the default is \code{FALSE}).}
   \item{progbar}{logical to specify whether a progress bar should be shown (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   For \code{"rma.uni"} objects, the model specified via \code{x} must be a model without moderators (i.e., either an equal- or a random-effects model).

   If argument \code{order} is not specified, the studies are added according to their order in the original dataset.

   When a variable is specified for \code{order}, the variable is assumed to be of the same length as the original dataset that was used in the model fitting (and if the \code{data} argument was used in the original model fit, then the variable will be searched for within this data frame first). Any subsetting and removal of studies with missing values that was applied during the model fitting is also automatically applied to the variable specified via the \code{order} argument.

   By default, studies are added one at a time. However, if a variable is specified for the \code{order} argument and \code{collapse=TRUE}, then studies with the same value of the \code{order} variable are added simultaneously.
}
\value{
   An object of class \code{c("list.rma","cumul.rma")}. The object is a list containing the following components:
   \item{k}{number of studies included in the analysis.}
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
   \item{order}{values of the \code{order} variable (if specified).}
   \item{\dots}{other arguments.}

   When the model was fitted with \code{test="t"}, \code{test="knha"}, \code{test="hksj"}, or \code{test="adhoc"}, then \code{zval} is called \code{tval} in the object that is returned by the function.

   The object is formatted and printed with the \code{\link[=print.list.rma]{print}} function. To format the results as a data frame, one can use the \code{\link[=as.data.frame.list.rma]{as.data.frame}} function. A forest plot showing the results from the cumulative meta-analysis can be obtained with \code{\link[=forest.cumul.rma]{forest}}. Alternatively, \code{\link[=plot.cumul.rma]{plot}} can also be used to visualize the results.
}
\note{
   When using the \code{transf} option, the transformation is applied to the estimated coefficients and the corresponding interval bounds. The standard errors are then set equal to \code{NA} and are omitted from the printed output.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Chalmers, T. C., & Lau, J. (1993). Meta-analytic stimulus for changes in clinical trials. \emph{Statistical Methods in Medical Research}, \bold{2}(2), 161--172. \verb{https://doi.org/10.1177/096228029300200204}

   Lau, J., Schmid, C. H., & Chalmers, T. C. (1995). Cumulative meta-analysis of clinical trials builds evidence for exemplary medical care. \emph{Journal of Clinical Epidemiology}, \bold{48}(1), 45--57. \verb{https://doi.org/10.1016/0895-4356(94)00106-z}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link[=forest.cumul.rma]{forest}} for a function to draw cumulative forest plots and \code{\link[=plot.cumul.rma]{plot}} for a different visualization of the cumulative results.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
              data=dat.bcg, slab=paste0(author, ", ", year))

### fit random-effects model
res <- rma(yi, vi, data=dat, digits=3)

### cumulative meta-analysis (in the order of publication year)
cumul(res, order=year)
cumul(res, order=year, transf=exp)

### add studies with the same publication year simultaneously
cumul(res, order=year, transf=exp, collapse=TRUE)

### meta-analysis of the (log) risk ratios using the Mantel-Haenszel method
res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
              data=dat.bcg, slab=paste0(author, ", ", year), digits=3)

### cumulative meta-analysis
cumul(res, order=year)
cumul(res, order=year, transf=exp)

### add studies with the same publication year simultaneously
cumul(res, order=year, transf=exp, collapse=TRUE)

### meta-analysis of the (log) odds ratios using Peto's method
res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg,
                data=dat.bcg, slab=paste0(author, ", ", year), digits=3)

### cumulative meta-analysis
cumul(res, order=year)
cumul(res, order=year, transf=exp)

### add studies with the same publication year simultaneously
cumul(res, order=year, transf=exp, collapse=TRUE)

### make the first log risk ratio missing and fit the model without study 2;
### then the variable specified via 'order' should still be of the same length
### as the original dataset; subsetting and removal of studies with missing
### values is automatically done by the cumul() function
dat$yi[1] <- NA
res <- rma(yi, vi, data=dat, subset=-2, digits=3)
cumul(res, transf=exp, order=year)
}
\keyword{methods}
