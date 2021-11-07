\name{residuals.rma}
\alias{residuals}
\alias{rstandard}
\alias{rstudent}
\alias{residuals.rma}
\alias{rstandard.rma.uni}
\alias{rstandard.rma.mh}
\alias{rstandard.rma.mv}
\alias{rstandard.rma.peto}
\alias{rstudent.rma.uni}
\alias{rstudent.rma.mh}
\alias{rstudent.rma.mv}
\alias{rstudent.rma.peto}
\title{Residual Values based on 'rma' Objects}
\description{
   The \code{residuals}, \code{rstandard}, and \code{rstudent} functions compute residuals, corresponding standard errors, and standardized residuals for models fitted with the \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.mv}} functions. \loadmathjax
}
\usage{
\method{residuals}{rma}(object, type="response", \dots)

\method{rstandard}{rma.uni}(model, digits, type="marginal", \dots)
\method{rstandard}{rma.mh}(model, digits, \dots)
\method{rstandard}{rma.peto}(model, digits, \dots)
\method{rstandard}{rma.mv}(model, digits, cluster, \dots)

\method{rstudent}{rma.uni}(model, digits, progbar=FALSE, \dots)
\method{rstudent}{rma.mh}(model, digits, progbar=FALSE, \dots)
\method{rstudent}{rma.peto}(model, digits, progbar=FALSE, \dots)
\method{rstudent}{rma.mv}(model, digits, progbar=FALSE, cluster,
         reestimate=TRUE, parallel="no", ncpus=1, cl, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"} (for \code{residuals}).}
   \item{type}{the type of residuals which should be returned. For \code{residuals}, the alternatives are: \code{"response"} (default), \code{"rstandard"}, \code{"rstudent"}, and \code{"pearson"}. For \code{rstandard.rma.uni}, the alternatives are: \code{"marginal"} (default) and \code{"conditional"}. See \sQuote{Details}.}
   \item{model}{an object of class \code{"rma"} (for \code{residuals}) or an object of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, or \code{"rma.mv"} (for \code{rstandard} and \code{rstudent}).}
   \item{cluster}{optional vector to specify a clustering variable to use for computing cluster-level multivariate standardized residuals (only for \code{"rma.mv"} objects).}
   \item{reestimate}{logical to specify whether variance/correlation components should be re-estimated after deletion of the \mjeqn{i\textrm{th}}{ith} case when computing externally standardized residuals for \code{"rma.mv"} objects (the default is \code{TRUE}).}
   \item{parallel}{character string to specify whether parallel processing should be used (the default is \code{"no"}). For parallel processing, set to either \code{"snow"} or \code{"multicore"}. See \sQuote{Details}.}
   \item{ncpus}{integer to specify the number of processes to use in the parallel processing.}
   \item{cl}{optional cluster to use if \code{parallel="snow"}. If unspecified, a cluster on the local machine is created for the duration of the call.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{progbar}{logical to specify whether a progress bar should be shown (only for \code{rstudent}) (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   The observed residuals (obtained with \code{residuals}) are simply equal to the \sQuote{observed - fitted} values. These can be obtained with \code{residuals(object)} (using the default \code{type="response"}).

   Dividing the observed residuals by the model-implied standard errors of the observed effect sizes or outcomes yields Pearson (or semi-standardized) residuals. These can be obtained with \code{residuals(object, type="pearson")}.

   Dividing the observed residuals by their corresponding standard errors yields (internally) standardized residuals. These can be obtained with \code{rstandard(model)} or \code{residuals(object, type="rstandard")}.

   With \code{rstudent(model)} (or \code{residuals(object, type="rstudent")}), one can obtain the externally standardized residuals (also called standardized deleted residuals or (externally) studentized residuals). The externally standardized residual for the \mjeqn{i\textrm{th}}{ith} case is obtained by deleting the \mjeqn{i\textrm{th}}{ith} case from the dataset, fitting the model based on the remaining cases, calculating the predicted value for the \mjeqn{i\textrm{th}}{ith} case based on the fitted model, taking the difference between the observed and the predicted value for the \mjeqn{i\textrm{th}}{ith} case (which yields the deleted residual), and then standardizing the deleted residual based on its standard error.

   If a particular case fits the model, its standardized residual follows (asymptotically) a standard normal distribution. A large standardized residual for a case therefore may suggest that the case does not fit the assumed model (i.e., it may be an outlier).

   For \code{"rma.uni"} objects, \code{rstandard(model, type="conditional")} computes conditional residuals, which are the deviations of the observed effect sizes or outcomes from the best linear unbiased predictions (BLUPs) of the study-specific true effect sizes or outcomes (see \code{\link[=blup.rma.uni]{blup}}).

   For \code{"rma.mv"} objects, one can specify a clustering variable (via the \code{cluster} argument). If specified, \code{rstandard(model)} and \code{rstudent(model)} also compute cluster-level multivariate (internally or externally) standardized residuals. If all outcomes within a cluster fit the model, then the multivariate standardized residual for the cluster follows (asymptotically) a chi-square distribution with \mjseqn{k_i} degrees of freedom (where \mjseqn{k_i} denotes the number of outcomes within the cluster).

   See also \code{\link{influence.rma.uni}} and \code{\link{influence.rma.mv}} for other leave-one-out diagnostics that are useful for detecting influential cases in models fitted with the \code{\link{rma.uni}} and \code{\link{rma.mv}} functions.
}
\value{
   Either a vector with the residuals of the requested type (for \code{residuals}) or an object of class \code{"list.rma"}, which is a list containing the following components:
   \item{resid}{observed residuals (for \code{rstandard}) or deleted residuals (for \code{rstudent}).}
   \item{se}{corresponding standard errors.}
   \item{z}{standardized residuals (internally standardized for \code{rstandard} or externally standardized for \code{rstudent}).}

   When a clustering variable is specified for \code{"rma.mv"} objects, the returned object is a list with the first element (named \code{obs}) as described above and a second element (named \code{cluster} of class \code{"list.rma"} with:
   \item{X2}{cluster-level multivariate standardized residuals.}
   \item{k}{number of observed effect sizes or outcomes within the clusters.}

   The object is formatted and printed with \code{\link{print.list.rma}}.
}
\note{
   Right now, the externally standardized residuals (obtained with \code{rstudent}) are calculated by refitting the model \mjseqn{k} times (where \mjseqn{k} is the number of cases). Depending on how large \mjseqn{k} is, it may take a few moments to finish the calculations. For complex models fitted with \code{\link{rma.mv}}, this can become computationally expensive.

   On machines with multiple cores, one can usually speed things up by delegating the model fitting to separate worker processes, that is, by setting \code{parallel="snow"} or \code{parallel="multicore"} and \code{ncpus} to some value larger than 1 (only for objects of class \code{"rma.mv"}). Parallel processing makes use of the \code{\link[parallel]{parallel}} package, using the \code{\link[parallel]{makePSOCKcluster}} and \code{\link[parallel]{parLapply}} functions when \code{parallel="snow"} or using \code{\link[parallel]{mclapply}} when \code{parallel="multicore"} (the latter only works on Unix/Linux-alikes). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.

   Alternatively (or in addition to using parallel processing), one can also set \code{reestimate=FALSE}, in which case any variance/correlation components in the model are not re-estimated after deleting the \mjeqn{i\textrm{th}}{ith} case from the dataset. Doing so only yields an approximation to the externally standardized residuals (and the cluster-level multivariate standardized residuals) that ignores the influence of the \mjeqn{i\textrm{th}}{ith} case on the variance/correlation components, but is considerably faster (and often yields similar results).

   It may not be possible to fit the model after deletion of the \mjeqn{i\textrm{th}}{ith} case from the dataset. This will result in \code{NA} values for that case when calling \code{rstudent}.

   Also, for \code{"rma.mv"} objects with a clustering variable specified, it may not be possible to compute the cluster-level multivariate standardized residual for a particular cluster (if the var-cov matrix of the residuals within a cluster is not of full rank). This will result in \code{NA} for that cluster.

   For objects of class \code{"rma.mh"} and \code{"rma.peto"}, \code{rstandard} actually computes Pearson (or semi-standardized) residuals.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}(2), 112--125. \verb{https://doi.org/10.1002/jrsm.11}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}, \code{\link{influence.rma.uni}}, \code{\link{influence.rma.mv}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model
res <- rma(yi, vi, data=dat)

### compute the studentized residuals
rstudent(res)

### fit mixed-effects model with absolute latitude as moderator
res <- rma(yi, vi, mods = ~ ablat, data=dat)

### compute the studentized residuals
rstudent(res)
}
\keyword{models}
