\name{influence.rma.mv}
\alias{influence.rma.mv}
\alias{cooks.distance.rma.mv}
\alias{hatvalues.rma.mv}
\title{Outlier and Influential Case Diagnostics for 'rma.mv' Objects}
\description{
   The functions compute various outlier and influential case diagnostics (some of which indicate the influence of deleting one case/study at a time on the model fit and the fitted/residual values) for objects of class \code{"rma.mv"}. At the moment, only Cook's distances and hat values can be obtained.
}
\usage{
\method{cooks.distance}{rma.mv}(model, progbar=FALSE, cluster,
               reestimate=TRUE, parallel="no", ncpus=1, cl=NULL, \dots)
\method{hatvalues}{rma.mv}(model, type="diagonal", \dots)
}
\arguments{
   \item{model}{an object of class \code{"rma.mv"}.}
   \item{progbar}{logical indicating whether a progress bar should be shown (the default is \code{FALSE}). Ignored when using parallel processing.}
   \item{cluster}{optional vector specifying a clustering variable to use for computing the Cook's distances. If not specified, Cook's distances are computed for all individual observed outcomes.}
   \item{reestimate}{logical indicating whether variance/correlation components should be re-estimated after deletion of the \eqn{i}th study/cluster (the default is \code{TRUE}).}
   \item{parallel}{character string indicating whether parallel processing should be used (the default is \code{"no"}). For parallel processing, set to either \code{"snow"} or \code{"multicore"}. See \sQuote{Details}.}
   \item{ncpus}{integer specifying the number of processes to be used in the parallel processing.}
   \item{cl}{optional snow cluster to use if \code{parallel="snow"}. If not supplied, a cluster on the local machine is created for the duration of the call.}
   \item{type}{character string indicating whether to return only the diagonal of the hat matrix (\code{"diagonal"}) or the entire hat matrix (\code{"matrix"}).}
   \item{\dots}{other arguments.}
}
\details{
   Cook's distance for the \eqn{i}th study/cluster can be interpreted as the Mahalanobis distance between the entire set of predicted values once with the \eqn{i}th study/cluster included and once with the \eqn{i}th study/cluster excluded from the model fitting.
}
\value{
   The \code{cooks.distance} function returns a vector. The \code{hatvalues} function returns either a vector with the diagonal elements of the hat matrix or the entire hat matrix.
}
\note{
   Right now, leave-one-out diagnostics (at this moment, only Cook's distance) are calculated by refitting the model \eqn{k} times (where \eqn{k} is the number of studies/clusters). Depending on how large \eqn{k} is, it may take a few moments to finish the calculations. On machines with multiple cores, one can usually speed things up by delegating the model fitting to separate worker processes, that is, by setting \code{parallel="snow"} or \code{parallel="multicore"} and \code{ncpus} to some value larger than 1. Parallel processing makes use of the \code{\link[parallel]{parallel}} package, using the \code{\link[parallel]{makePSOCKcluster}} and \code{\link[parallel]{parLapply}} functions when \code{parallel="snow"} or using \code{\link[parallel]{mclapply}} when \code{parallel="multicore"} (the latter only works on Unix/Linux-alikes). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.

   Alternatively (or in addition to using parallel processing), one can also set \code{reestimate=FALSE}, in which case any variance/correlation components in the model are not re-estimated after deleting the \eqn{i}th study/cluster from the dataset. Doing so only yields an approximation to the Cook's distances that ignores the influence of the \eqn{i}th study/cluster on the variance/correlation components, but is considerably faster (and often yields similar results).

   It may not be possible to fit the model after deletion of the \eqn{i}th study/cluster from the dataset. This will result in \code{NA} values for that study/cluster.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). \emph{Regression diagnostics}. New York: Wiley.

   Cook, R. D., & Weisberg, S. (1982). \emph{Residuals and influence in regression}. London: Chapman and Hall.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Viechtbauer, W., & Cheung, M. W.-L. (2010). Outlier and influence diagnostics for meta-analysis. \emph{Research Synthesis Methods}, \bold{1}, 112--125.
}
\examples{
### load data
dat <- get(data(dat.konstantopoulos2011))

### multilevel random-effects model
res <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat)
print(res, digits=3)

### Cook's distances for each observed outcome
x <- cooks.distance(res)
x
plot(x, type="o", pch=19, xlab="Observed Outcome", ylab="Cook's Distance")

### Cook's distances for each district
x <- cooks.distance(res, cluster=dat$district)
x
plot(x, type="o", pch=19, xlab="District", ylab="Cook's Distance", xaxt="n")
axis(side=1, at=seq_along(x), labels=as.numeric(names(x)))

### hat values
hatvalues(res)
}
\keyword{models}
