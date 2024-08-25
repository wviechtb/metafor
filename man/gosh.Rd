\name{gosh}
\alias{gosh}
\alias{gosh.rma}
\title{GOSH Plots for 'rma' Objects}
\description{
   Function to create GOSH plots for objects of class \code{"rma"}. \loadmathjax
}
\usage{
gosh(x, \dots)

\method{gosh}{rma}(x, subsets, progbar=TRUE, parallel="no", ncpus=1, cl, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{subsets}{optional integer to specify the number of subsets.}
   \item{progbar}{logical to specify whether a progress bar should be shown (the default is \code{TRUE}).}
   \item{parallel}{character string to specify whether parallel processing should be used (the default is \code{"no"}). For parallel processing, set to either \code{"snow"} or \code{"multicore"}. See \sQuote{Note}.}
   \item{ncpus}{integer to specify the number of processes to use in the parallel processing.}
   \item{cl}{optional cluster to use if \code{parallel="snow"}. If unspecified, a cluster on the local machine is created for the duration of the call.}
   \item{\dots}{other arguments.}
}
\details{
   The model specified via \code{x} must be a model fitted with either the \code{\link{rma.uni}}, \code{\link{rma.mh}}, or \code{\link{rma.peto}} functions.

   Olkin et al. (2012) proposed the GOSH (graphical display of study heterogeneity) plot, which is based on examining the results of an equal-effects model in all possible subsets of size \mjseqn{1, \ldots, k} of the \mjseqn{k} studies included in a meta-analysis. In a homogeneous set of studies, the model estimates obtained this way should form a roughly symmetric, contiguous, and unimodal distribution. On the other hand, when the distribution is multimodal, then this suggests the presence of heterogeneity, possibly due to outliers and/or distinct subgroups of studies. Plotting the estimates against some measure of heterogeneity (e.g., \mjseqn{I^2}, \mjseqn{H^2}, or the \mjseqn{Q}-statistic) can also help to reveal subclusters, which are indicative of heterogeneity. The same type of plot can be produced by first fitting an equal-effects model with either the \code{\link{rma.uni}} (using \code{method="EE"}), \code{\link{rma.mh}}, or \code{\link{rma.peto}} functions and then passing the fitted model object to the \code{gosh} function and then plotting the results.

   For models fitted with the \code{\link{rma.uni}} function (which may be random-effects or mixed-effects meta-regressions models), the idea underlying this type of plot can be generalized (Viechtbauer, 2021) by examining the distribution of all model coefficients, plotting them against each other, and against some measure of (residual) heterogeneity (including the estimate of \mjseqn{\tau^2} or its square root).

   Note that for models without moderators, application of the method requires fitting a total of \mjseqn{2^k - 1} models, which could be an excessively large number when \mjseqn{k} is large. For example, for \mjseqn{k=10}, there are only 1023 possible subsets, but for \mjseqn{k=20}, this number already grows to 1,048,575. For even larger \mjseqn{k}, it may become computationally infeasible to consider all possible subsets. Instead, we can then examine (a sufficiently large number of) random subsets.

   By default, if the number of possible subsets is \mjseqn{\le 10^6}, the function will consider all possible subsets and otherwise \mjseqn{10^6} random subsets. One can use the \code{subsets} argument to specify a different number of subsets to consider. If \code{subsets} is specified and it is actually larger than the number of possible subsets, then the function automatically only considers the possible subsets and does not use random subsets.

   When \code{x} is an equal-effects model or a random-effects model fitted using \code{method="DL"}, provisions have been made to speed up the model fitting to the various subsets. For random-effects models using some other estimator of \mjseqn{\tau^2} (especially an iterative one like \code{method="REML"}), the computations will be considerably slower.
}
\value{
   An object of class \code{"gosh.rma"}. The object is a list containing the following components:
   \item{res}{a data frame with the results for each subset (including various heterogeneity statistics and the model coefficient(s)).}
   \item{incl}{a matrix indicating which studies were included in which subset.}
   \item{\dots}{some additional elements/values.}

   The results can be printed with the \code{\link[=print.gosh.rma]{print}} function and plotted with the \code{\link[=plot.gosh.rma]{plot}} function.
}
\note{
   On machines with multiple cores, one can try to speed things up by delegating the model fitting to separate worker processes, that is, by setting \code{parallel="snow"} or \code{parallel="multicore"} and \code{ncpus} to some value larger than 1. Parallel processing makes use of the \code{\link[parallel]{parallel}} package, using the \code{\link[parallel]{makePSOCKcluster}} and \code{\link[parallel]{parLapply}} functions when \code{parallel="snow"} or using \code{\link[parallel]{mclapply}} when \code{parallel="multicore"} (the latter only works on Unix/Linux-alikes). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Olkin, I., Dahabreh, I. J., & Trikalinos, T. A. (2012). GOSH - a graphical display of study heterogeneity. \emph{Research Synthesis Methods}, \bold{3}(3), 214--223. \verb{https://doi.org/10.1002/jrsm.1053}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W. (2021). Model checking in meta-analysis. In C. H. Schmid, T. Stijnen, & I. R. White (Eds.), \emph{Handbook of meta-analysis} (pp. 219--254). Boca Raton, FL: CRC Press. \verb{https://doi.org/10.1201/9781315119403}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, and \code{\link{rma.peto}} for functions to fit models for which GOSH plots can be drawn.

   \code{\link[=influence.rma.uni]{influence}} for other model diagnostics.
}
\examples{
### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001)

### meta-analysis of all trials including ISIS-4 using an equal-effects model
res <- rma(yi, vi, data=dat, method="EE")

### fit FE model to all possible subsets (65535 models)
\dontrun{
sav <- gosh(res, progbar=FALSE)
sav

### create GOSH plot
### red points for subsets that include and blue points
### for subsets that exclude study 16 (the ISIS-4 trial)
plot(sav, out=16, breaks=100)
}
}
\keyword{methods}
