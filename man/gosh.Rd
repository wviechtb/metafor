\name{gosh}
\alias{gosh}
\alias{gosh.rma}
\title{GOSH Plots for 'rma' Objects}
\description{Function to create GOSH plots for objects of class \code{"rma"}.}
\usage{
gosh(x, \dots)

\method{gosh}{rma}(x, subsets, progbar=TRUE, parallel="no", ncpus=1, cl=NULL, \ldots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{subsets}{optional integer to specify the number of subsets.}
   \item{progbar}{logical indicating whether a progress bar should be shown (the default is \code{TRUE}).}
   \item{parallel}{character string indicating whether parallel processing should be used (the default is \code{"no"}). For parallel processing, set to either \code{"snow"} or \code{"multicore"}. See \sQuote{Details}.}
   \item{ncpus}{integer specifying the number of processes to be used in the parallel operation.}
   \item{cl}{optional snow cluster to use if \code{parallel="snow"}. If not supplied, a cluster on the local machine is created for the duration of the call.}
   \item{\dots}{other arguments.}
}
\details{
   The model specified via \code{x} must be a model fitted with either the \code{\link{rma.uni}}, \code{\link{rma.mh}}, or \code{\link{rma.peto}} function.

   Olkin et al. (2012) proposed the GOSH (graphical display of study heterogeneity) plot, which is based on examining the results of a fixed-effects model in all possible subsets of size \eqn{1, \ldots, k} of the \eqn{k} studies included in a meta-analysis. In a homogeneous set of studies, the model estimates obtained this way should form a roughly symmetric, contiguous, and unimodal distribution. On the other hand, when the distribution is multimodal, then this suggests the presence of heterogeneity, possibly due to outliers and/or distinct subgroupings of studies. Plotting the estimates against some measure of heterogeneity (e.g., \if{latex}{\eqn{I^2}}\if{html}{\eqn{I²}}\if{text}{\eqn{I^2}}, \if{latex}{\eqn{H^2}}\if{html}{\eqn{H²}}\if{text}{\eqn{H^2}}, or the \eqn{Q}-statistic) can also help to reveal subclusters, which are indicative of heterogeneity. The same type of plot can be produced by first fitting a fixed-effects model with either the \code{\link{rma.uni}} (using \code{method="FE"}), \code{\link{rma.mh}}, or \code{\link{rma.peto}} functions and then passing the fitted model object to the \code{gosh} function and then plotting the results.

   For models fitted with the \code{\link{rma.uni}} function (which may involve moderators and/or may be random/mixed-effects models), the idea underlying this type of plot can be generalized by examining the distribution of all model coefficients, plotting them against each other, and against some measure of (residual) heterogeneity (including the estimate of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}}).

   Note that for models without moderators, application of the method requires fitting a total of \eqn{2^k - 1} models, which could be an excessively large number when \eqn{k} is large. For example, for \eqn{k=10}, there are only 1023 possible subsets, but for \eqn{k=20}, this number has already grown to 1048575. For even larger \eqn{k}, it may become computationally infeasible to consider all possible subsets (and even storage of the results could become an issue). Instead, we can then examine (a sufficiently large number of) random subsets.

   By default, if the number of possible subsets is \eqn{\le 10^6}, the function will consider all possible subsets and otherwise \eqn{10^6} random subsets. One can use the \code{subsets} argument to specify a different number of subsets to consider. If \code{subsets} is specified and it is actually larger than the number of possible subsets, then the function automatically only considers the possible subsets and does not use random subsets.

   On machines with multiple cores, one can try to speed things up by delegating the model fitting to separate worker processes, that is, by setting \code{parallel="snow"} or \code{parallel="multicore"} and \code{ncpus} to some value larger than 1. Parallel processing makes use of the \code{\link[parallel]{parallel}} package, using the \code{\link[parallel]{makePSOCKcluster}} and \code{\link[parallel]{parLapply}} functions when \code{parallel="snow"} or using the \code{\link[parallel]{mclapply}} function when \code{parallel="multicore"} (the latter only works on Unix/Linux-alikes). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine. Note that no progress bar will be shown when using parallel processing.
}
\value{
   An object of class \code{"gosh.rma"}. The object is a list containing the following components:
   \item{res}{a data frame with the results for each subset (various heterogeneity statistics and the model coefficient(s)).}
   \item{incl}{a matrix indicating which studies were included in which subset.}
   \item{\dots}{some additional elements/values.}

   The results can be plotted with the \code{\link{plot.gosh.rma}} function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Olkin, I., Dahabreh, I. J., & Trikalinos, T. A. (2012). GOSH - a graphical display of study heterogeneity. \emph{Research Synthesis Methods}, \bold{3}, 214--223.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{influence.rma.uni}},
}
\examples{
### meta-analysis of all trials including ISIS-4 using a FE model
res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, method="FE")

### fit FE model to all possible subsets (65535 models)
sav <- gosh(res)

### create GOSH plot
### red points for subsets that include and blue points
### for subsets that exclude study 16 (the ISIS-4 trial)
plot(sav, out=16, breaks=100, adjust=.3)
}
\keyword{methods}
