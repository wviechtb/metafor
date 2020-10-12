\name{trimfill}
\alias{trimfill}
\alias{trimfill.rma.uni}
\title{Trim and Fill Analysis for 'rma.uni' Objects}
\description{
   Carry out a trim and fill analysis for objects of class \code{"rma.uni"}. \loadmathjax
}
\usage{
trimfill(x, \dots)

\method{trimfill}{rma.uni}(x, side, estimator="L0", maxiter=100, verbose=FALSE, ilim, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}.}
   \item{side}{either \code{"left"} or \code{"right"}, indicating on which side of the funnel plot the missing studies should be imputed. If left undefined, the side is chosen within the function depending on the results of Egger's regression test (see \code{\link{regtest}} for details on this test).}
   \item{estimator}{either \code{"L0"}, \code{"R0"}, or \code{"Q0"}, indicating the estimator to use for estimating the number of missing studies (the default is \code{"L0"}).}
   \item{maxiter}{integer indicating the maximum number of iterations to use for the trim and fill method (the default is \code{100}).}
   \item{verbose}{logical indicating whether output should be generated on the progress of the iterative algorithm used as part of the trim and fill method (the default is \code{FALSE}).}
   \item{ilim}{limits for the imputed values. If unspecified, no limits are used.}
   \item{\dots}{other arguments.}
}
\details{
   The trim and fill method is a nonparametric (rank-based) data augmentation technique proposed by Duval and Tweedie (2000a, 2000b; see also Duval, 2005). The method can be used to estimate the number of studies missing from a meta-analysis due to suppression of the most extreme results on one side of the funnel plot. The method then augments the observed data so that the funnel plot is more symmetric and recomputes the summary estimate based on the complete data. The trim and fill method can only be used in the context of a fixed- or random-effects model (i.e., in models without moderators). The method should not be regarded as a way of yielding a more \sQuote{valid} estimate of the overall effect or outcome, but as a way of examining the sensitivity of the results to one particular selection mechanism (i.e., one particular form of publication bias).
}
\value{
   An object of class \code{c("rma.uni.trimfill","rma.uni","rma")}. The object is a list containing the same components as objects created by \code{\link{rma.uni}}, except that the data are augmented by the trim and fill method. The following components are also added:
   \item{k0}{estimated number of missing studies.}
   \item{side}{either \code{"left"} or \code{"right"}, indicating on which side of the funnel plot the missing studies (if any) were imputed.}
   \item{se.k0}{standard error of k0.}
   \item{p.k0}{p-value for the test of \mjeqn{\mbox{H}_0}{H_0}: no missing studies on the chosen side (only when \code{estimator="R0"}; \code{NA} otherwise).}
   \item{yi}{the observed effect sizes or outcomes plus the augmented values (if there are any).}
   \item{vi}{the corresponding sampling variances}
   \item{fill}{a logical vector indicating which of the values in \code{yi} are the observed (\code{FALSE}) and the augmented (\code{TRUE}) data.}

   The results of the fitted model after the data augmentation are printed with the \code{\link{print.rma.uni}} function. Calling \code{\link{funnel.rma}} on the object provides a funnel plot of the observed and augmented data.
}
\note{
   Three different estimators for the number of missing studies were proposed by Duval and Tweedie (2000a, 2000b). Based on these articles and Duval (2005), \code{"R0"} and \code{"L0"} are recommended. An advantage of estimator \code{"R0"} is that it provides a test of the null hypothesis that the number of missing studies (on the chosen side) is zero.

   If the outcome measure used for the analysis is bounded (e.g., correlations are bounded between \code{-1} and \code{1}, proportions are bounded between \code{0} and \code{1}), one can use the \code{ilim} argument to enforce those limits when imputing values (imputed values cannot exceed those bounds then).

   The model used during the trim and fill procedure is the same as used by the original model object. Hence, if a fixed-effects model is passed to the function, then a fixed-effects model is also used during the trim and fill procedure and the results provided are also based on a fixed-effects model. This would be a \sQuote{fixed-fixed} approach. Similarly, if a random-effects model is passed to the function, then the same model is used as part of the trim and fill procedure and for the final analysis. This would be a \sQuote{random-random} approach. However, one can also easily fit a different model for the final analysis than was used for the trim and fill procedure. See \sQuote{Examples} for an illustration of a \sQuote{fixed-random} approach.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Duval, S. J., & Tweedie, R. L. (2000a). Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. \emph{Biometrics}, \bold{56}(2), 455--463. \verb{https://doi.org/10.1111/j.0006-341x.2000.00455.x}

   Duval, S. J., & Tweedie, R. L. (2000b). A nonparametric "trim and fill" method of accounting for publication bias in meta-analysis. \emph{Journal of the American Statistical Association}, \bold{95}(449), 89--98. \verb{https://doi.org/10.1080/01621459.2000.10473905}

   Duval, S. J. (2005). The trim and fill method. In H. R. Rothstein, A. J. Sutton, & M. Borenstein (Eds.) \emph{Publication bias in meta-analysis: Prevention, assessment, and adjustments} (pp. 127--144). Chichester, England: Wiley.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{funnel.rma}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### meta-analysis of the log risk ratios using a fixed-effects model
res <- rma(yi, vi, data=dat, method="FE")
res.tf <- trimfill(res)
res.tf
funnel(res.tf, legend=TRUE, cex=1.2)

### estimator "R0" also provides test
res.tf <- trimfill(res, estimator="R0")
res.tf

### meta-analysis of the log risk ratios using a random-effects model
res <- rma(yi, vi, data=dat)
res.tf <- trimfill(res)
res.tf
funnel(res.tf, legend=TRUE, cex=1.2)

### the examples above are fixed-fixed and random-random approaches

### illustration of a fixed-random approach
res <- rma(yi, vi, data=dat, method="FE")
res.tf <- trimfill(res)
filled <- data.frame(yi = res.tf$yi, vi = res.tf$vi, fill = res.tf$fill)
filled
rma(yi, vi, data=filled)
}
\keyword{models}
