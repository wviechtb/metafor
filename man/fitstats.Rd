\name{fitstats}
\alias{fitstats}
\alias{fitstats.rma}
\alias{logLik.rma}
\alias{deviance.rma}
\alias{AIC.rma}
\alias{BIC.rma}
\alias{nobs.rma}
\alias{df.residual.rma}
\title{Fit Statistics and Information Criteria for 'rma' Objects}
\description{
   Functions to extract the log-likelihood, deviance, AIC, BIC, and AICc values from objects of class \code{"rma"}. \loadmathjax
}
\usage{
fitstats(object, \dots)

\method{fitstats}{rma}(object, \dots, REML)

\method{logLik}{rma}(object, REML, \dots)
\method{deviance}{rma}(object, REML, \dots)

\method{AIC}{rma}(object, \dots, k=2, correct=FALSE)
\method{BIC}{rma}(object, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{\dots}{optionally more fitted model objects (only for \code{fitstats()}, \code{AIC()}, and \code{BIC()}).}
   \item{REML}{logical to specify whether the regular or restricted likelihood function should be used to obtain the fit statistics and information criteria. Defaults to the method of estimation used (i.e., \code{TRUE} if \code{object} was fitted with \code{method="REML"} and \code{FALSE} otherwise).}
   \item{k}{numeric value to specify the penalty per parameter. The default (\code{k=2}) is the classical AIC. See \code{\link{AIC}} for more details.}
   \item{correct}{logical to specify whether the regular (default) or corrected (i.e., AICc) should be extracted.}
}
\value{
   For \code{fitstats}, a data frame with the (restricted) log-likelihood, deviance, AIC, BIC, and AICc values for each model passed to the function.

   For \code{logLik}, an object of class \code{"logLik"}, providing the (restricted) log-likelihood of the model evaluated at the estimated coefficient(s).

   For \code{deviance}, a numeric value with the corresponding deviance.

   For \code{AIC} and \code{BIC}, either a numeric value with the corresponding AIC, AICc, or BIC or a data frame with rows corresponding to the models and columns representing the number of parameters in the model (\code{df}) and the AIC, AICc, or BIC.
}
\note{
   Variance components in the model (e.g., \mjseqn{\tau^2} in random/mixed-effects models fitted with \code{\link{rma.uni}}) are counted as additional parameters in the calculation of the AIC, BIC, and AICc. Also, the fixed effects are counted as parameters in the calculation of the AIC, BIC, and AICc even when using REML estimation.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which fit statistics and information criteria can be extracted.

   \code{\link{anova.rma}} for a function to conduct likelihood ratio tests.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res1 <- rma(yi, vi, data=dat, method="ML")

### mixed-effects model with absolute latitude and publication year as moderators
res2 <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="ML")

### compare fit statistics
fitstats(res1, res2)

### log-likelihoods
logLik(res1)
logLik(res2)

### deviances
deviance(res1)
deviance(res2)

### AIC, AICc, and BIC values
AIC(res1, res2)
AIC(res1, res2, correct=TRUE)
BIC(res1, res2)
}
\keyword{models}
