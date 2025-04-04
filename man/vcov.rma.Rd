\name{vcov.rma}
\alias{vcov}
\alias{vcov.rma}
\title{Extract Various Types of Variance-Covariance Matrices from 'rma' Objects}
\description{
   Function to extract various types of variance-covariance matrices from objects of class \code{"rma"}. By default, the variance-covariance matrix of the fixed effects is returned. \loadmathjax
}
\usage{
\method{vcov}{rma}(object, type="fixed", \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{type}{character string to specify the type of variance-covariance matrix to return: \code{type="fixed"} returns the variance-covariance matrix of the fixed effects (the default), \code{type="obs"} returns the marginal variance-covariance matrix of the observed effect sizes or outcomes, \code{type="fitted"} returns the variance-covariance matrix of the fitted values, \code{type="resid"} returns the variance-covariance matrix of the residuals.}
   \item{\dots}{other arguments.}
}
\details{
   Note that \code{type="obs"} currently only works for object of class \code{"rma.uni"} and \code{"rma.mv"}.

   For objects of class \code{"rma.uni"}, the marginal variance-covariance matrix of the observed effect sizes or outcomes is a diagonal matrix with \mjeqn{\hat{\tau}^2 + v_i}{\tau^2 + v_i} along the diagonal, where \mjeqn{\hat{\tau}^2}{\tau^2} is the estimated amount of (residual) heterogeneity (set to 0 in equal-effects models) and \mjseqn{v_i} is the sampling variance of the \mjeqn{i\text{th}}{ith} study.

   For objects of class \code{"rma.mv"}, the structure can be more complex and depends on the random effects included in the model.
}
\value{
   A matrix corresponding to the requested variance-covariance matrix.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which the various types of variance-covariance matrices can be extracted.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### var-cov matrix of the fixed effects (i.e., the model coefficients)
vcov(res)

### marginal var-cov matrix of the observed log risk ratios
round(vcov(res, type="obs"), 3)

### var-cov matrix of the fitted values
round(vcov(res, type="fitted"), 3)

### var-cov matrix of the residuals
round(vcov(res, type="resid"), 3)
}
\keyword{models}
