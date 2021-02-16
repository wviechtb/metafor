\name{vcov.rma}
\alias{vcov}
\alias{vcov.rma}
\title{Extract Various Types of Variance-Covariance Matrices from 'rma' Objects}
\description{
   The function extracts various types of variance-covariance matrices from objects of class \code{"rma"}. By default, the variance-covariance matrix of the parameter estimates (fixed effects) is returned. \loadmathjax
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

   For objects of class \code{"rma.uni"}, the marginal variance-covariance matrix of the observed effect sizes or outcomes is just a diagonal matrix with \mjeqn{\hat{\tau}^2 + v_i}{\tau^2 + v_i} along the diagonal, where \mjeqn{\hat{\tau}^2}{\tau^2} is the estimated amount of (residual) heterogeneity (set to 0 in fixed-effects models) and \mjseqn{v_i} is the sampling variance of the \mjseqn{i}th study.

   For objects of class \code{"rma.mv"}, the structure can be more complex and depends on the types of random effects included in the model.
}
\value{
   A matrix corresponding to the requested variance-covariance matrix.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}
}
\examples{
### meta-analysis of the log risk ratios using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year, data=dat.bcg)

### var-cov matrix of the fixed effects (i.e., the model coefficients)
vcov(res)

### marginal var-cov matrix of the observed log risk ratios
vcov(res, type="obs")

### var-cov matrix of the fitted values
vcov(res, type="fitted")

### var-cov matrix of the residuals
vcov(res, type="resid")
}
\keyword{models}
