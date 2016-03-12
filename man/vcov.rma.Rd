\name{vcov.rma}
\alias{vcov.rma}
\title{Extract the Variance-Covariance Matrix of the Parameter Estimates from 'rma' Objects}
\description{
   The function extracts the (estimated) variance-covariance matrix of the (fixed effects) parameter estimates from objects of class \code{"rma"} or the marginal variance-covariance matrix of the observed effect sizes or outcomes.
}
\usage{
\method{vcov}{rma}(object, type="fixed", \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"}.}
   \item{type}{character string indicating whether to return the variance-covariance matrix of the fixed effects (\code{type="fixed"}) or the marginal variance-covariance matrix of the observed effect sizes or outcomes (\code{type="obs"}).}
   \item{\dots}{other arguments.}
}
\details{
   Note that \code{type="obs"} currently only works for object of class \code{"rma.uni"} and \code{"rma.mv"}.

   For objects of class \code{"rma.uni"}, the marginal variance-covariance matrix is just a diagonal matrix with \if{latex}{\eqn{\tau^2 + v_i}}\if{html}{\eqn{\tau² + vᵢ}}\if{text}{\eqn{\tau^2 + v_i}} along the diagonal, where \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} is the estimated amount of (residual) heterogeneity (set to 0 in fixed-effects models) and \if{latex}{\eqn{v_i}}\if{html}{\eqn{vᵢ}}\if{text}{\eqn{v_i}} is the sampling variance of the \eqn{i}th study.

   For objects of class \code{"rma.mv"}, the structure can be more complex and depends on the types of random effects included in the model.
}
\value{
   A matrix corresponding to the requested variance-covariance matrix.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}
}
\examples{
### meta-analysis of the log relative risks using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year, data=dat.bcg)
vcov(res)

### marginal var-cov matrix of the observed log relative risks
vcov(res, type="obs")
}
\keyword{models}
