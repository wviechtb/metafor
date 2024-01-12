\name{coef.rma}
\alias{coef}
\alias{coef.rma}
\alias{coef.summary.rma}
\title{Extract the Model Coefficients and Coefficient Table from 'rma' and 'summary.rma' Objects}
\description{
   Function to extract the estimated model coefficients from objects of class \code{"rma"}. For objects of class \code{"summary.rma"}, the model coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds are extracted.
}
\usage{
\method{coef}{rma}(object, \dots)
\method{coef}{summary.rma}(object, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma"} or \code{"summary.rma"}.}
   \item{\dots}{other arguments.}
}
\value{
   Either a vector with the estimated model coefficient(s) or a data frame with the following elements:
   \item{estimate}{estimated model coefficient(s).}
   \item{se}{corresponding standard error(s).}
   \item{zval}{corresponding test statistic(s).}
   \item{pval}{corresponding p-value(s).}
   \item{ci.lb}{corresponding lower bound of the confidence interval(s).}
   \item{ci.ub}{corresponding upper bound of the confidence interval(s).}

   When the model was fitted with \code{test="t"}, \code{test="knha"}, \code{test="hksj"}, or \code{test="adhoc"}, then \code{zval} is called \code{tval} in the data frame that is returned by the function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for functions to fit models for which model coefficients/tables can be extracted.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit mixed-effects model with absolute latitude and publication year as moderators
res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

### extract model coefficients
coef(res)

### extract model coefficient table
coef(summary(res))
}
\keyword{models}
