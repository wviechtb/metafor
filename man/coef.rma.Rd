\name{coef.rma}
\alias{coef}
\alias{coef.rma}
\alias{coef.summary.rma}
\title{Extract the Model Coefficients and Coefficient Table from 'rma' and 'summary.rma' Objects}
\description{
   The \code{coef} function extracts the estimated model coefficients from objects of class \code{"rma"}. For objects of class \code{"summary.rma"}, the model coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds are extracted.
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

   When the model was fitted with the Knapp and Hartung (2003) method (i.e., \code{test="knha"} in the \code{\link{rma.uni}} function) or with \code{test="t"} in the \code{\link{rma.glmm}} and \code{\link{rma.mv}} functions, then \code{zval} is called \code{tval} in the data frame that is returned by the function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}
}
\examples{
### meta-analysis of the log risk ratios using a mixed-effects model
### with two moderators (absolute latitude and publication year)
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat + year, data=dat.bcg)
coef(res)
coef(summary(res))
}
\keyword{models}
