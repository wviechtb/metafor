\name{hc}
\alias{hc}
\alias{hc.rma.uni}
\title{Meta-Analysis based on the Method by Henmi and Copas (2010)}
\description{
   The function can be used to obtain an estimate of the average true effect and corresponding confidence interval under a random-effects model using the method described by Henmi and Copas (2010).
}
\usage{
hc(object, \dots)

\method{hc}{rma.uni}(object, digits, transf, targs, control, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"}.}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the estimated true effect and the corresponding interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{control}{list of control values for the iterative algorithm. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{other arguments.}
}
\details{
   The model specified via \code{object} must be a model without moderators (i.e., either a fixed- or a random-effects model and not a fixed-effects with moderators or mixed-effects model).

   When using the usual method for fitting a random-effects model (i.e., weighted estimation with inverse-variance weights), the weights assigned to smaller and larger studies become more uniform as the amount of heterogeneity increases. As a consequence, the estimated average effect could become increasingly biased under certain forms of publication bias (where smaller studies on one side of the funnel plot are missing). The method by Henmi and Copas (2010) tries to counteract this problem by providing an estimate of the average true effect that is based on inverse-variance weights as used under a fixed-effects model (which do not take the amount of heterogeneity into consideration). The amount of heterogeneity is still estimated (with the DerSimonian-Laird estimator) and incorporated into the standard error of the estimated average effect and the corresponding confidence interval.

   Currently, there is only a method for handling objects of class \code{"rma.uni"} with the \code{hc} function. It therefore provides a method for conducting a sensitivity analysis after the model has been fitted with the \code{\link{rma.uni}} function.
}
\value{
   An object of class \code{"hc.rma.uni"}. The object is a list containing the following components:
   \item{beta}{estimated average true effect.}
   \item{se}{corresponding standard error.}
   \item{ci.lb}{lower bound of the confidence intervals for the average true effect.}
   \item{ci.ub}{upper bound of the confidence intervals for the average true effect.}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{\link{print.hc.rma.uni}} function.
}
\note{
   The method makes use of the \code{\link{uniroot}} function. By default, the desired accuracy is set equal to \code{.Machine$double.eps^0.25} and the maximum number of iterations to \code{1000}. The desired accuracy (\code{tol}) and the maximum number of iterations (\code{maxiter}) can be adjusted with the \code{control} argument (i.e., \code{control=list(tol=value, maxiter=value)}).
}
\author{
   Original code by Henmi and Copas (2010). Corrected for typos by Michael Dewey (\email{lists@dewey.myzen.co.uk}). Incorporated into the package with some small adjustments for consistency with the other functions in the package by Wolfgang Viechtbauer (\email{wvb@metafor-project.org}).
}
\references{
   Henmi, M., & Copas, J. B. (2010). Confidence intervals for random effects meta-analysis and robustness to publication bias. \emph{Statistics in Medicine}, \bold{29}, 2969--2983.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{https://www.jstatsoft.org/v036/i03}.
}
\seealso{
   \code{\link{rma.uni}}
}
\examples{
### meta-analysis based on log odds ratios
res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.lee2004)
res

### funnel plot as in Henmi and Copas (2010)
funnel(res, yaxis="seinv", refline=0, xlim=c(-3,3), ylim=c(.5,3.5), steps=7, digits=1)

### use method by Henmi and Copas (2010) as a sensitivity analysis
hc(res)

### back-transform results to odds ratio scale
hc(res, transf=exp)
}
\keyword{htest}
