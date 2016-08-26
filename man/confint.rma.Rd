\name{confint.rma.uni}
\alias{confint.rma.uni}
\alias{confint.rma.mh}
\alias{confint.rma.peto}
\alias{confint.rma.glmm}
\alias{confint.rma.mv}
\title{Confidence Intervals for 'rma' Objects}
\description{
   The function calculates confidence intervals for the model coefficients and/or for variance/correlation parameters in the model (the latter only for \code{\link{rma.uni}} and \code{\link{rma.mv}} objects).
}
\usage{
\method{confint}{rma.uni}(object, parm, level, fixed=FALSE, random=TRUE,
        digits, transf, targs, verbose=FALSE, control, \dots)

\method{confint}{rma.mh}(object, parm, level, digits, transf, targs, \dots)

\method{confint}{rma.peto}(object, parm, level, digits, transf, targs, \dots)

\method{confint}{rma.glmm}(object, parm, level, digits, transf, targs, \dots)

\method{confint}{rma.mv}(object, parm, level, fixed=FALSE, sigma2, tau2, rho, gamma2, phi,
        digits, transf, targs, verbose=FALSE, control, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, or \code{"rma.mv"}. The method is not yet implemented for objects of class \code{"rma.glmm"}.}
   \item{parm}{this argument is here for compatability with the generic function \code{\link{confint}}, but is (currently) ignored.}
   \item{fixed}{logical indicating whether confidence intervals for the model coefficients should be returned (the default is \code{FALSE}).}
   \item{random}{logical indicating whether a confidence interval for the amount of (residual) heterogeneity should be returned (the default is \code{TRUE}) (only relevant for \code{"rma.uni"} objects).}
   \item{sigma2}{integer specifying for which \if{latex}{\eqn{\sigma^2}}\if{html}{\eqn{\sigma²}}\if{text}{\eqn{\sigma^2}} value a confidence interval should be obtained (only relevant for \code{"rma.mv"} objects).}
   \item{tau2}{integer specifying for which \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value a confidence interval should be obtained (only relevant for \code{"rma.mv"} objects).}
   \item{rho}{integer specifying for which \eqn{\rho} value the confidence interval should be obtained (only relevant for \code{"rma.mv"} objects).}
   \item{gamma2}{integer specifying for which \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}} value a confidence interval should be obtained (only relevant for \code{"rma.mv"} objects).}
   \item{phi}{integer specifying for which \eqn{\phi} value a confidence interval should be obtained (only relevant for \code{"rma.mv"} objects).}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (if unspecified, the default is to take the value from the object).}
   \item{digits}{integer specifying the number of decimal places to which the results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{transf}{optional argument specifying the name of a function that should be used to transform the model coefficients and interval bounds (e.g., \code{transf=exp}; see also \link{transf}). If unspecified, no transformation is used.}
   \item{targs}{optional arguments needed by the function specified under \code{transf}.}
   \item{verbose}{logical indicating whether output should be generated on the progress of the iterative algorithms used to obtain the confidence intervals (the default is \code{FALSE}). See \sQuote{Note}.}
   \item{control}{list of control values for the iterative algorithm. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{other arguments.}
}
\details{
   The confidence intervals for the model coefficients are simply the usual Wald-type intervals which are also shown when printing the fitted object.

   For objects of class \code{"rma.uni"} obtained with the \code{\link{rma.uni}} function, the confidence interval for the amount of (residual) heterogeneity (i.e., \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}}) is obtained iteratively either via the Q-profile method or via the generalized Q-statistic method (Hartung and Knapp, 2005; Viechtbauer, 2007; Jackson, 2013; Jackson et al., 2014). The latter is automatically used when the model was fitted with \code{method="GENQ"}, the former is used in all other cases. Either method provides an exact confidence interval for \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} in random- and mixed-effects models. The square root of the interval bounds is also returned for easier interpretation. Confidence intervals for \if{latex}{\eqn{I^2}}\if{html}{\eqn{I²}}\if{text}{\eqn{I^2}} and \if{latex}{\eqn{H^2}}\if{html}{\eqn{H²}}\if{text}{\eqn{H^2}} are also provided (Higgins & Thompson, 2002). Since \if{latex}{\eqn{I^2}}\if{html}{\eqn{I²}}\if{text}{\eqn{I^2}} and \if{latex}{\eqn{H^2}}\if{html}{\eqn{H²}}\if{text}{\eqn{H^2}} are just monotonic transformations of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} (for details, see \code{\link{print.rma.uni}}), the confidence intervals for \if{latex}{\eqn{I^2}}\if{html}{\eqn{I²}}\if{text}{\eqn{I^2}} and \if{latex}{\eqn{H^2}}\if{html}{\eqn{H²}}\if{text}{\eqn{H^2}} are also exact.

   Usually, the estimate of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} from the random/mixed-effects model will fall within the confidence interval provided by the Q-profile method. However, this is not guaranteed. Depending on the method used to estimate \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} and the width of the confidence interval, it can happen that the confidence interval does not actually contain the estimate (trying to explain this to reviewers can be tricky). However, using the empirical Bayes or Paule-Mandel estimator of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} when fitting the model (i.e., using \code{method="EB"} or \code{method="PM"}) guarantees that the estimate of \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} falls within the confidence interval. When \code{method="GENQ"} was used to fit the model, the corresponding CI obtained via the generalized Q-statistic method is also guaranteed to contain the estimate \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}}.

   For objects of class \code{"rma.mv"} obtained with the \code{\link{rma.mv}} function, confidence intervals are obtained by default for all (non-fixed) variance and correlation components of the model. Alternatively, one can use the \code{sigma2}, \code{tau2}, \code{rho}, \code{gamma2}, or \code{phi} arguments to specify for which variance/correlation parameter a confidence interval should be obtained. Only one of these arguments can be used at a time. A single integer is used to specify the number of the parameter. The function then provides a profile likelihood confidence interval for the specified parameter. It is a good idea to examine the corresponding profile likelihood plot (via the \code{\link{profile.rma.mv}} function) to make sure that the bounds obtained are sensible.

   The methods used to find a confidence interval for a variance/correlation component are iterative and require the use of the \code{\link{uniroot}} function. By default, the desired accuracy is set equal to \code{.Machine$double.eps^0.25} and the maximum number of iterations to \code{1000}. The desired accuracy (\code{tol}) and maximum number of iterations (\code{maxiter}) can be adjusted with \code{control=list(tol=value, maxiter=value)}, but the defaults should be adequate for most purposes. If \code{verbose=TRUE}, output is generated on the progress of the iterative algorithm. This is especially useful when model fitting is slow, in which case finding the confidence interval bounds can also take considerable amounts of time.

   More problematic when using the \code{\link{uniroot}} function is setting appropriate end points of the interval to be searched for the confidence interval bounds. The function tries to set some sensible defaults for the end points, but it may happen that the function is only able to determine that a bound is below/above a certain limit (this is indicated in the output accordingly with \code{<} or \code{>} signs). It can also happen that the model cannot be fitted or does not converge especially at the extremes of the interval to be searched. This will result in missing (\code{NA}) bounds and corresponding warnings. It may then be necessary to adjust the end points manually (see \sQuote{Note}).

   Finally, it is also possible that the lower and upper confidence interval bounds for a variance component both fall below zero. Since both bounds then fall outside of the parameter space, the confidence interval then consists of the null/empty set. Alternatively, one could interpret this as a CI with bounds \eqn{[0,0]} or as indicating 'highly (or overly) homogeneous' data.
}
\value{
   An object of class \code{"confint.rma"}. The object is a list with either one or two elements (named \code{fixed} and \code{random}) with the following elements:
   \item{estimate}{estimate of the model coefficient or variance component.}
   \item{ci.lb}{lower bound of the confidence interval.}
   \item{ci.ub}{upper bound of the confidence interval.}

   The results are formated and printed with the \code{\link{print.confint.rma}} function.
}
\note{
   For objects of class \code{"rma.uni"} with \code{random=TRUE}, the function provides a confidence interval for \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}}. The upper end point of the interval to be searched for the confidence interval bounds is set to 100 (which should be large enough for most cases). The upper end point can be adjusted with \code{control=list(tau2.max=value)}. One can also adjust the lower end point of the interval to be searched with \code{control=list(tau2.min=value)} (the default is to take that value from the object, which is 0 by default). You should only play around with this value if you know what you are doing.

   For objects of class \code{"rma.mv"}, the function provides confidence intervals for the variance/correlation parameters in the model. For variance components, the lower end point of the interval to be searched is set to 0 and the upper end point to the larger of 10 and 100 times the value of the component. For correlations, the function tries to set the lower end point to a sensible default depending on the type of variance structure chosen, while the upper end point is set to 1. One can adjust the lower and/or upper end points with \code{control=list(vc.min=value, vc.max=value)}. Also, the function tries to adjust the lower/upper end points when the model does not converge at these extremes (the end points are then moved closer to the estimated value of the component). The total number of tries for setting/adjusting the end points is determined via \code{control=list(eptries=<value>)}, with the default being 10 tries.

   The Q-profile and generalized Q-statistic methods are both exact under the assumptions of the random- and mixed-effects models (i.e., normally distributed observed and true outcomes and known sampling variances). In practice, these assumptions are usually only approximately true, turning confidence intervals for \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} also into approximations. Profile likelihood confidence intervals are not exact by construction and rely on the asymptotic behavior of the likelihood ratio statistic, so they may be inaccurate in small samples.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Hartung, J., & Knapp, G. (2005). On confidence intervals for the among-group variance in the one-way random effects model with unequal error variances. \emph{Journal of Statistical Planning and Inference}, \bold{127}, 157--177.

   Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine}, \bold{21}, 1539--1558.

   Jackson, D. (2013). Confidence intervals for the between-study variance in random effects meta-analysis using generalised Cochran heterogeneity statistics. \emph{Research Synthesis Methods}, \bold{4}, 220--229.

   Jackson, D., Turner, R., Rhodes, K., & Viechtbauer, W. (2014). Methods for calculating confidence and credible intervals for the residual between-study variance in random effects meta-regression models. \emph{BMC Medical Research Methodology}, \bold{14}, 103.

   Mantel, N., & Haenszel, W. (1959). Statistical aspects of the analysis of data from retrospective studies of disease. \emph{Journal of the National Cancer Institute}, \bold{22}, 719--748.

   Viechtbauer, W. (2007). Confidence intervals for the amount of heterogeneity in meta-analysis. \emph{Statistics in Medicine}, \bold{26}, 37--52.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.

   Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}, 335--371.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, \code{\link{rma.mv}}
}
\examples{
### meta-analysis of the log relative risks using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### confidence interval for the total amount of heterogeneity
confint(res)

### mixed-effects model with absolute latitude in the model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
           mods = ~ ablat, data=dat.bcg)

### confidence interval for the residual amount of heterogeneity
confint(res)

### multilevel random-effects model
res <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat.konstantopoulos2011)

### profile plots and confidence intervals for the variance components
\dontrun{
par(mfrow=c(2,1))
profile(res, sigma2=1)
abline(h=logLik(res) - qchisq(.95, df=1)/2, lty="dotted")
sav <- confint(res, sigma2=1)
sav
abline(v=sav$random[1,2:3], lty="dotted")
profile(res, sigma2=2)
abline(h=logLik(res) - qchisq(.95, df=1)/2, lty="dotted")
sav <- confint(res, sigma2=2)
sav
abline(v=sav$random[1,2:3], lty="dotted")}

### multivariate parameterization of the model
res <- rma.mv(yi, vi, random = ~ factor(school) | district, data=dat.konstantopoulos2011)

### profile plots and confidence intervals for the variance component and correlation
\dontrun{
par(mfrow=c(2,1))
profile(res, tau2=1)
abline(h=logLik(res) - qchisq(.95, df=1)/2, lty="dotted")
sav <- confint(res, tau2=1)
sav
abline(v=sav$random[1,2:3], lty="dotted")
profile(res, rho=1)
abline(h=logLik(res) - qchisq(.95, df=1)/2, lty="dotted")
sav <- confint(res, rho=1)
sav
abline(v=sav$random[1,2:3], lty="dotted")}
}
\keyword{models}
