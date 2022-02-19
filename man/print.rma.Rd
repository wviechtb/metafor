\name{print.rma}
\alias{print.rma}
\alias{print.rma.uni}
\alias{print.rma.mh}
\alias{print.rma.peto}
\alias{print.rma.glmm}
\alias{print.rma.mv}
\alias{summary}
\alias{summary.rma}
\alias{print.summary.rma}
\title{Print and Summary Methods for 'rma' Objects}
\description{
   Print and summary methods for objects of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, \code{"rma.glmm"}, \code{"rma.glmm"}, and \code{"rma.mv"}. \loadmathjax
}
\usage{
\method{print}{rma.uni}(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)

\method{print}{rma.mh}(x, digits, showfit=FALSE, \dots)

\method{print}{rma.peto}(x, digits, showfit=FALSE, \dots)

\method{print}{rma.glmm}(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)

\method{print}{rma.mv}(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)

\method{summary}{rma}(object, digits, showfit=TRUE, \dots)

\method{print}{summary.rma}(x, digits, showfit=TRUE, signif.stars=getOption("show.signif.stars"),
      signif.legend=signif.stars, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}, \code{"rma.mh"}, \code{"rma.peto"}, \code{"rma.glmm"}, \code{"rma.mv"}, or \code{"summary.rma"} (for \code{print}).}
   \item{object}{an object of class \code{"rma"} (for \code{summary}).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object. See also \link[=misc-options]{here} for further details on how to control the number of digits in the output.}
   \item{showfit}{logical to specify whether the fit statistics and information criteria should be printed (the default is \code{FALSE} for \code{print} and \code{TRUE} for \code{summary}).}
   \item{signif.stars}{logical to specify whether p-values should be encoded visually with \sQuote{significance stars}. Defaults to the \code{show.signif.stars} slot of \code{\link{options}}.}
   \item{signif.legend}{logical to specify whether the legend for the \sQuote{significance stars} should be printed. Defaults to the value for \code{signif.stars}.}
   \item{\dots}{other arguments.}
}
\details{
   The output includes:
   \itemize{
   \item the log-likelihood, deviance, AIC, BIC, and AICc value (when setting \code{showfit=TRUE} or by default for \code{summary}).
   \item for objects of class \code{"rma.uni"} and \code{"rma.glmm"}, the amount of (residual) heterogeneity in the random/mixed-effects model (i.e., the estimate of \mjseqn{\tau^2} and its square root). Suppressed for equal-effects models. The (asymptotic) standard error of the estimate of \mjseqn{\tau^2} is also provided (where possible).
   \item for objects of \code{"rma.mv"}, a table providing information about the variance components and correlations in the model. For \mjseqn{\sigma^2} components, the estimate and its square root are provided, in addition to the number of values/levels, whether the component was fixed or estimated, and the name of the grouping variable/factor. If the \code{R} argument was used to specify known correlation matrices, this is also indicated. For models with an \sQuote{\code{~ inner | outer}} formula term, the name of the inner and outer grouping variable/factor are given and the number of values/levels of these variables/factors. In addition, for each \mjseqn{\tau^2} component, the estimate and its square root are provided, the number of effects or outcomes observed at each level of the inner grouping variable/factor (only for \code{struct="HCS"}, \code{struct="DIAG"}, \code{struct="HAR"}, and \code{struct="UN"}), and whether the component was fixed or estimated. Finally, either the estimate of \mjseqn{\rho} (for \code{struct="CS"}, \code{struct="AR"}, \code{struct="CAR"}, \code{struct="HAR"}, or \code{struct="HCS"}) or the entire estimated correlation matrix (for \code{struct="UN"}) between the levels of the inner grouping variable/factor is provided, again with information whether a particular correlation was fixed or estimated, and how often each combination of levels of the inner grouping variable/factor was observed across the levels of the outer grouping variable/factor. If there is a second \sQuote{\code{~ inner | outer}} formula term, the same information as described above will be provided, but now for the \mjseqn{\gamma^2} and \mjseqn{\phi} components.
   \item the \mjseqn{I^2} statistic, which estimates (in percent) how much of the total variability in the observed effect sizes or outcomes (which is composed of heterogeneity plus sampling variability) can be attributed to heterogeneity among the true effects. For a meta-regression model, \mjseqn{I^2} estimates how much of the unaccounted variability (which is composed of residual heterogeneity plus sampling variability) can be attributed to residual heterogeneity. See \sQuote{Note} for how \mjseqn{I^2} is computed.
   \item the \mjseqn{H^2} statistic, which estimates the ratio of the total amount of variability in the observed effect sizes or outcomes to the amount of sampling variability. For a meta-regression model, \mjseqn{H^2} estimates the ratio of the unaccounted variability in the observed effect sizes or outcomes to the amount of sampling variability. See \sQuote{Note} for how \mjseqn{H^2} is computed.
   \item for objects of class \code{"rma.uni"}, the \mjseqn{R^2} statistic, which estimates the amount of heterogeneity accounted for by the moderators included in the model and can be regarded as a pseudo \mjseqn{R^2} statistic (Raudenbush, 2009). Only provided when fitting a model including moderators. This is suppressed (and set to \code{NULL}) for models without moderators or if the model does not contain an intercept. See \sQuote{Note} for how \mjseqn{R^2} is computed.
   \item for objects of class \code{"rma.glmm"}, the amount of study level variability (only when using a model that models study level differences as a random effect).
   \item the results of the test for (residual) heterogeneity. This is the usual \mjseqn{Q}-test for heterogeneity when not including moderators in the model and the \mjseqn{Q_E}-test for residual heterogeneity when moderators are included. For objects of class \code{"rma.glmm"}, the results from a Wald-type test and a likelihood ratio test are provided (see \code{\link{rma.glmm}} for more details).
   \item the results of the omnibus (Wald-type) test of the coefficients in the model (the indices of the coefficients tested are also indicated). Suppressed if the model includes only one coefficient (e.g., only an intercept, like in the equal- and random-effects models).
   \item a table with the estimated coefficients, corresponding standard errors, test statistics, p-values, and confidence interval bounds.
   \item the Cochran-Mantel-Haenszel test and Tarone's test for heterogeneity (only when analyzing odds ratios using the Mantel-Haenszel method, i.e., \code{"rma.mh"}).
   }
   See also \link[=misc-options]{here} for details on the option to create styled/colored output with the help of the \href{https://cran.r-project.org/package=crayon}{crayon} package.
}
\value{
   The \code{print} functions do not return an object. The \code{summary} function returns the object passed to it (with additional class \code{"summary.rma"}).
}
\note{
   For random-effects models, the \mjseqn{I^2} statistic is computed with \mjdeqn{I^2 = 100\\\\\\\% \times \frac{\hat{\tau}^2}{\hat{\tau}^2 + \tilde{v}},}{I^2 = 100\\\% hat(\tau)^2 / (hat(\tau)^2 + v),} where \mjeqn{\hat{\tau}^2}{hat(\tau)^2} is the estimated value of \mjseqn{\tau^2} and \mjdeqn{\tilde{v} = \frac{(k-1) \sum w_i}{(\sum w_i)^2 - \sum w_i^2},}{v = ((k-1) \sum w_i) / ((\sum w_i)^2 - \sum w_i^2),} where \mjseqn{w_i = 1 / v_i} is the inverse of the sampling variance of the \mjeqn{i\textrm{th}}{ith} study (\mjeqn{\tilde{v}}{v} is equation 9 in Higgins & Thompson, 2002, and can be regarded as the \sQuote{typical} within-study variance of the observed effect sizes or outcomes). The \mjseqn{H^2} statistic is computed with \mjdeqn{H^2 = \frac{\hat{\tau}^2 + \tilde{v}}{\tilde{v}}.}{H^2 = (hat(\tau)^2 + v) / v.} Analogous equations are used for mixed-effects models.

   Therefore, depending on the estimator of \mjseqn{\tau^2} used, the values of \mjseqn{I^2} and \mjseqn{H^2} will change. For random-effects models, \mjseqn{I^2} and \mjseqn{H^2} are often computed with \mjseqn{I^2 = (Q-(k-1))/Q} and \mjseqn{H^2 = Q/(k-1)}, where \mjseqn{Q} denotes the statistic of the test for heterogeneity and \mjseqn{k} the number of studies (i.e., observed effect sizes or outcomes) included in the meta-analysis. The equations used in the metafor package to compute these statistics are more general and have the advantage that the values of \mjseqn{I^2} and \mjseqn{H^2} will be consistent with the estimated value of \mjseqn{\tau^2} (i.e., if \mjeqn{\hat{\tau}^2 = 0}{hat(\tau)^2 = 0}, then \mjseqn{I^2 = 0} and \mjseqn{H^2 = 1} and if \mjeqn{\hat{\tau}^2 > 0}{hat(\tau)^2 > 0}, then \mjseqn{I^2 > 0} and \mjseqn{H^2 > 1}).

   The two definitions of \mjseqn{I^2} and \mjseqn{H^2} actually coincide when using the DerSimonian-Laird estimator of \mjseqn{\tau^2} (i.e., the commonly used equations are actually special cases of the more general definitions given above). Therefore, if you prefer the more conventional definitions of these statistics, use \code{method="DL"} when fitting the random/mixed-effects model with the \code{\link{rma.uni}} function. The conventional definitions are also automatically used when fitting an equal-effects models.

   For mixed-effects models, the pseudo \mjseqn{R^2} statistic (Raudenbush, 2009) is computed with \mjdeqn{R^2 = \frac{\hat{\tau}_{RE}^2 - \hat{\tau}_{ME}^2}{\hat{\tau}_{RE}^2},}{R^2 = (hat(\tau)^2_RE - hat(\tau)^2_ME) / hat(\tau)^2_RE,} where \mjeqn{\hat{\tau}_{RE}^2}{hat(\tau)^2_RE} denotes the estimated value of \mjseqn{\tau^2} based on the random-effects model (i.e., the total amount of heterogeneity) and \mjeqn{\hat{\tau}_{ME}^2}{hat(\tau)^2_ME} denotes the estimated value of \mjseqn{\tau^2} based on the mixed-effects model (i.e., the residual amount of heterogeneity). It can happen that \mjeqn{\hat{\tau}_{RE}^2 < \hat{\tau}_{ME}^2}{hat(\tau)^2_RE < hat(\tau)^2_ME}, in which case \mjseqn{R^2} is set to zero (and also if \mjeqn{\hat{\tau}_{RE}^2 = 0}{hat(\tau)^2_RE = 0}). Again, the value of \mjseqn{R^2} will change depending on the estimator of \mjseqn{\tau^2} used. This statistic is only computed when the mixed-effects model includes an intercept (so that the random-effects model is clearly nested within the mixed-effects model). You can also use the \code{\link[=anova.rma]{anova}} function to compute \mjseqn{R^2} for any two models that are known to be nested. Note that the pseudo \mjseqn{R^2} statistic may not be very accurate unless \mjseqn{k} is large (Lopez-Lopez et al., 2014).

   For fixed-effects with moderators models, the \mjseqn{R^2} statistic is simply the standard \mjseqn{R^2} statistic (also known as the \sQuote{coefficient of determination}) computed based on weighted least squares estimation. To be precise, the so-called \sQuote{adjusted} \mjseqn{R^2} statistic is provided, since \mjseqn{k} is often relatively small in meta-analyses, in which case the adjustment is relevant.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity in a meta-analysis. \emph{Statistics in Medicine}, \bold{21}(11), 1539--1558. \verb{https://doi.org/10.1002/sim.1186}

   \enc{López-López}{Lopez-Lopez}, J. A., \enc{Marín-Martínez}{Marin-Martinez}, F., \enc{Sánchez-Meca}{Sanchez-Meca}, J., Van den Noortgate, W., & Viechtbauer, W. (2014). Estimation of the predictive power of the model in mixed-effects meta-regression: A simulation study. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{67}(1), 30--48. \verb{https://doi.org/10.1111/bmsp.12002}

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for the corresponding model fitting functions.
}
\keyword{print}
