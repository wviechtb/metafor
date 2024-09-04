\name{anova.rma}
\alias{anova}
\alias{anova.rma}
\title{Likelihood Ratio and Wald-Type Tests for 'rma' Objects}
\description{
   For two (nested) models of class \code{"rma.uni"} or \code{"rma.mv"}, the function provides a full versus reduced model comparison in terms of model fit statistics and a likelihood ratio test. When a single model is specified, a Wald-type test of one or more model coefficients or linear combinations thereof is carried out. \loadmathjax
}
\usage{
\method{anova}{rma}(object, object2, btt, X, att, Z, rhs, adjust, digits, refit=FALSE, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"} or \code{"rma.mv"}.}
   \item{object2}{an (optional) object of class \code{"rma.uni"} or \code{"rma.mv"}. Only relevant when conducting a model comparison and likelihood ratio test. See \sQuote{Details}.}
   \item{btt}{optional vector of indices (or list thereof) to specify which coefficients should be included in the Wald-type test. Can also be a string to \code{\link{grep}} for. See \sQuote{Details}.}
   \item{X}{optional numeric vector or matrix to specify one or more linear combinations of the coefficients in the model that should be tested. See \sQuote{Details}.}
   \item{att}{optional vector of indices (or list thereof) to specify which scale coefficients should be included in the Wald-type test. Can also be a string to \code{\link{grep}} for. See \sQuote{Details}. Only relevant for location-scale models (see \code{\link{rma.uni}}).}
   \item{Z}{optional numeric vector or matrix to specify one or more linear combinations of the scale coefficients in the model that should be tested. See \sQuote{Details}. Only relevant for location-scale models (see \code{\link{rma.uni}}).}
   \item{rhs}{optional scalar or vector of values for the right-hand side of the null hypothesis when testing a set of coefficients (via \code{btt} or \code{att}) or linear combinations thereof (via \code{X} or \code{Z}). If unspecified, this defaults to a vector of zeros of the appropriate length. See \sQuote{Details}.}
   \item{adjust}{optional argument to specify (as a character string) a method for adjusting the p-values of Wald-type tests for multiple testing. See \code{\link{p.adjust}} for possible options. Can be abbreviated. Can also be a logical and if \code{TRUE}, then a Bonferroni correction is used.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{refit}{logical to specify whether models fitted with REML estimation and differing in their fixed effects should be refitted with ML estimation when conducting a likelihood ratio test (the default is \code{FALSE}).}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used in three different ways:

   \enumerate{

   \item When a single model is specified (via argument \code{object}), the function provides a Wald-type test of one or more model coefficients, that is, \mjdeqn{\text{H}_0{:}\; \beta_{j \in \texttt{btt}} = 0,}{H_0: \beta_{j ∈ btt} = 0,} where \mjeqn{\beta_{j \in \texttt{btt}}}{\beta_{j ∈ btt}} is the set of coefficients to be tested (by default whether the set of coefficients is significantly different from zero, but one can specify a different value under the null hypothesis via argument \code{rhs}).

   In particular, for equal- or random-effects models (i.e., models without moderators), this is just the test of the single coefficient of the model (i.e., \mjeqn{\text{H}_0{:}\; \theta = 0}{H_0: \theta = 0} or \mjeqn{\text{H}_0{:}\; \mu = 0}{H_0: \mu = 0}). For models including moderators, an omnibus test of all model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all coefficients in the model including the first.

   Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} (\sQuote{betas to test}) argument. For example, with \code{btt=c(3,4)}, only the third and fourth coefficients from the model are included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model). Instead of specifying the coefficient numbers, one can specify a string for \code{btt}. In that case, \code{\link{grep}} will be used to search for all coefficient names that match the string (and hence, one can use regular expressions to fine-tune the search for matching strings). Using the \code{btt} argument, one can for example select all coefficients corresponding to a particular factor to test if the factor as a whole is significant. One can also specify a list of indices/strings, in which case tests of all list elements will be conducted. See \sQuote{Examples}.

   For location-scale models fitted with the \code{\link{rma.uni}} function, one can use the \code{att} argument in an analogous manner to specify the indices of the scale coefficients to test (i.e., \mjeqn{\text{H}_0{:}\; \alpha_{j \in \texttt{att}} = 0}{H_0: \alpha_{j ∈ att} = 0}, where \mjeqn{\alpha_{j \in \texttt{att}}}{\alpha_{j ∈ att}} is the set of coefficients to be tested).

   \item When a single model is specified (via argument \code{object}), one can use the \code{X} argument\mjseqn{^1} to specify a linear combination of the coefficients in the model that should be tested using a Wald-type test, that is, \mjdeqn{\text{H}_0{:}\; \tilde{x} \beta = 0,}{H_0: x \beta = 0,} where \mjeqn{\tilde{x}}{x} is a (row) vector of the same length as there are coefficients in the model (by default whether the linear combination is significantly different from zero, but one can specify a different value under the null hypothesis via argument \code{rhs}). One can also specify a matrix of linear combinations via the \code{X} argument to test \mjdeqn{\text{H}_0{:}\; \tilde{X} \beta = 0,}{H_0: X \beta = 0,} where each row of \mjeqn{\tilde{X}}{X} defines a particular linear combination to be tested (if \code{rhs} is used, then it should either be a scalar or of the same length as the number of combinations to be tested). If the matrix is of full rank, an omnibus Wald-type test of all linear combinations is also provided. Linear combinations can also be obtained with the \code{\link[=predict.rma]{predict}} function, which provides corresponding confidence intervals. See also the \code{\link{pairmat}} function for constructing a matrix of pairwise contrasts for testing the levels of a categorical moderator against each other.

   For location-scale models fitted with the \code{\link{rma.uni}} function, one can use the \code{Z} argument in an analogous manner to specify one or multiple linear combinations of the scale coefficients in the model that should be tested (i.e., \mjeqn{\text{H}_0{:}\; \tilde{Z} \alpha = 0}{H_0: Z \alpha = 0}).

   \item When specifying two models for comparison (via arguments \code{object} and \code{object2}), the function provides a likelihood ratio test (LRT) comparing the two models. The two models must be based on the same set of data, must be of the same class, and should be nested for the LRT to make sense. Also, LRTs are not meaningful when using REML estimation and the two models differ in terms of their fixed effects (setting \code{refit=TRUE} automatically refits the two models using ML estimation). Also, the theory underlying LRTs is only really applicable when comparing models that were fitted with ML/REML estimation, so if some other estimation method was used to fit the two models, the results should be treated with caution.

   }

   ---------

   \mjseqn{^1} This argument used to be called \code{L}, but was renamed to \code{X} (but using \code{L} in place of \code{X} still works).

}
\value{
   An object of class \code{"anova.rma"}. When a single model is specified (without any further arguments or together with the \code{btt} or \code{att} argument), the object is a list containing the following components:
   \item{QM}{test statistic of the Wald-type test of the model coefficients.}
   \item{QMdf}{corresponding degrees of freedom.}
   \item{QMp}{corresponding p-value.}
   \item{btt}{indices of the coefficients tested by the Wald-type test.}
   \item{k}{number of outcomes included in the analysis.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the Wald-type test.}
   \item{\dots}{some additional elements/values.}
   When \code{btt} or \code{att} was a list, then the object is a list of class \code{"list.anova.rma"}, where each element is an \code{"anova.rma"} object as described above.

   When argument \code{X} is used, the object is a list containing the following components:
   \item{QM}{test statistic of the omnibus Wald-type test of all linear combinations.}
   \item{QMdf}{corresponding degrees of freedom.}
   \item{QMp}{corresponding p-value.}
   \item{hyp}{description of the linear combinations tested.}
   \item{Xb}{values of the linear combinations.}
   \item{se}{standard errors of the linear combinations.}
   \item{zval}{test statistics of the linear combinations.}
   \item{pval}{corresponding p-values.}

   When two models are specified, the object is a list containing the following components:
   \item{fit.stats.f}{log-likelihood, deviance, AIC, BIC, and AICc for the full model.}
   \item{fit.stats.r}{log-likelihood, deviance, AIC, BIC, and AICc for the reduced model.}
   \item{parms.f}{number of parameters in the full model.}
   \item{parms.r}{number of parameters in the reduced model.}
   \item{LRT}{likelihood ratio test statistic.}
   \item{pval}{corresponding p-value.}
   \item{QE.f}{test statistic of the test for (residual) heterogeneity from the full model.}
   \item{QE.r}{test statistic of the test for (residual) heterogeneity from the reduced model.}
   \item{tau2.f}{estimated \mjseqn{\tau^2} value from the full model. \code{NA} for \code{"rma.mv"} objects.}
   \item{tau2.r}{estimated \mjseqn{\tau^2} value from the reduced model. \code{NA} for \code{"rma.mv"} objects.}
   \item{R2}{amount (in percent) of the heterogeneity in the reduced model that is accounted for in the full model (\code{NA} for \code{"rma.mv"} objects). This can be regarded as a pseudo \mjseqn{R^2} statistic (Raudenbush, 2009). Note that the value may not be very accurate unless \mjseqn{k} is large (Lopez-Lopez et al., 2014).}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{\link[=print.anova.rma]{print}} function. To format the results as a data frame, one can use the \code{\link[=as.data.frame.anova.rma]{as.data.frame}} function.
}
\note{
   The function can also be used to conduct a likelihood ratio test (LRT) for the amount of (residual) heterogeneity in random- and mixed-effects models. The full model should then be fitted with either \code{method="ML"} or \code{method="REML"} and the reduced model with \code{method="EE"} (or with \code{tau2=0}). The p-value for the test is based on a chi-square distribution with 1 degree of freedom, but actually needs to be adjusted for the fact that the parameter (i.e., \mjseqn{\tau^2}) falls on the boundary of the parameter space under the null hypothesis (see Viechtbauer, 2007, for more details).

   LRTs for variance components in more complex models (as fitted with the \code{\link{rma.mv}} function) can also be conducted in this manner (see \sQuote{Examples}).
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to meta-analysis with random effects. \emph{Statistics in Medicine}, \bold{15}(6), 619--629. \verb{https://doi.org/10.1002/(sici)1097-0258(19960330)15:6\%3C619::aid-sim188\%3E3.0.co;2-a}

   Huizenga, H. M., Visser, I., & Dolan, C. V. (2011). Testing overall and moderator effects in random effects meta-regression. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{64}(1), 1--19. \verb{https://doi.org/10.1348/000711010X522687}

   \enc{López-López}{Lopez-Lopez}, J. A., \enc{Marín-Martínez}{Marin-Martinez}, F., \enc{Sánchez-Meca}{Sanchez-Meca}, J., Van den Noortgate, W., & Viechtbauer, W. (2014). Estimation of the predictive power of the model in mixed-effects meta-regression: A simulation study. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{67}(1), 30--48. \verb{https://doi.org/10.1111/bmsp.12002}

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2007). Hypothesis tests for population heterogeneity in meta-analysis. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{60}(1), 29--60. \verb{https://doi.org/10.1348/000711005X64042}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., & \enc{López-López}{Lopez-Lopez}, J. A. (2022). Location-scale models for meta-analysis. \emph{Research Synthesis Methods}. \bold{13}(6), 697--715. \verb{https://doi.org/10.1002/jrsm.1562}
}
\seealso{
   \code{\link{rma.uni}} and \code{\link{rma.mv}} for functions to fit models for which likelihood ratio and Wald-type tests can be conducted.

   \code{\link[=print.anova.rma]{print}} for the print method and \code{\link[=as.data.frame.anova.rma]{as.data.frame}} for the method to format the results as a data frame.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model
res1 <- rma(yi, vi, data=dat, method="ML")
res1

### fit mixed-effects model with two moderators (absolute latitude and publication year)
res2 <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="ML")
res2

### Wald-type test of the two moderators
anova(res2)

### alternative way of specifying the same test
anova(res2, X=rbind(c(0,1,0), c(0,0,1)))

### corresponding likelihood ratio test
anova(res1, res2)

### Wald-type test of a linear combination
anova(res2, X=c(1,35,1970))

### use predict() to obtain the same linear combination (with its CI)
predict(res2, newmods=c(35,1970))

### Wald-type tests of several linear combinations
anova(res2, X=cbind(1,seq(0,60,by=10),1970))

### adjust for multiple testing with the Bonferroni method
anova(res2, X=cbind(1,seq(0,60,by=10),1970), adjust="bonf")

### mixed-effects model with three moderators
res3 <- rma(yi, vi, mods = ~ ablat + year + alloc, data=dat, method="ML")
res3

### Wald-type test of the 'alloc' factor
anova(res3, btt=4:5)

### instead of specifying the coefficient numbers, grep for "alloc"
anova(res3, btt="alloc")

### specify a list for the 'btt' argument
anova(res3, btt=list(2,3,4:5))

### adjust for multiple testing with the Bonferroni method
anova(res3, btt=list(2,3,4:5), adjust="bonf")

############################################################################

### an example of doing LRTs of variance components in more complex models
dat <- dat.konstantopoulos2011
res <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat)

### likelihood ratio test of the district-level variance component
res0 <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat, sigma2=c(0,NA))
anova(res, res0)

### likelihood ratio test of the school-level variance component
res0 <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat, sigma2=c(NA,0))
anova(res, res0)

### likelihood ratio test of both variance components simultaneously
res0 <- rma.mv(yi, vi, data=dat)
anova(res, res0)

############################################################################

### an example illustrating a workflow involving cluster-robust inference
dat <- dat.assink2016

### assume that the effect sizes within studies are correlated with rho=0.6
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)

### fit multilevel model using this approximate V matrix
res <- rma.mv(yi, V, random = ~ 1 | study/esid, data=dat)
res

### likelihood ratio tests of the two variance components
res0 <- rma.mv(yi, V, random = ~ 1 | study/esid, data=dat, sigma2=c(0,NA))
anova(res, res0)
res0 <- rma.mv(yi, V, random = ~ 1 | study/esid, data=dat, sigma2=c(NA,0))
anova(res, res0)

### use cluster-robust methods for inferences about the fixed effects
sav <- robust(res, cluster=study, clubSandwich=TRUE)
sav

### examine if 'deltype' is a potential moderator
res <- rma.mv(yi, V, mods = ~ deltype, random = ~ 1 | study/esid, data=dat)
sav <- robust(res, cluster=study, clubSandwich=TRUE)
sav

### note: the (denominator) dfs for the omnibus F-test are very low, so the results
### of this test may not be trustworthy; consider using cluster wild bootstrapping
\dontrun{
library(wildmeta)
Wald_test_cwb(res, constraints=constrain_zero(2:3), R=1000, seed=1234)
}
}
\keyword{models}
