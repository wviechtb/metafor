\name{anova.rma}
\alias{anova}
\alias{anova.rma}
\title{Likelihood Ratio and Wald-Type Tests for 'rma' Objects}
\description{
   For two (nested) models of class \code{"rma.uni"} or \code{"rma.mv"}, the function provides a full versus reduced model comparison in terms of model fit statistics and a likelihood ratio test. When a single model is specified, a Wald-type test of one or more model coefficients or linear combinations thereof is carried out. \loadmathjax
}
\usage{
\method{anova}{rma}(object, object2, btt, X, att, Z, digits, \dots)
}
\arguments{
   \item{object}{an object of class \code{"rma.uni"} or \code{"rma.mv"}.}
   \item{object2}{an (optional) object of class \code{"rma.uni"} or \code{"rma.mv"}. See \sQuote{Details}.}
   \item{btt}{optional vector of indices to specify which coefficients should be included in the Wald-type test. Can also be a string to grep for. See \sQuote{Details}.}
   \item{X}{optional numeric vector or matrix to specify one or more linear combinations of the coefficients in the model that should be tested. See \sQuote{Details}.}
   \item{att}{optional vector of indices to specify which scale coefficients should be included in the Wald-type test. Can also be a string to grep for. See \sQuote{Details}. Only relevant for location-scale models (see \code{\link{rma}}).}
   \item{Z}{optional numeric vector or matrix to specify one or more linear combinations of the scale coefficients in the model that should be tested. See \sQuote{Details}. Only relevant for location-scale models (see \code{\link{rma}}).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   When a single model is specified, the function provides Wald-type tests of one or more model coefficients or linear combinations thereof.

   In particular, for a fixed- or random-effects model (i.e., a model without moderators), this is just the test of the single coefficient of the model. For models including moderators, an omnibus test of all the model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first.

   Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, with \code{btt=c(3,4)}, only the third and fourth coefficient from the model would be included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model). Instead of specifying the coefficient numbers, one can specify a string for \code{btt}. In that case, \code{\link{grep}} will be used to search for all coefficient names that match the string. Using the \code{btt} argument, one can for example select all coefficients corresponding to a particular factor to test if the factor as a whole is significant. See \sQuote{Examples}.

   Instead, one can use the \code{X} argument to specify a linear combination of the coefficients in the model that should be tested (i.e., whether the linear combination is significantly different from zero). If a matrix of linear combinations is specified, each row defines a particular linear combination to be tested. If the matrix is of full rank, an omnibus Wald-type test of all linear combinations is also provided. Linear combinations can also be obtained with the \code{\link{predict.rma}} function, which will provide corresponding confidence intervals.

   When specifying two models for comparison, the function provides a likelihood ratio test (LRT) comparing the two models. The two models must be based on the same set of data, must be of the same class, and should be nested for the LRT to make sense. Note that LRTs are not meaningful when using REML estimation and the two models have different fixed effects.

   For location-scale models fitted with the \code{\link{rma}} function, one can use \code{att} to specify the indices of the scale coefficients to test. Similarly, one can use the \code{Z} argument to specify one or multiple linear combinations of the scale coefficients in the model that should be tested.
}
\value{
   An object of class \code{"anova.rma"}. When a single model is specified (without any further arguments or together with the \code{btt} argument), the object is a list containing the following components:
   \item{QM}{test statistic of the Wald-type test of the model coefficients.}
   \item{QMp}{corresponding p-value.}
   \item{btt}{indices of the coefficients tested by the Wald-type test.}
   \item{k}{number of outcomes included in the analysis.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the Wald-type test.}
   \item{\dots}{some additional elements/values.}
   When argument \code{X} is used, the object is a list containing the following components:
   \item{QM}{test statistic of the omnibus Wald-type test of all linear combinations.}
   \item{QMp}{corresponding p-value.}
   \item{hyp}{description of the linear combinations tested.}
   \item{Xb}{values of the linear combinations.}
   \item{se}{standard errors of the linear combinations.}
   \item{zval}{test statistics of the linear combinations.}
   \item{pval}{corresponding p-values.}
   When two models are specified, the object is a list containing the following components:
   \item{fit.stats.f}{log-likelihood, deviance, AIC, BIC, and AICc for the full model.}
   \item{fit.stats.r}{log-likelihood, deviance, AIC, BIC, and AICc for the reduced model.}
   \item{p.f}{number of parameters in the full model.}
   \item{p.r}{number of parameters in the reduced model.}
   \item{LRT}{likelihood ratio test statistic.}
   \item{pval}{corresponding p-value.}
   \item{QE.f}{test statistic of the test for (residual) heterogeneity from the full model.}
   \item{QE.r}{test statistic of the test for (residual) heterogeneity from the reduced model.}
   \item{tau2.f}{estimated \mjseqn{\tau^2} value from the full model. \code{NA} for \code{"rma.mv"} objects.}
   \item{tau2.r}{estimated \mjseqn{\tau^2} value from the reduced model. \code{NA} for \code{"rma.mv"} objects.}
   \item{R2}{amount (in percent) of the heterogeneity in the reduced model that is accounted for in the full model (\code{NA} for fixed-effects models or for \code{"rma.mv"} objects). This can be regarded as a pseudo \mjseqn{R^2} statistic (Raudenbush, 2009). Note that the value may not be very accurate unless \mjseqn{k} is large (Lopez-Lopez et al., 2014).}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{\link{print.anova.rma}} function.
}
\note{
   The function can also be used to conduct a likelihood ratio test (LRT) for the amount of (residual) heterogeneity in random- and mixed-effects models. The full model should then be fitted with either \code{method="ML"} or \code{method="REML"} and the reduced model with \code{method="FE"}. The p-value for the test is based on a chi-square distribution with 1 degree of freedom, but actually needs to be adjusted for the fact that the parameter (i.e., \mjseqn{\tau^2}) falls on the boundary of the parameter space under the null hypothesis (see Viechtbauer, 2007, for more details).

   LRTs for variance components in more complex models (as fitted with the \code{\link{rma.mv}} function) can also be conducted in this manner (see \sQuote{Examples}).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to meta-analysis with random effects. \emph{Statistics in Medicine}, \bold{15}(6), 619--629. \verb{https://doi.org/10.1002/(sici)1097-0258(19960330)15:6\%3C619::aid-sim188\%3E3.0.co;2-a}

   Huizenga, H. M., Visser, I., & Dolan, C. V. (2011). Testing overall and moderator effects in random effects meta-regression. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{64}(1), 1--19. \verb{https://doi.org/10.1348/000711010X522687}

   \enc{López-López}{Lopez-Lopez}, J. A., \enc{Marín-Martínez}{Marin-Martinez}, F., \enc{Sánchez-Meca}{Sanchez-Meca}, J., Van den Noortgate, W., & Viechtbauer, W. (2014). Estimation of the predictive power of the model in mixed-effects meta-regression: A simulation study. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{67}(1), 30--48. \verb{https://doi.org/10.1111/bmsp.12002}

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Viechtbauer, W. (2007). Hypothesis tests for population heterogeneity in meta-analysis. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{60}(1), 29--60. \verb{https://doi.org/10.1348/000711005X64042}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mv}}, \code{\link{print.anova.rma}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model
res1 <- rma(yi, vi, data=dat, method="ML")

### mixed-effects model with two moderators (absolute latitude and publication year)
res2 <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="ML")

### Wald-type test of the two moderators
anova(res2)

### alternative way of specifying the same test
anova(res2, X=rbind(c(0,1,0), c(0,0,1)))

### corresponding likelihood ratio test
anova(res1, res2)

### test of a linear combination
anova(res2, X=c(1,35,1970))

### use predict() to obtain the same linear combination (with its CI)
predict(res2, newmods=c(35,1970))

### mixed-effects model with three moderators
res3 <- rma(yi, vi, mods = ~ ablat + year + alloc, data=dat, method="ML")
res3

### test the 'alloc' factor
anova(res3, btt=4:5)

### instead of specifying the coefficient numbers, grep for "alloc"
anova(res3, btt="alloc")

############################################################################

### an example of doing LRTs of variance components in more complex models
dat <- dat.konstantopoulos2011
res <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat)

### test the district-level variance component
res0 <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat, sigma2=c(0,NA))
anova(res, res0)

### test the school-level variance component
res0 <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat, sigma2=c(NA,0))
anova(res, res0)

### test both variance components simultaneously
res0 <- rma.mv(yi, vi, data=dat)
anova(res, res0)
}
\keyword{models}
