\name{rma.uni}
\alias{rma.uni}
\alias{rma}
\title{Meta-Analysis via Linear (Mixed-Effects) Models}
\description{
   Function to fit the meta-analytic fixed- and random/mixed-effects models with or without moderators via linear (mixed-effects) models. See the documentation of the \pkg{\link{metafor-package}} for more details on these models. \loadmathjax
}
\usage{
rma.uni(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
        m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, r2i, ni, mods,
        measure="GEN", intercept=TRUE, data, slab, subset,
        add=1/2, to="only0", drop00=FALSE, vtype="LS",
        method="REML", weighted=TRUE, test="z",
        level=95, digits, btt, tau2, verbose=FALSE, control, \dots)
rma(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
        m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, r2i, ni, mods,
        measure="GEN", intercept=TRUE, data, slab, subset,
        add=1/2, to="only0", drop00=FALSE, vtype="LS",
        method="REML", weighted=TRUE, test="z",
        level=95, digits, btt, tau2, verbose=FALSE, control, \dots)
}
\arguments{
   \item{yi}{vector of length \mjseqn{k} with the observed effect sizes or outcomes. See \sQuote{Details}.}
   \item{vi}{vector of length \mjseqn{k} with the corresponding sampling variances. See \sQuote{Details}.}
   \item{sei}{vector of length \mjseqn{k} with the corresponding standard errors (only relevant when not using \code{vi}). See \sQuote{Details}.}
   \item{weights}{optional argument to specify a vector of length \mjseqn{k} with user-defined weights. See \sQuote{Details}.}
   \item{ai}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{bi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ci}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{di}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{x1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{x2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{t1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{t2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{m1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{m2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{sd1i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{sd2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{xi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{mi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ri}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ti}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{sdi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{r2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ni}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \mjseqn{k} rows and as many columns as there are moderator variables. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{measure}{character string to specify the type of data supplied to the function. When \code{measure="GEN"} (default), the observed effect sizes or outcomes and corresponding sampling variances (or standard errors) should be supplied to the function via the \code{yi}, \code{vi}, and \code{sei} arguments (only one of the two, \code{vi} or \code{sei}, needs to be specified). Alternatively, one can set \code{measure} to one of the effect size or outcome measures described under the documentation for the \code{\link{escalc}} function in which case one must specify the required data via the appropriate arguments.}
   \item{intercept}{logical to specify whether an intercept should be added to the model (the default is \code{TRUE}). Ignored when \code{mods} is a formula.}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \mjseqn{k} studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the analysis.}
   \item{add}{see the documentation of the \code{\link{escalc}} function.}
   \item{to}{see the documentation of the \code{\link{escalc}} function.}
   \item{drop00}{see the documentation of the \code{\link{escalc}} function.}
   \item{vtype}{see the documentation of the \code{\link{escalc}} function.}
   \item{method}{character string to specify whether a fixed- or a random/mixed-effects model should be fitted. A fixed-effects model (with or without moderators) is fitted when using \code{method="FE"}. Random/mixed-effects models are fitted by setting \code{method} equal to one of the following: \code{"DL"}, \code{"HE"}, \code{"SJ"}, \code{"ML"}, \code{"REML"}, \code{"EB"}, \code{"HS"}, \code{"HSk"}, or \code{"GENQ"}. Default is \code{"REML"}. See \sQuote{Details}.}
   \item{weighted}{logical to specify whether weighted (default) or unweighted estimation should be used to fit the model.}
   \item{test}{character string to specify how test statistics and confidence intervals for the fixed effects should be computed. By default (\code{test="z"}), Wald-type tests and CIs are obtained, which are based on a standard normal distribution. When \code{test="knha"}, the method by Knapp and Hartung (2003) is used for adjusting test statistics and confidence intervals. See \sQuote{Details}.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (if unspecified, the default is 4).}
   \item{btt}{optional vector of indices to specify which coefficients to include in the omnibus test of moderators. Can also be a string to grep for. See \sQuote{Details}.}
   \item{tau2}{optional numeric value to specify the amount of (residual) heterogeneity in a random- or mixed-effects model (instead of estimating it). Useful for sensitivity analyses (e.g., for plotting results as a function of \mjseqn{\tau^2}). When unspecified, the value of \mjseqn{\tau^2} is estimated from the data.}
   \item{verbose}{logical to specify whether output should be generated on the progress of the model fitting (the default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{control}{optional list of control values for the iterative estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{additional arguments.}
}
\details{

   \bold{Specifying the Data}

   The function can be used in conjunction with any of the usual effect size or outcome measures used in meta-analyses (e.g., log risk ratios, log odds ratios, risk differences, mean differences, standardized mean differences, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation, and so on). Simply specify the observed effect sizes or outcomes via the \code{yi} argument and the corresponding sampling variances via the \code{vi} argument. Instead of specifying \code{vi}, one can specify the standard errors (the square root of the sampling variances) via the \code{sei} argument. The \code{\link{escalc}} function can be used to compute a wide variety of effect size or outcome measures (and the corresponding sampling variances) based on summary statistics.

   Alternatively, the function can automatically calculate the values of a chosen effect size or outcome measure (and the corresponding sampling variances) when supplied with the necessary data. The \code{\link{escalc}} function describes which effect size or outcome measures are currently implemented and what data/arguments should then be specified/used. The \code{measure} argument should then be set to the desired effect size or outcome measure.

   \bold{Specifying the Model}

   The function can be used to fit fixed- and random/mixed-effects models, as well as meta-regression models including moderators (the difference between the various models is described in detail in the introductory \pkg{\link{metafor-package}} help file).

   Assuming the observed effect sizes or outcomes and corresponding sampling variances are supplied via \code{yi} and \code{vi}, a \emph{fixed-effects model} can be fitted with \code{rma(yi, vi, method="FE")}. Weighted estimation (with inverse-variance weights) is used by default. User-defined weights can be supplied via the \code{weights} argument. Unweighted estimation can be used by setting \code{weighted=FALSE} (same as setting the weights equal to a constant).

   A \emph{random-effects model} can be fitted with the same code but setting \code{method} to one of the various estimators for the amount of heterogeneity:
   \itemize{
   \item \code{method="DL"} = DerSimonian-Laird estimator
   \item \code{method="HE"} = Hedges estimator
   \item \code{method="HS"} = Hunter-Schmidt estimator
   \item \code{method="HSk"} = Hunter-Schmidt estimator with k-correction
   \item \code{method="SJ"} = Sidik-Jonkman estimator
   \item \code{method="ML"} = maximum-likelihood estimator
   \item \code{method="REML"} = restricted maximum-likelihood estimator
   \item \code{method="EB"} = empirical Bayes estimator
   \item \code{method="PM"} = Paule-Mandel estimator
   \item \code{method="GENQ"} = generalized Q-statistic estimator
   } For a description of the various estimators, see Brannick et al. (2019), DerSimonian and Kacker (2007), Raudenbush (2009), Viechtbauer (2005), and Viechtbauer et al. (2015). Note that the Hedges estimator is also called the \sQuote{variance component estimator} or \sQuote{Cochran estimator}, the Sidik-Jonkman estimator is also called the \sQuote{model error variance estimator}, and the empirical Bayes estimator is actually identical to the Paule-Mandel estimator (Paule & Mandel, 1982). Finally, the generalized Q-statistic estimator is a general method-of-moments estimator (DerSimonian & Kacker, 2007) requiring the specification of weights (the HE and DL estimators are just special cases with equal and inverse variance weights, respectively).

   One or more moderators can be included in these models via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate model matrix (i.e., \mjseqn{X}) with \mjseqn{k} rows and as many columns as there are moderator variables (e.g., \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for three moderator variables). The intercept is added to the model matrix by default unless \code{intercept=FALSE}.

   Alternatively, one can use the standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept is controlled by the specified formula (e.g., \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept).

   A \emph{fixed-effects with moderators model} is then fitted by setting \code{method="FE"}, while a \emph{mixed-effects model} is fitted by specifying one of the estimators for the amount of (residual) heterogeneity given earlier.

   When the observed effect sizes or outcomes and corresponding sampling variances are supplied via the \code{yi} and \code{vi} (or \code{sei}) arguments, one can also specify moderators via the \code{yi} argument (e.g., \code{rma(yi ~ mod1 + mod2 + mod3, vi)}). In that case, the \code{mods} argument is ignored and the inclusion/exclusion of the intercept again is controlled by the specified formula.

   \bold{Omnibus Test of Moderators}

   For models including moderators, an omnibus test of all model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, with \code{btt=c(3,4)}, only the third and fourth coefficient from the model would be included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model). Instead of specifying the coefficient numbers, one can specify a string for \code{btt}. In that case, \code{\link{grep}} will be used to search for all coefficient names that match the string. The omnibus test is called the \mjseqn{Q_M}-test and follows, under the assumptions of the model, a chi-square distribution with \mjseqn{m} degrees of freedom (with \mjseqn{m} denoting the number of coefficients tested) under the null hypothesis (that the true value of all coefficients tested is equal to 0).

   \bold{Categorical Moderators}

   Categorical moderator variables can be included in the model via the \code{mods} argument in the same way that appropriately (dummy) coded categorical independent variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to let \R handle the coding automatically. An example to illustrate these different approaches is provided below.

   \bold{Knapp & Hartung Adjustment}

   By default, tests of individual coefficients in the model (and the corresponding confidence intervals) are based on a standard normal distribution, while the omnibus test is based on a chi-square distribution (see above). The Knapp and Hartung (2003) method (\code{test="knha"}) is an adjustment to the standard errors of the estimated coefficients, which helps to account for the uncertainty in the estimate of the amount of (residual) heterogeneity and leads to different reference distributions. Tests of individual coefficients and confidence intervals are then based on a t-distribution with \mjseqn{k-p} degrees of freedom, while the omnibus test statistic then uses an F-distribution with \mjseqn{m} and \mjseqn{k-p} degrees of freedom (with \mjseqn{p} denoting the total number of model coefficients including the intercept if it is present). The Knapp and Hartung (2003) adjustment is only meant to be used in the context of random- or mixed-effects models.

   \bold{Test for (Residual) Heterogeneity}

   A test for (residual) heterogeneity is automatically carried out by the function. Without moderators in the model, this is simply Cochran's \mjseqn{Q}-test (Cochran, 1954), which tests whether the variability in the observed effect sizes or outcomes is larger than would be expected based on sampling variability alone. A significant test suggests that the true effects or outcomes are heterogeneous. When moderators are included in the model, this is the \mjseqn{Q_E}-test for residual heterogeneity, which tests whether the variability in the observed effect sizes or outcomes not accounted for by the moderators included in the model is larger than would be expected based on sampling variability alone.
}
\value{
   An object of class \code{c("rma.uni","rma")}. The object is a list containing the following components:
   \item{beta}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{tau2}{estimated amount of (residual) heterogeneity. Always \code{0} when \code{method="FE"}.}
   \item{se.tau2}{standard error of the estimated amount of (residual) heterogeneity.}
   \item{k}{number of studies included in the analysis.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of moderators.}
   \item{QE}{test statistic of the test for (residual) heterogeneity.}
   \item{QEp}{corresponding p-value.}
   \item{QM}{test statistic of the omnibus test of moderators.}
   \item{QMp}{corresponding p-value.}
   \item{I2}{value of \mjseqn{I^2}. See \code{\link{print.rma.uni}} for more details.}
   \item{H2}{value of \mjseqn{H^2}. See \code{\link{print.rma.uni}} for more details.}
   \item{R2}{value of \mjseqn{R^2}. See \code{\link{print.rma.uni}} for more details.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, vi, X}{the vector of outcomes, the corresponding sampling variances, and the model matrix.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values under the unrestricted and restricted likelihood.}
   \item{\dots}{some additional elements/values.}
}
\section{Methods}{
   The results of the fitted model are formatted and printed with the \code{\link{print.rma.uni}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them). Full versus reduced model comparisons in terms of fit statistics and likelihoods can be obtained with \code{\link{anova.rma}}. Wald-type tests for sets of model coefficients or linear combinations thereof can be obtained with the same function. Permutation tests for the model coefficient(s) can be obtained with \code{\link{permutest.rma.uni}}. Tests and confidence intervals based on (cluster) robust methods can be obtained with \code{\link{robust.rma.uni}}.

   Predicted/fitted values can be obtained with \code{\link{predict.rma}} and \code{\link{fitted.rma}}. For best linear unbiased predictions, see \code{\link{blup.rma.uni}} and \code{\link{ranef.rma.uni}}.

   The \code{\link{residuals.rma}}, \code{\link{rstandard.rma.uni}}, and \code{\link{rstudent.rma.uni}} functions extract raw and standardized residuals. Additional case diagnostics (e.g., to determine influential studies) can be obtained with the \code{\link{influence.rma.uni}} function. For models without moderators, leave-one-out diagnostics can also be obtained with \code{\link{leave1out.rma.uni}}. For models with moderators, variance inflation factors can be obtained with \code{\link{vif.rma}}.

   A confidence interval for the amount of (residual) heterogeneity in the random/mixed-effects model can be obtained with \code{\link{confint.rma.uni}}.

   Forest, funnel, radial, \enc{L'Abbé}{L'Abbe}, and Baujat plots can be obtained with \code{\link{forest.rma}}, \code{\link{funnel.rma}}, \code{\link{radial.rma}}, \code{\link{labbe.rma}}, and \code{\link{baujat}} (radial and \enc{L'Abbé}{L'Abbe} plots only for models without moderators). The \code{\link{qqnorm.rma.uni}} function provides normal QQ plots of the standardized residuals. One can also just call \code{\link{plot.rma.uni}} on the fitted model object to obtain various plots at once. For random/mixed-effects models, the \code{\link{profile.rma.uni}} function can be used to obtain a plot of the (restricted) log-likelihood as a function of \mjseqn{\tau^2}.

   Tests for funnel plot asymmetry (which may be indicative of publication bias) can be obtained with \code{\link{ranktest.rma}} and \code{\link{regtest.rma}}. For models without moderators, the \code{\link{trimfill.rma.uni}} method can be used to carry out a trim and fill analysis and \code{\link{hc.rma.uni}} provides a random-effects model analysis that is more robust to publication bias (based on the method by Henmi & Copas, 2010). The test of \sQuote{excess significance} can be carried out with the \code{\link{tes.rma}} function. Selection models can be fitted with the \code{\link{selmodel}} function.

   For models without moderators, a cumulative meta-analysis (i.e., adding one observation at a time) can be obtained with \code{\link{cumul.rma.uni}}.

   Other extractor functions include \code{\link{coef.rma}}, \code{\link{vcov.rma}}, \code{\link{logLik.rma}}, \code{\link{deviance.rma}}, \code{\link{AIC.rma}}, \code{\link{BIC.rma}}, \code{\link{hatvalues.rma.uni}}, and \code{\link{weights.rma.uni}}.
}
\note{
   While the HS, HSk, HE, DL, SJ, and GENQ estimators of \mjseqn{\tau^2} are based on closed-form solutions, the ML, REML, and EB estimators must be obtained numerically. For this, the function makes use of the Fisher scoring algorithm, which is robust to poor starting values and usually converges quickly (Harville, 1977; Jennrich & Sampson, 1976). By default, the starting value is set equal to the value of the Hedges (HE) estimator and the algorithm terminates when the change in the estimated value of \mjseqn{\tau^2} is smaller than \mjeqn{10^{-5}}{10^(-5)} from one iteration to the next. The maximum number of iterations is 100 by default (which should be sufficient in most cases). Information on the progress of the algorithm can be obtained by setting \code{verbose=TRUE}. One can also set \code{verbose} to an integer (\code{verbose=2} yields even more information and \code{verbose=3} also sets \code{option(warn=1)} temporarily).

   A different starting value, threshold, and maximum number of iterations can be specified via the \code{control} argument by setting \code{control=list(tau2.init=value, threshold=value, maxiter=value)}. The step length of the Fisher scoring algorithm can also be manually adjusted by a desired factor with \code{control=list(stepadj=value)} (values below 1 will reduce the step length). If using \code{verbose=TRUE} shows the estimate jumping around erratically (or cycling through a few values), decreasing the step length (and increasing the maximum number of iterations) can often help with convergence (e.g., \code{control=list(stepadj=0.5, maxiter=1000)}).

   The PM estimator also involves an iterative algorithm, which makes use of the \code{\link{uniroot}} function. By default, the desired accuracy (\code{tol}) is set equal to \code{.Machine$double.eps^0.25} and the maximum number of iterations (\code{maxiter}) to \code{100} (as above). The upper bound of the interval searched (\code{tau2.max}) is set to 100 (which should be large enough for most cases). These values can be adjusted with \code{control=list(tol=value, maxiter=value, tau2.max=value)}.

   All of the heterogeneity estimators except SJ can in principle yield negative estimates for the amount of (residual) heterogeneity. However, negative estimates of \mjseqn{\tau^2} are outside of the parameter space. For the HS, HSk, HE, DL, and GENQ estimators, negative estimates are therefore truncated to zero. For the ML, REML, and EB estimators, the Fisher scoring algorithm makes use of step halving (Jennrich & Sampson, 1976) to guarantee a non-negative estimate. Finally, for the PM estimator, the lower bound of the interval searched is set by default to zero. For those brave enough to step into risky territory, there is the option to set the lower bound for all these estimators equal to some other value besides zero (even a negative one) with \code{control=list(tau2.min=value)}, but the lowest value permitted is \code{-min(vi)} (to ensure that the marginal variances are always non-negative).

   The Hunter-Schmidt estimator for the amount of heterogeneity is defined in Hunter and Schmidt (1990) only in the context of the random-effects model when analyzing correlation coefficients. A general version of this estimator for random- and mixed-effects models not specific to any particular outcome measure is described in Viechtbauer (2005) and Viechtbauer et al. (2015) and is implemented here.

   The Sidik-Jonkman estimator starts with a crude estimate of \mjseqn{\tau^2}, which is then updated as described in Sidik and Jonkman (2005b, 2007). If, instead of the crude estimate, one wants to use a better a priori estimate, one can do so by passing this value via \code{control=list(tau2.init=value)}.

   Outcomes with non-positive sampling variances are problematic. If a sampling variance is equal to zero, then its weight will be \mjseqn{1/0} for fixed-effects models when using weighted estimation. Switching to unweighted estimation is a possible solution then. For random/mixed-effects model, some estimators of \mjseqn{\tau^2} are undefined when there is at least one sampling variance equal to zero. Other estimators may work, but it may still be necessary to switch to unweighted model fitting, especially when the estimate of \mjseqn{\tau^2} turns out to be zero.

   When including moderators in the model, it is possible that the model matrix is not of full rank (i.e., there is a linear relationship between the moderator variables included in the model). The function automatically tries to reduce the model matrix to full rank by removing redundant predictors, but if this fails the model cannot be fitted and an error will be issued. Deleting (redundant) moderator variables from the model as needed should solve this problem.

   Finally, some general words of caution about the assumptions underlying the models:
   \itemize{
   \item The sampling variances (i.e., the \code{vi} values) are treated as if they are known constants. This (usually) implies that the distributions of the test statistics and corresponding confidence intervals are only exact and have nominal coverage when the within-study sample sizes are large (i.e., when the error in the sampling variance estimates is small). Certain outcome measures (e.g., the arcsine square root transformed risk difference and the Fisher r-to-z transformed correlation coefficient) are based on variance stabilizing transformations that also help to make the assumption of known sampling variances much more reasonable.
   \item When fitting a mixed/random-effects model, \mjseqn{\tau^2} is estimated and then treated as a known constant thereafter. This ignores the uncertainty in the estimate of \mjseqn{\tau^2}. As a consequence, the standard errors of the parameter estimates tend to be too small, yielding test statistics that are too large and confidence intervals that are not wide enough. The Knapp and Hartung (2003) adjustment can be used to counter this problem, yielding test statistics and confidence intervals whose properties are closer to nominal.
   \item Most effect size or outcome measures are not exactly normally distributed as assumed under the various models. However, the normal approximation usually becomes more accurate for most effect size or outcome measures as the within-study sample sizes increase. Therefore, sufficiently large within-study sample sizes are (usually) needed to be certain that the tests and confidence intervals have nominal levels/coverage. Again, certain outcome measures (e.g., the Fisher r-to-z transformed correlation coefficient) may be preferable from this perspective as well.
   }
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Berkey, C. S., Hoaglin, D. C., Mosteller, F., & Colditz, G. A. (1995). A random-effects regression model for meta-analysis. \emph{Statistics in Medicine}, \bold{14}(4), 395--411. \verb{https://doi.org/10.1002/sim.4780140406}

   Brannick, M. T., Potter, S. M., Benitez, B., & Morris, S. B. (2019). Bias and precision of alternate estimators in meta-analysis: Benefits of blending Schmidt–Hunter and Hedges approaches. \emph{Organizational Research Methods}, \bold{22}(2), 490--514. \verb{https://doi.org/10.1177/1094428117741966}

   Cochran, W. G. (1954). The combination of estimates from different experiments. \emph{Biometrics}, \bold{10}(1), 101--129. \verb{https://doi.org/10.2307/3001666}

   DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials. \emph{Controlled Clinical Trials}, \bold{7}(3), 177--188. \verb{https://doi.org/10.1016/0197-2456(86)90046-2}

   DerSimonian, R., & Kacker, R. (2007). Random-effects model for meta-analysis of clinical trials: An update. \emph{Contemporary Clinical Trials}, \bold{28}(2), 105--114. \verb{https://doi.org/10.1016/j.cct.2006.04.004}

   Harville, D. A. (1977). Maximum likelihood approaches to variance component estimation and to related problems. \emph{Journal of the American Statistical Association}, \bold{72}(358), 320--338. \verb{https://doi.org/10.2307/2286796}

   Hedges, L. V. (1983). A random effects model for effect sizes. \emph{Psychological Bulletin}, \bold{93}(2), 388--395. \verb{https://doi.org/10.1037/0033-2909.93.2.388}

   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Henmi, M., & Copas, J. B. (2010). Confidence intervals for random effects meta-analysis and robustness to publication bias. \emph{Statistics in Medicine}, \bold{29}(29), 2969--2983. \verb{https://doi.org/10.1002/sim.4029}

   Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis: Correcting error and bias in research findings} (2nd ed.). Thousand Oaks, CA: Sage.

   Jennrich, R. I., & Sampson, P. F. (1976). Newton-Raphson and related algorithms for maximum likelihood variance component estimation. \emph{Technometrics}, \bold{18}(1), 11--17. \verb{https://doi.org/10.2307/1267911}

   Knapp, G., & Hartung, J. (2003). Improved tests for a random effects meta-regression with a single covariate. \emph{Statistics in Medicine}, \bold{22}(17), 2693--2710. \verb{https://doi.org/10.1002/sim.1482}

   Morris, C. N. (1983). Parametric empirical Bayes inference: Theory and applications. \emph{Journal of the American Statistical Association}, \bold{78}(381), 47--55. \verb{https://doi.org/10.2307/2287098}

   Paule, R. C., & Mandel, J. (1982). Consensus values and weighting factors. \emph{Journal of Research of the National Bureau of Standards}, \bold{87}(5), 377--385. \verb{https://doi.org/10.6028/jres.087.022}

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Sidik, K., & Jonkman, J. N. (2005a). A note on variance estimation in random effects meta-regression. \emph{Journal of Biopharmaceutical Statistics}, \bold{15}(5), 823--838. \verb{https://doi.org/10.1081/BIP-200067915}

   Sidik, K., & Jonkman, J. N. (2005b). Simple heterogeneity variance estimation for meta-analysis. \emph{Journal of the Royal Statistical Society, Series C}, \bold{54}(2), 367--384. \verb{https://doi.org/10.1111/j.1467-9876.2005.00489.x}

   Sidik, K., & Jonkman, J. N. (2007). A comparison of heterogeneity variance estimators in combining results of studies. \emph{Statistics in Medicine}, \bold{26}(9), 1964--1981. \verb{https://doi.org/10.1002/sim.2688}

   Viechtbauer, W. (2005). Bias and efficiency of meta-analytic variance estimators in the random-effects model. \emph{Journal of Educational and Behavioral Statistics}, \bold{30}(3), 261--293. \verb{https://doi.org/10.3102/10769986030003261}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., \enc{López-López}{Lopez-Lopez}, J. A., \enc{Sánchez-Meca}{Sanchez-Meca}, J., & \enc{Marín-Martínez}{Marin-Martinez}, F. (2015). A comparison of procedures to test for moderators in mixed-effects meta-regression models. \emph{Psychological Methods}, \bold{20}(3), 360--374. \verb{https://doi.org/10.1037/met0000023}
}
\seealso{
   \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for other model fitting functions.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model, using log risk ratios and variances as input
### note: method="REML" is the default, so one could leave this out
rma(yi, vi, data=dat, method="REML")

### random-effects model, using log risk ratios and standard errors as input
### note: the second argument of rma() is for the *variances*, so we use the
### named argument 'sei' to supply the standard errors to the function
dat$sei <- sqrt(dat$vi)
rma(yi, sei=sei, data=dat)

### supplying the 2x2 table cell frequencies directly to the rma() function
rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### mixed-effects model with two moderators (absolute latitude and publication year)
rma(yi, vi, mods=cbind(ablat, year), data=dat)

### using a model formula to specify the same model
rma(yi, vi, mods = ~ ablat + year, data=dat)

### using a model formula as part of the yi argument
rma(yi ~ ablat + year, vi, data=dat)

### manual dummy coding of the allocation factor
alloc.random     <- ifelse(dat$alloc == "random",     1, 0)
alloc.alternate  <- ifelse(dat$alloc == "alternate",  1, 0)
alloc.systematic <- ifelse(dat$alloc == "systematic", 1, 0)

### test the allocation factor (in the presence of the other moderators)
### note: 'alternate' is the reference level of the allocation factor,
###       since this is the dummy/level we leave out of the model
### note: the intercept is the first coefficient, so with btt=2:3 we test
###       coefficients 2 and 3, corresponding to the coefficients for the
###       allocation factor
rma(yi, vi, mods = ~ alloc.random + alloc.systematic + year + ablat, data=dat, btt=2:3)

### using a model formula to specify the same model
rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat, btt=2:3)

### test all pairwise differences with Holm's method (using the 'multcomp' package if installed)
res <- rma(yi, vi, mods = ~ factor(alloc) - 1, data=dat)
res
if (require(multcomp))
   summary(glht(res, linfct=contrMat(c("alternate"=1,"random"=1,"systematic"=1),
           type="Tukey")), test=adjusted("holm"))

### subgrouping versus using a single model with a factor (subgrouping provides
### an estimate of tau^2 within each subgroup, but the number of studies in each
### subgroup is quite small; the model with the allocation factor provides a
### single estimate of tau^2 based on a larger number of studies, but assumes
### that tau^2 is the same within each subgroup)
res.a <- rma(yi, vi, data=dat, subset=(alloc=="alternate"))
res.r <- rma(yi, vi, data=dat, subset=(alloc=="random"))
res.s <- rma(yi, vi, data=dat, subset=(alloc=="systematic"))
res.a
res.r
res.s
res <- rma(yi, vi, mods = ~ factor(alloc) - 1, data=dat)
res

### demonstrating that Q_E + Q_M = Q_Total for fixed-effects models
### note: this does not work for random/mixed-effects models, since Q_E and
### Q_Total are calculated under the assumption that tau^2 = 0, while the
### calculation of Q_M incorporates the estimate of tau^2
res <- rma(yi, vi, data=dat, method="FE")
res ### this gives Q_Total
res <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="FE")
res ### this gives Q_E and Q_M
res$QE + res$QM

### decomposition of Q_E into subgroup Q-values
res <- rma(yi, vi, mods = ~ factor(alloc), data=dat)
res

res.a <- rma(yi, vi, data=dat, subset=(alloc=="alternate"))
res.r <- rma(yi, vi, data=dat, subset=(alloc=="random"))
res.s <- rma(yi, vi, data=dat, subset=(alloc=="systematic"))

res.a$QE ### Q-value within subgroup "alternate"
res.r$QE ### Q-value within subgroup "random"
res.s$QE ### Q-value within subgroup "systematic"

res$QE
res.a$QE + res.r$QE + res.s$QE
}
\keyword{models}
