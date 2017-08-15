\name{rma.glmm}
\alias{rma.glmm}
\title{Meta-Analysis via Generalized Linear (Mixed-Effects) Models}
\description{Function to fit meta-analytic fixed- and random/mixed-effects models with or without moderators via generalized linear (mixed-effects) models. See below and the documentation of the \pkg{\link{metafor-package}} for more details on these models.}
\usage{
rma.glmm(ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, ti, ni,
         mods, measure, intercept=TRUE, data, slab, subset,
         add=1/2, to="only0", drop00=TRUE, vtype="LS",
         model="UM.FS", method="ML", test="z",
         level=95, digits=4, btt, nAGQ=7, verbose=FALSE, control, \dots)
}
\arguments{
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
   \item{xi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{mi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ti}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ni}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and as many columns as there are moderator variables. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{measure}{character string indicating the outcome measure to use for the meta-analysis. Possible options are the odds ratio (\code{"OR"}), the incidence rate ratio (\code{"IRR"}), the logit transformed proportion (\code{"PLO"}), or the log transformed incidence rate (\code{"IRLN"}).}
   \item{intercept}{logical indicating whether an intercept term should be added to the model (the default is \code{TRUE}).}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \eqn{k} studies.}
   \item{subset}{optional vector indicating the subset of studies that should be used for the analysis. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the observations to include.}
   \item{add}{non-negative number indicating the amount to add to zero cells, counts, or frequencies when calculating the observed outcomes of the individual studies. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{character string indicating when the values under \code{add} should be added (either \code{"only0"}, \code{"all"}, \code{"if0all"}, or \code{"none"}). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{drop00}{logical indicating whether studies with no cases/events (or only cases) in both groups should be dropped. See the documentation of the \code{\link{escalc}} function for more details.}
   \item{vtype}{character string indicating the type of sampling variances to calculate when calculating the observed outcomes. See the documentation of the \code{\link{escalc}} function for more details.}
   \item{model}{character string specifying the general model type to use for the analysis (either \code{"UM.FS"} (the default), \code{"UM.RS"}, \code{"CM.EL"}, or \code{"CM.AL"}). See \sQuote{Details}.}
   \item{method}{character string specifying whether a fixed- or a random/mixed-effects model should be fitted. A fixed-effects model (with or without moderators) is fitted when using \code{method="FE"}. Random/mixed-effects models are fitted by setting \code{method="ML"} (the default). See \sQuote{Details}.}
   \item{test}{character string specifying how test statistics and confidence intervals for the fixed effects should be computed. By default (\code{test="z"}), Wald-type tests and CIs are obtained, which are based on a standard normal distribution. When \code{test="t"}, a t-distribution is used instead. See \sQuote{Details}.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is 95).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is 4).}
   \item{btt}{optional vector of indices specifying which coefficients to include in the omnibus test of moderators. See \sQuote{Details}.}
   \item{nAGQ}{positive integer specifying the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. The default is 7. Setting this to 1 corresponds to the Laplacian approximation. See \sQuote{Note}.}
   \item{verbose}{logical indicating whether output should be generated on the progress of the model fitting (the default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{control}{optional list of control values for the estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{additional arguments.}
}
\details{

   \bold{Specifying the Data}

   The function can be used in conjunction with the following effect size or outcome measures:
   \itemize{
   \item \code{measure="OR"} for odds ratios (analyzed in log units)
   \item \code{measure="IRR"} for incidence rate ratios (analyzed in log units)
   \item \code{measure="PLO"} for logit transformed proportions (i.e., log odds)
   \item \code{measure="IRLN"} for log transformed incidence rates.
   }
   The \code{\link{escalc}} function describes the data/arguments that should be specified/used for these measures.

   \bold{Specifying the Model}

   A variety of model types are available when analyzing \eqn{2 \times 2}{2x2} table data (i.e., when \code{measure="OR"}) or two-group event count data (i.e., when \code{measure="IRR"}):
   \itemize{
   \item \code{model="UM.FS"} for an unconditional generalized linear mixed-effects model with fixed study effects
   \item \code{model="UM.RS"} for an unconditional generalized linear mixed-effects model with random study effects
   \item \code{model="CM.AL"} for a conditional generalized linear mixed-effects model (approximate likelihood)
   \item \code{model="CM.EL"} for a conditional generalized linear mixed-effects model (exact likelihood).
   }
   For \code{measure="OR"}, models \code{"UM.FS"} and \code{"UM.RS"} are essentially (mixed-effects) logistic regression models, while for \code{measure="IRR"}, these models are (mixed-effects) Poisson regression models. A choice must be made on how to model study level variability (i.e., differences in outcomes across studies irrespective of group membership). One can choose between using fixed study effects (which means that \eqn{k} dummy variables are added to the model) or random study effects (which means that random effects corresponding to the levels of the study factor are added to the model).

   The conditional model (\code{model="CM.EL"}) avoids having to model study level variability by conditioning on the total numbers of cases/events in each study. For \code{measure="OR"}, this leads to a non-central hypergeometric distribution for the data within each study and the corresponding model is then a (mixed-effects) conditional logistic model. Fitting this model can be difficult and computationally expensive. When the number of cases in each study is small relative to the group sizes, one can approximate the exact likelihood by a binomial distribution, which leads to a regular (mixed-effects) logistic regression model (\code{model="CM.AL"}). For \code{measure="IRR"}, the conditional model leads directly to a binomial distribution for the data within each study and the resulting model is again a (mixed-effects) logistic regression model (no approximate likelihood model is needed here).

   When analyzing proportions (i.e., \code{measure="PLO"}) or incidence rates (i.e., \code{measure="IRLN"}) of individual groups, the model type is always a (mixed-effects) logistic or Poisson regression model, respectively (i.e., the \code{model} argument is not relevant here).

   Aside from choosing the general model type, one has to decide whether to fit a fixed- or random-effects model to the data. A \emph{fixed-effects model} is fitted by setting \code{method="FE"}. A \emph{random-effects model} is fitted by setting \code{method="ML"} (the default). Note that random-effects models with dichotomous data are often referred to as \sQuote{binomial-normal} models in the meta-analytic literature. Analogously, for event count data, such models could be referred to as \sQuote{Poisson-normal} models.

   One or more moderators can be included in all of these models via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate model matrix (i.e., \eqn{\mathbf{X}}{X}) with \eqn{k} rows and as many columns as there are moderator variables (e.g., \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for three moderator variables). The intercept is added to the model matrix by default unless \code{intercept=FALSE}.

   Alternatively, one can use the standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept term is controlled by the specified formula (e.g., \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept term). With moderators, a \emph{fixed-effects with moderators model} is then fitted by setting \code{method="FE"}, while a \emph{mixed-effects model} is fitted by setting \code{method="ML"}.

   \bold{Fixed-, Saturated-, and Random/Mixed-Effects Models}

   When fitting a particular model, actually up to three different models are fitted within the function:
   \itemize{
   \item the fixed-effects model (i.e., where \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} is set to 0),
   \item the saturated model (i.e., the model with a deviance of 0), and
   \item the random/mixed-effects model (i.e., where \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} is estimated) (only if \code{method="ML"}).
   } The saturated model is obtained by adding as many dummy variables to the model as needed so that the model deviance is equal to zero. Even when \code{method="ML"}, the fixed-effects and saturated models are fitted, as they are used to compute the test statistics for the Wald-type and likelihood ratio tests for (residual) heterogeneity (see below).

   \bold{Omnibus Test of Parameters}

   For models including moderators, an omnibus test of all the model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, with \code{btt=c(3,4)}, only the third and fourth coefficient from the model would be included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model).

   \bold{Categorical Moderators}

   Categorical moderator variables can be included in the model via the \code{mods} argument in the same way that appropriately (dummy) coded categorical independent variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to let \R handle the coding automatically.

   \bold{Tests and Confidence Intervals}

   By default, the test statistics of the individual coefficients in the model (and the corresponding confidence intervals) are based on a standard normal distribution, while the omnibus test is based on a chi-square distribution with \eqn{m} degrees of freedom (\eqn{m} being the number of coefficients tested). As an alternative, one can set \code{test="t"}, which slightly mimics the Knapp and Hartung (2003) method by using a t-distribution with \eqn{k-p} degrees of freedom for tests of individual coefficients and confidence intervals and an F-distribution with \eqn{m} and \eqn{k-p} degrees of freedom (\eqn{p} being the total number of model coefficients including the intercept if it is present) for the omnibus test statistic.

   \bold{Tests for (Residual) Heterogeneity}

   Two different tests for (residual) heterogeneity are automatically carried out by the function. The first is a Wald-type test, which tests the coefficients corresponding to the dummy variables added in the saturated model for significance. The second is a likelihood ratio test, which tests the same set of coefficients, but does so by computing \eqn{-2} times the difference in the log-likelihood of the fixed-effects and the saturated model. These two tests are not identical for the types of models fitted by the \code{rma.glmm} function and may even lead to conflicting conclusions.

   \bold{Observed Outcomes of the Individual Studies}

   The various models do not require the calculation of the observed outcomes of the individual studies (e.g., the observed odds ratios of the \eqn{k} studies) and directly make use of the table/event counts. Zero cells/events are not a problem (except in extreme cases, such as when one of the two outcomes never occurs or when there are no events in any of the studies). Therefore, it is unnecessary to add some constant to the cell/event counts when there are zero cells/events.

   However, for plotting and various other functions, it is necessary to calculate the observed outcomes for the \eqn{k} studies. Here, zero cells/events can be problematic, so adding a constant value to the cell/event counts ensures that all \eqn{k} values can be calculated. The \code{add} and \code{to} arguments are used to specify what value should be added to the cell/event counts and under what circumstances when calculating the observed outcomes. The documentation of the \code{\link{escalc}} function explains how the \code{add} and \code{to} arguments work. Note that \code{drop00} is set to \code{TRUE} by default, since studies where \code{ai=ci=0} or \code{bi=di=0} or studies where \code{x1i=x2i=0} are uninformative about the size of the effect.
}
\value{
   An object of class \code{c("rma.glmm","rma")}. The object is a list containing the following components:
   \item{beta}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{p-values for the test statistics.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{tau2}{estimated amount of (residual) heterogeneity. Always \code{0} when \code{method="FE"}.}
   \item{sigma2}{estimated amount of study level variability (only for \code{model="UM.RS"}).}
   \item{k}{number of studies included in the model.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of coefficients.}
   \item{QE.Wld}{Wald-type test statistic for the test of (residual) heterogeneity.}
   \item{QEp.Wld}{p-value for the Wald-type test of (residual) heterogeneity.}
   \item{QE.LRT}{likelihood ratio test statistic for the test of (residual) heterogeneity.}
   \item{QEp.LRT}{p-value for the likelihood ratio test of (residual) heterogeneity.}
   \item{QM}{test statistic for the omnibus test of coefficients.}
   \item{QMp}{p-value for the omnibus test of coefficients.}
   \item{I2}{value of \if{latex}{\eqn{I^2}}\if{html}{\eqn{I²}}\if{text}{\eqn{I^2}}.}
   \item{H2}{value of \if{latex}{\eqn{H^2}}\if{html}{\eqn{H²}}\if{text}{\eqn{H^2}}.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, vi, X}{the vector of outcomes, the corresponding sampling variances, and the model matrix.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values.}
   \item{\dots}{some additional elements/values.}

   The results of the fitted model are formatted and printed with the \code{\link{print.rma.glmm}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them).
}
\note{
   Fitting the various types of models requires several different iterative algorithms:
   \itemize{
   \item For \code{model="UM.FS"} and \code{model="CM.AL"}, iteratively reweighted least squares (IWLS) as implemented in the \code{\link{glm}} function is used for fitting the fixed-effects and the saturated models. For \code{method="ML"}, adaptive Gauss-Hermite quadrature as implemented in the \code{\link[lme4]{glmer}} function is used. The same applies when \code{model="CM.EL"} is used in combination with \code{measure="IRR"} or when \code{measure="PLO"} or \code{measure="IRLN"} (regardless of the model type).
   \item For \code{model="UM.RS"}, adaptive Gauss-Hermite quadrature as implemented in the \code{\link[lme4]{glmer}} function is used to fit all of the models.
   \item For \code{model="CM.EL"} and \code{measure="OR"}, the quasi-Newton method (\code{"BFGS"}) as implemented in the \code{\link{optim}} function is used by default for fitting the fixed-effects and the saturated models. For \code{method="ML"}, the same algorithm is used, together with adaptive quadrature as implemented in the \code{\link{integrate}} function (for the integration over the density of the non-central hypergeometric distribution). Standard errors of the parameter estimates are obtained by inverting the Hessian, which is numerically approximated using the \code{\link[numDeriv]{hessian}} function.
   } When \code{model="CM.EL"} and \code{measure="OR"}, actually \code{model="CM.AL"} is used first to obtain starting values for \code{\link{optim}}, so either 4 (if \code{method="FE"}) or 6 (if \code{method="ML"}) models need to be fitted in total.

   Various control parameters can be adjusted via the \code{control} argument:
   \itemize{
      \item \code{optimizer} is set by default to \code{\link[stats:optim]{"optim"}}, but can be set to \code{\link[stats:nlminb]{"nlminb"}} or one of the optimizers from the \pkg{minqa} package (i.e., \code{\link[minqa:bobyqa]{"bobyqa"}}, \code{\link[minqa:newuoa]{"newuoa"}}, or \code{\link[minqa:uobyqa]{"uobyqa"}}),
      \item \code{optmethod} is used to set the \code{method} argument for \code{\link{optim}} (the default is \code{"BFGS"}),
      \item \code{optCtrl} is a list of named arguments to be passed on to the \code{control} argument of the chosen optimizer,
      \item \code{glmCtrl} is a list of named arguments to be passed on to the \code{control} argument of the \code{\link{glm}} function,
      \item \code{glmerCtrl} is a list of named arguments to be passed on to the \code{control} argument of the \code{\link[lme4]{glmer}} function, and
      \item \code{intCtrl} is a list of named arguments (i.e., \code{rel.tol} and \code{subdivisions}) to be passed on to the \code{\link{integrate}} function.
      \item \code{hessianCtrl} is a list of named arguments to be passed on to the \code{method.args} argument of the \code{\link[numDeriv]{hessian}} function. For some borderline cases, it may be necessary to bump up the \code{r} argument to a higher number to get sufficient accuracy when approximating the Hessian numerically (the default is \code{control=list(hessianCtrl=list(r=16))}).
   } Also, for \code{\link[lme4]{glmer}}, the \code{nAGQ} argument is used to specify the number of quadrature points. The default value is 7, which should provide sufficient accuracy in the evaluation of the log-likelihood in most cases, but at the expense of speed. Setting this to 1 corresponds to the Laplacian approximation (which is faster, but less accurate).

   Information on the progress of the various algorithms is obtained by setting \code{verbose=TRUE} or with \code{control=list(verbose=TRUE)}. Since fitting the various models can be computationally expensive, this option is quite useful to determine how the model fitting is progressing.

   For \code{model="CM.EL"} and \code{measure="OR"}, optimization involves repeated calculation of the density of the non-central hypergeometric distribution. When \code{method="ML"}, this also requires integration over the same density. This is currently implemented in a rather brute-force manner and may not be numerically stable, especially when models with moderators are fitted. Stability can be improved by scaling the moderators in a similar manner (i.e., don't use a moderator that is coded 0 and 1, while another uses values in the 1000s). For models with an intercept and moderators, the function actually rescales (non-dummy) variables to z-scores during the model fitting (results are given after back-scaling, so this should be transparent to the user). For models without an intercept, this is not done, so sensitivity analyses are highly recommended here (to ensure that the results do not depend on the scaling of the moderators).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}

   Code for computing the density of the non-central hypergeometric distribution from package \pkg{MCMCpack}, which in turn is based on Liao and Rosen (2001).
}
\references{
   Agresti, A. (2002). \emph{Categorical data analysis} (2nd. ed). Hoboken, NJ: Wiley.

   Bagos, P. G., & Nikolopoulos, G. K. (2009). Mixed-effects Poisson regression models for meta-analysis of follow-up studies with constant or varying durations. \emph{The International Journal of Biostatistics}, \bold{5}(1), article 21.

   van Houwelingen, H. C., Zwinderman, K. H., & Stijnen, T. (1993). A bivariate approach to meta-analysis. \emph{Statistics in Medicine}, \bold{12}, 2273--2284.

   Liao, J. G., & Rosen, O. (2001). Fast and stable algorithms for computing and sampling from the noncentral hypergeometric distribution. \emph{American Statistician}, \bold{55}, 366--369.

   Stijnen, T., Hamza, T. H., & Ozdemir, P. (2010). Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. \emph{Statistics in Medicine}, \bold{29}, 3046--3067.

   Turner, R. M., Omar, R. Z., Yang, M., Goldstein, H., & Thompson, S. G. (2000). A multilevel model framework for meta-analysis of clinical trials with binary outcomes. \emph{Statistics in Medicine}, \bold{19}, 3417--3432.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.mv}} for other model fitting functions.

   \code{\link{dat.nielweise2007}}, \code{\link{dat.nielweise2008}}, \code{\link{dat.collins1985a}}, and \code{link{dat.pritz1997}} for further examples of the use of the \code{rma.glmm} function.
}
\examples{
### random-effects model using rma.uni() (standard RE model analysis)
rma(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
    data=dat.bcg, method="ML", verbose=TRUE)

### random-effects models using rma.glmm() (require 'lme4' package)

### unconditional model with fixed study effects
\dontrun{
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
         data=dat.bcg, model="UM.FS", verbose=TRUE)}

### unconditional model with random study effects
\dontrun{
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
         data=dat.bcg, model="UM.RS", verbose=TRUE)}

### conditional model with approximate likelihood
\dontrun{
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
         data=dat.bcg, model="CM.AL", verbose=TRUE)}

### conditional model with exact likelihood (be patient!)
### note: fitting this model is very slow, so be patient
\dontrun{
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
         data=dat.bcg, model="CM.EL", verbose=TRUE)}
}
\keyword{models}
