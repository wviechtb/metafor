\name{rma.glmm}
\alias{rma.glmm}
\title{Meta-Analysis via Generalized Linear (Mixed-Effects) Models}
\description{
   Function to fit meta-analytic equal-, fixed-, and random-effects models and (mixed-effects) meta-regression models using a generalized linear (mixed-effects) model framework. See below and the introduction to the \pkg{\link{metafor-package}} for more details on these models. \loadmathjax
}
\usage{
rma.glmm(ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i, xi, mi, ti, ni,
         mods, measure, intercept=TRUE, data, slab, subset,
         add=1/2, to="only0", drop00=TRUE, vtype="LS",
         model="UM.FS", method="ML", coding=1/2, cor=FALSE, test="z",
         level=95, btt, nAGQ=7, verbose=FALSE, digits, control, \dots)
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
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \mjseqn{k} rows and as many columns as there are moderator variables. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{measure}{character string to specify the outcome measure to use for the meta-analysis. Possible options are \code{"OR"} for the (log transformed) odds ratio, \code{"IRR"} for the (log transformed) incidence rate ratio, \code{"PLO"} for the (logit transformed) proportion, or \code{"IRLN"} for the (log transformed) incidence rate.}
   \item{intercept}{logical to specify whether an intercept should be added to the model (the default is \code{TRUE}).}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \mjseqn{k} studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the analysis.}
   \item{add}{non-negative number to specify the amount to add to zero cells, counts, or frequencies when calculating the observed effect sizes or outcomes of the individual studies. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{character string to specify when the values under \code{add} should be added (either \code{"only0"}, \code{"all"}, \code{"if0all"}, or \code{"none"}). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{drop00}{logical to specify whether studies with no cases/events (or only cases) in both groups should be dropped. See the documentation of the \code{\link{escalc}} function for more details.}
   \item{vtype}{character string to specify the type of sampling variances to calculate when calculating the observed effect sizes or outcomes. See the documentation of the \code{\link{escalc}} function for more details.}
   \item{model}{character string to specify the general model type to use for the analysis. Either \code{"UM.FS"} (the default), \code{"UM.RS"}, \code{"CM.EL"}, or \code{"CM.AL"}. See \sQuote{Details}.}
   \item{method}{character string to specify whether an equal- or a random-effects model should be fitted. An equal-effects model is fitted when using \code{method="EE"}. A random-effects model is fitted by setting \code{method="ML"} (the default). See \sQuote{Details}.}
   \item{coding}{numeric scalar to indicate how the group variable should be coded in the random effects structure for random/mixed-effects models (the default is \code{1/2}). See \sQuote{Note}.}
   \item{cor}{logical to indicate whether the random study effects should be allowed to be correlated with the random group effects for random/mixed-effects models when \code{model="UM.RS"} (the default is \code{FALSE}). See \sQuote{Note}.}
   \item{test}{character string to specify how test statistics and confidence intervals for the fixed effects should be computed. By default (\code{test="z"}), Wald-type tests and CIs are obtained, which are based on a standard normal distribution. When \code{test="t"}, a t-distribution is used instead. See \sQuote{Details} and also \link[=misc-recs]{here} for some recommended practices.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95; see \link[=misc-options]{here} for details).}
   \item{btt}{optional vector of indices to specify which coefficients to include in the omnibus test of moderators. Can also be a string to \code{\link{grep}} for. See \sQuote{Details}.}
   \item{nAGQ}{positive integer to specify the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood. The default is 7. Setting this to 1 corresponds to the Laplacian approximation. See \sQuote{Note}.}
   \item{verbose}{logical to specify whether output should be generated on the progress of the model fitting (the default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is 4. See also \link[=misc-options]{here} for further details on how to control the number of digits in the output.}
   \item{control}{optional list of control values for the estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{additional arguments.}
}
\details{

   \subsection{Specifying the Data}{

      The function can be used in combination with the following effect sizes or outcome measures:
      \itemize{
      \item \code{measure="OR"} for (log transformed) odds ratios,
      \item \code{measure="IRR"} for (log transformed) incidence rate ratios,
      \item \code{measure="PLO"} for (logit transformed) proportions (i.e., log odds),
      \item \code{measure="IRLN"} for (log transformed) incidence rates.
      }
      The \code{\link{escalc}} function describes the data/arguments that should be specified/used for these measures.

   }

   \subsection{Specifying the Model}{

      A variety of model types are available when analyzing \mjeqn{2 \times 2}{2x2} table data (i.e., when \code{measure="OR"}) or two-group event count data (i.e., when \code{measure="IRR"}):
      \itemize{
      \item \code{model="UM.FS"} for an unconditional generalized linear mixed-effects model with fixed study effects,
      \item \code{model="UM.RS"} for an unconditional generalized linear mixed-effects model with random study effects,
      \item \code{model="CM.AL"} for a conditional generalized linear mixed-effects model (approximate likelihood),
      \item \code{model="CM.EL"} for a conditional generalized linear mixed-effects model (exact likelihood).
      }
      For \code{measure="OR"}, models \code{"UM.FS"} and \code{"UM.RS"} are essentially (mixed-effects) logistic regression models, while for \code{measure="IRR"}, these models are (mixed-effects) Poisson regression models. The difference between \code{"UM.FS"} and \code{"UM.RS"} is how study level variability (i.e., differences in outcomes across studies irrespective of group membership) is modeled. One can choose between using fixed study effects (which means that \mjseqn{k} dummy variables are added to the model) or random study effects (which means that random effects corresponding to the levels of the study factor are added to the model).

      The conditional model (\code{model="CM.EL"}) avoids having to model study level variability by conditioning on the total numbers of cases/events in each study. For \code{measure="OR"}, this leads to a non-central hypergeometric distribution for the data within each study and the corresponding model is then a (mixed-effects) conditional logistic model. Fitting this model can be difficult and computationally expensive. When the number of cases in each study is small relative to the group sizes, one can approximate the exact likelihood by a binomial distribution, which leads to a regular (mixed-effects) logistic regression model (\code{model="CM.AL"}). For \code{measure="IRR"}, the conditional model leads directly to a binomial distribution for the data within each study and the resulting model is again a (mixed-effects) logistic regression model (no approximate likelihood model is needed here).

      When analyzing proportions (i.e., \code{measure="PLO"}) or incidence rates (i.e., \code{measure="IRLN"}) of individual groups, the model type is always a (mixed-effects) logistic or Poisson regression model, respectively (i.e., the \code{model} argument is not relevant here).

      Aside from choosing the general model type, one has to decide whether to fit an equal- or a random-effects model to the data. An \emph{equal-effects model} is fitted by setting \code{method="EE"}. A \emph{random-effects model} is fitted by setting \code{method="ML"} (the default). Note that random-effects models with dichotomous data are often referred to as \sQuote{binomial-normal} models in the meta-analytic literature. Analogously, for event count data, such models could be referred to as \sQuote{Poisson-normal} models.

      One or more moderators can be included in a model via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate model matrix (i.e., \mjseqn{X}) with \mjseqn{k} rows and as many columns as there are moderator variables (e.g., \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for three moderator variables). The intercept is added to the model matrix by default unless \code{intercept=FALSE}.

      Alternatively, one can use standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept is controlled by the specified formula (e.g., \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept).

   }

   \subsection{Equal-, Saturated-, and Random/Mixed-Effects Models}{

      When fitting a particular model, actually up to three different models are fitted within the function:
      \itemize{
      \item the equal-effects model (i.e., where \mjseqn{\tau^2} is set to 0),
      \item the saturated model (i.e., the model with a deviance of 0), and
      \item the random/mixed-effects model (i.e., where \mjseqn{\tau^2} is estimated) (only if \code{method="ML"}).
      } The saturated model is obtained by adding as many dummy variables to the model as needed so that the model deviance is equal to zero. Even when \code{method="ML"}, the equal- and saturated models are also fitted, as they are used to compute the test statistics for the Wald-type and likelihood ratio tests for (residual) heterogeneity (see below).

   }

   \subsection{Omnibus Test of Moderators}{

      For models including moderators, an omnibus test of all model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} (\sQuote{betas to test}) argument (i.e., to test \mjseqn{\mbox{H}_0{:}\; \beta_{j \in \texttt{btt}} = 0}, where \mjseqn{\beta_{j \in \texttt{btt}}} is the set of coefficients to be tested). For example, with \code{btt=c(3,4)}, only the third and fourth coefficients from the model are included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model). Instead of specifying the coefficient numbers, one can specify a string for \code{btt}. In that case, \code{\link{grep}} will be used to search for all coefficient names that match the string. The omnibus test is called the \mjseqn{Q_M}-test and follows asymptotically a chi-square distribution with \mjseqn{m} degrees of freedom (with \mjseqn{m} denoting the number of coefficients tested) under the null hypothesis (that the true value of all coefficients tested is equal to 0).

   }

   \subsection{Categorical Moderators}{

      Categorical moderator variables can be included in the model via the \code{mods} argument in the same way that appropriately (dummy) coded categorical variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to automate the coding (note that string/character variables in a model formula are automatically converted to factors).

   }

   \subsection{Tests and Confidence Intervals}{

      By default, tests of individual coefficients in the model (and the corresponding confidence intervals) are based on a standard normal distribution, while the omnibus test is based on a chi-square distribution (see above). As an alternative, one can set \code{test="t"}, in which case tests of individual coefficients and confidence intervals are based on a t-distribution with \mjseqn{k-p} degrees of freedom, while the omnibus test then uses an F-distribution with \mjseqn{m} and \mjseqn{k-p} degrees of freedom (with \mjseqn{k} denoting the total number of estimates included in the analysis and \mjseqn{p} the total number of model coefficients including the intercept if it is present). Note that \code{test="t"} is not the same as \code{test="knha"} in \code{\link{rma.uni}}, as no adjustment to the standard errors of the estimated coefficients is made.

   }

   \subsection{Tests for (Residual) Heterogeneity}{

      Two different tests for (residual) heterogeneity are automatically carried out by the function. The first is a Wald-type test, which tests the coefficients corresponding to the dummy variables added in the saturated model for significance. The second is a likelihood ratio test, which tests the same set of coefficients, but does so by computing \mjseqn{-2} times the difference in the log-likelihoods of the equal-effects and the saturated models. These two tests are not identical for the types of models fitted by the \code{rma.glmm} function and may even lead to conflicting conclusions.

   }

   \subsection{Observed Effect Sizes or Outcomes of the Individual Studies}{

      The various models do not require the calculation of the observed effect sizes or outcomes of the individual studies (e.g., the observed log odds ratios of the \mjseqn{k} studies) and directly make use of the cell/event counts. Zero cells/events are not a problem (except in extreme cases, such as when one of the two outcomes never occurs or when there are no events in any of the studies). Therefore, it is unnecessary to add some constant to the cell/event counts when there are zero cells/events.

      However, for plotting and various other functions, it is necessary to calculate the observed effect sizes or outcomes for the \mjseqn{k} studies. Here, zero cells/events can be problematic, so adding a constant value to the cell/event counts ensures that all \mjseqn{k} values can be calculated. The \code{add} and \code{to} arguments are used to specify what value should be added to the cell/event counts and under what circumstances when calculating the observed effect sizes or outcomes. The documentation of the \code{\link{escalc}} function explains how the \code{add} and \code{to} arguments work. Note that \code{drop00} is set to \code{TRUE} by default, since studies where \code{ai=ci=0} or \code{bi=di=0} or studies where \code{x1i=x2i=0} are uninformative about the size of the effect.

   }

}
\value{
   An object of class \code{c("rma.glmm","rma")}. The object is a list containing the following components:
   \item{beta}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{tau2}{estimated amount of (residual) heterogeneity. Always \code{0} when \code{method="EE"}.}
   \item{sigma2}{estimated amount of study level variability (only for \code{model="UM.RS"}).}
   \item{k}{number of studies included in the analysis.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of moderators.}
   \item{QE.Wld}{Wald-type test statistic of the test for (residual) heterogeneity.}
   \item{QEp.Wld}{corresponding p-value.}
   \item{QE.LRT}{likelihood ratio test statistic of the test for (residual) heterogeneity.}
   \item{QEp.LRT}{corresponding p-value.}
   \item{QM}{test statistic of the omnibus test of moderators.}
   \item{QMp}{corresponding p-value.}
   \item{I2}{value of \mjseqn{I^2}.}
   \item{H2}{value of \mjseqn{H^2}.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, vi, X}{the vector of outcomes, the corresponding sampling variances, and the model matrix.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values.}
   \item{\dots}{some additional elements/values.}
}
\section{Methods}{
   The results of the fitted model are formatted and printed with the \code{\link[=print.rma.glmm]{print}} function. If fit statistics should also be given, use \code{\link[=summary.rma]{summary}} (or use the \code{\link[=fitstats.rma]{fitstats}} function to extract them).
}
\note{
   When \code{measure="OR"} or \code{measure="IRR"}, \code{model="UM.FS"} or \code{model="UM.RS"}, and \code{method="ML"}, one has to choose a coding scheme for the group variable in the random effects structure. When \code{coding=1/2} (the default), the two groups are coded with \code{+1/2} and \code{-1/2} (i.e., contrast coding), which is invariant under group label switching.

   When \code{coding=1}, the first group is coded with \code{1} and the second group with \code{0}. Finally, when \code{coding=0}, the first group is coded with \code{0} and the second group with \code{1}. Note that these coding schemes are not invariant under group label switching.

   When \code{model="UM.RS"} and \code{method="ML"}, one has to decide whether the random study effects are allowed to be correlated with the random group effects. By default (i.e., when \code{cor=FALSE}), no such correlation is allowed (which is typically an appropriate assumption when \code{coding=1/2}). When using a different coding scheme for the group variable (i.e., \code{coding=1} or \code{coding=0}), allowing the random study and group effects to be correlated (i.e., using \code{cor=TRUE}) is usually recommended.

   Fitting the various types of models requires several different iterative algorithms:
   \itemize{
   \item For \code{model="UM.FS"} and \code{model="CM.AL"}, iteratively reweighted least squares (IWLS) as implemented in the \code{\link{glm}} function is used for fitting the equal-effects and the saturated models. For \code{method="ML"}, adaptive Gauss-Hermite quadrature as implemented in the \code{\link[lme4]{glmer}} function is used. The same applies when \code{model="CM.EL"} is used in combination with \code{measure="IRR"} or when \code{measure="PLO"} or \code{measure="IRLN"} (regardless of the model type).
   \item For \code{model="UM.RS"}, adaptive Gauss-Hermite quadrature as implemented in the \code{\link[lme4]{glmer}} function is used to fit all of the models.
   \item For \code{model="CM.EL"} and \code{measure="OR"}, the quasi-Newton method optimizer as implemented in the \code{\link{nlminb}} function is used by default for fitting the equal-effects and the saturated models. For \code{method="ML"}, the same algorithm is used, together with adaptive quadrature as implemented in the \code{\link{integrate}} function (for the integration over the density of the non-central hypergeometric distribution). Standard errors of the parameter estimates are obtained by inverting the Hessian, which is numerically approximated using the \code{\link[numDeriv]{hessian}} function from the \code{numDeriv} package. One can also set \code{control=list(hesspack="pracma")} in which case the \code{\link[pracma]{hessian}} function from the \code{pracma} package is used instead for approximating the Hessian. When \mjseqn{\tau^2} is estimated to be smaller than \mjeqn{10^{-4}}{10^(-4)}, then \mjseqn{\tau^2} is effectively treated as zero for computing the standard errors (which helps to avoid numerical problems in approximating the Hessian). This cutoff can be adjusted via the \code{tau2tol} control argument (e.g., \code{control=list(tau2tol=1e-5)}).

   One can also chose a different optimizer from \code{\link{optim}} via the \code{control} argument (e.g., \code{control=list(optimizer="BFGS")} or \code{control=list(optimizer="Nelder-Mead")}). Besides \code{\link{nlminb}} and one of the methods from \code{\link{optim}}, one can also choose one of the optimizers from the \code{minqa} package (i.e., \code{\link[minqa]{uobyqa}}, \code{\link[minqa]{newuoa}}, or \code{\link[minqa]{bobyqa}}), one of the (derivative-free) algorithms from the \code{\link[nloptr]{nloptr}} package, the Newton-type algorithm implemented in \code{\link{nlm}}, the various algorithms implemented in the \code{dfoptim} package (\code{\link[dfoptim]{hjk}} for the Hooke-Jeeves, \code{\link[dfoptim]{nmk}} for the Nelder-Mead, and \code{\link[dfoptim]{mads}} for the Mesh Adaptive Direct Searches algorithm), the quasi-Newton type optimizers \code{\link[ucminf]{ucminf}} and \code{\link[lbfgsb3c]{lbfgsb3c}} and the subspace-searching simplex algorithm \code{\link[subplex]{subplex}} from the packages of the same name, the Barzilai-Borwein gradient decent method implemented in \code{\link[BB]{BBoptim}}, the \code{\link[Rcgmin]{Rcgmin}} and \code{\link[Rvmmin]{Rvmmin}} optimizers, or the parallelized version of the L-BFGS-B algorithm implemented in \code{\link[optimParallel]{optimParallel}} from the package of the same name.

   The optimizer name must be given as a character string (i.e., in quotes). Additional control parameters can be specified via the \code{optCtrl} elements of the \code{control} argument (e.g., \code{control=list(optCtrl=list(iter.max=1000, rel.tol=1e-8))}). For \code{\link[nloptr]{nloptr}}, the default is to use the BOBYQA implementation from that package with a relative convergence criterion of \code{1e-8} on the function value (i.e., log-likelihood), but this can be changed via the \code{algorithm} and \code{ftop_rel} arguments (e.g., \code{control=list(optimizer="nloptr", optCtrl=list(algorithm="NLOPT_LN_SBPLX", ftol_rel=1e-6))}). For \code{\link[optimParallel]{optimParallel}}, the control argument \code{ncpus} can be used to specify the number of cores to use for the parallelization (e.g., \code{control=list(optimizer="optimParallel", ncpus=2)}). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.} When \code{model="CM.EL"} and \code{measure="OR"}, actually \code{model="CM.AL"} is used first to obtain starting values for \code{\link{optim}}, so either 4 (if \code{method="EE"}) or 6 (if \code{method="ML"}) models need to be fitted in total.

   Various additional control parameters can be adjusted via the \code{control} argument:
   \itemize{
      \item \code{glmCtrl} is a list of named arguments to be passed on to the \code{control} argument of the \code{\link{glm}} function,
      \item \code{glmerCtrl} is a list of named arguments to be passed on to the \code{control} argument of the \code{\link[lme4]{glmer}} function,
      \item \code{intCtrl} is a list of named arguments (i.e., \code{rel.tol} and \code{subdivisions}) to be passed on to the \code{\link{integrate}} function, and
      \item \code{hessianCtrl} is a list of named arguments to be passed on to the \code{method.args} argument of the \code{\link[numDeriv]{hessian}} function. Most important is the \code{r} argument, which is set to 16 by default (i.e., \code{control=list(hessianCtrl=list(r=16))}). If the Hessian cannot be inverted, it may be necessary to adjust the \code{r} argument to a different number (e.g., try \code{r=4}, \code{r=6}, or \code{r=8}).
   } Also, for \code{\link[lme4]{glmer}}, the \code{nAGQ} argument is used to specify the number of quadrature points. The default value is 7, which should provide sufficient accuracy in the evaluation of the log-likelihood in most cases, but at the expense of speed. Setting this to 1 corresponds to the Laplacian approximation (which is faster, but less accurate). Note that \code{\link[lme4]{glmer}} does not allow values of \code{nAGQ > 1} when \code{model="UM.RS"} and \code{method="ML"}, so this value is automatically set to 1 for this model.

   Instead of \code{\link[lme4]{glmer}}, one can also choose to use \code{\link[GLMMadaptive]{mixed_model}} from the \code{GLMMadaptive} package or \code{\link[glmmTMB]{glmmTMB}} from the \code{glmmTMB} package for the model fitting. This is done by setting \code{control=list(package="GLMMadaptive")} or \code{control=list(package="glmmTMB")}, respectively.

   Information on the progress of the various algorithms can be obtained by setting \code{verbose=TRUE}. Since fitting the various models can be computationally expensive, this option is useful to determine how the model fitting is progressing. One can also set \code{verbose} to an integer (\code{verbose=2} yields even more information and \code{verbose=3} also sets \code{option(warn=1)} temporarily).

   For \code{model="CM.EL"} and \code{measure="OR"}, optimization involves repeated calculation of the density of the non-central hypergeometric distribution. When \code{method="ML"}, this also requires integration over the same density. This is currently implemented in a rather brute-force manner and may not be numerically stable, especially when models with moderators are fitted. Stability can be improved by scaling the moderators in a similar manner (i.e., don't use a moderator that is coded 0 and 1, while another uses values in the 1000s). For models with an intercept and moderators, the function actually rescales (non-dummy) variables to z-scores during the model fitting (results are given after back-scaling, so this should be transparent to the user). For models without an intercept, this is not done, so sensitivity analyses are highly recommended here (to ensure that the results do not depend on the scaling of the moderators). Also, if a warning is issued that the standard errors of the fixed effects are unusually small, one should try sensitivity analyses with different optimizers and/or adjusted settings for the \code{hessianCtrl} and \code{tau2tol} control arguments.

   Finally, there is also (experimental!) support for the following measures:
   \itemize{
   \item \code{measure="RR"} for log transformed risk ratios,
   \item \code{measure="RD"} for raw risk differences,
   \item \code{measure="PLN"} for log transformed proportions,
   \item \code{measure="PR"} for raw proportions,
   } (the first two only for models \code{"UM.FS"} and \code{"UM.RS"}) by using log and identity links for the binomial models. However, model fitting with these measures will often lead to numerical problems.

   Via the (undocumented) \code{link} argument, one can also directly adjust the link function that is used (by default, measures \code{"OR"} and \code{"PLO"} use a \code{"logit"} link, measures \code{"RR"} and \code{"PLN"} use a \code{"log"} link, measures \code{"RD"} and \code{"PR"} use an \code{"identity"} link, and measures \code{"IRR"} and \code{"IRLN"} use a \code{"log"} link). See \code{\link{family}} for alternative options. Changing these defaults is only recommended for users familiar with the consequences and the interpretation of the resulting estimates (when misused, the results could be meaningless).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}

   Code for computing the density of the non-central hypergeometric distribution comes from the \href{https://cran.r-project.org/package=MCMCpack}{MCMCpack} package, which in turn is based on Liao and Rosen (2001).
}
\references{
   Agresti, A. (2002). \emph{Categorical data analysis} (2nd. ed). Hoboken, NJ: Wiley.

   Bagos, P. G., & Nikolopoulos, G. K. (2009). Mixed-effects Poisson regression models for meta-analysis of follow-up studies with constant or varying durations. \emph{The International Journal of Biostatistics}, \bold{5}(1). \verb{https://doi.org/10.2202/1557-4679.1168}

   van Houwelingen, H. C., Zwinderman, K. H., & Stijnen, T. (1993). A bivariate approach to meta-analysis. \emph{Statistics in Medicine}, \bold{12}(24), 2273--2284. \verb{https://doi.org/10.1002/sim.4780122405}

   Jackson, D., Law, M., Stijnen, T., Viechtbauer, W., & White, I. R. (2018). A comparison of seven random-effects models for meta-analyses that estimate the summary odds ratio. \emph{Statistics in Medicine}, \bold{37}(7), 1059-1085. \verb{https://doi.org/10.1002/sim.7588}

   Liao, J. G., & Rosen, O. (2001). Fast and stable algorithms for computing and sampling from the noncentral hypergeometric distribution. \emph{American Statistician}, \bold{55}(4), 366--369. \verb{https://doi.org/10.1198/000313001753272547}

   Simmonds, M. C., & Higgins, J. P. T. (2016). A general framework for the use of logistic regression models in meta-analysis. \emph{Statistical Methods in Medical Research}, \bold{25}(6), 2858--2877. \verb{https://doi.org/10.1177/0962280214534409}

   Stijnen, T., Hamza, T. H., & Ozdemir, P. (2010). Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. \emph{Statistics in Medicine}, \bold{29}(29), 3046--3067. \verb{https://doi.org/10.1002/sim.4040}

   Turner, R. M., Omar, R. Z., Yang, M., Goldstein, H., & Thompson, S. G. (2000). A multilevel model framework for meta-analysis of clinical trials with binary outcomes. \emph{Statistics in Medicine}, \bold{19}(24), 3417--3432. \verb{https://doi.org/10.1002/1097-0258(20001230)19:24<3417::aid-sim614>3.0.co;2-l}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.mv}} for other model fitting functions.

   \code{\link[metadat]{dat.nielweise2007}}, \code{\link[metadat]{dat.nielweise2008}}, \code{\link[metadat]{dat.collins1985a}}, and \code{\link[metadat]{dat.pritz1997}} for further examples of the use of the \code{rma.glmm} function.
}
\examples{
############################################################################

### random-effects model using rma.uni() (standard RE model analysis)
rma(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, method="ML")

### random-effects models using rma.glmm() (requires 'lme4' package)

\dontrun{
### unconditional model with fixed study effects (the default)
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, model="UM.FS")

### unconditional model with random study effects
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, model="UM.RS")

### conditional model with approximate likelihood
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, model="CM.AL")

### conditional model with exact likelihood
### note: fitting this model may take a bit of time, so be patient
rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, model="CM.EL")
}

############################################################################

### try some alternative measures

\dontrun{
rma.glmm(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
rma.glmm(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
}

############################################################################

### meta-analysis of proportions

\dontrun{
dat <- dat.debruin2009

### binomial-normal model (with logit link) = mixed-effects logistic model
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat)
predict(res, transf=transf.ilogit)

### binomial-normal model with measure="PLN" (uses a log link)
res <- rma.glmm(measure="PLN", xi=xi, ni=ni, data=dat)
predict(res, transf=exp)

### binomial-normal model with measure="PR" (uses an identity link)
res <- rma.glmm(measure="PR", xi=xi, ni=ni, data=dat)
predict(res)

### binomial-normal model (with probit link) = mixed-effects probit model
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat, link="probit")
predict(res, transf=pnorm)

### further link functions that one could consider here
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat, link="cauchit")
predict(res, transf=pcauchy)
res <- rma.glmm(measure="PLO", xi=xi, ni=ni, data=dat, link="cloglog")
predict(res, transf=\(x) 1-exp(-exp(x)))
}

############################################################################
}
\keyword{models}
