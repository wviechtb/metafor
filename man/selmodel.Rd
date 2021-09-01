\name{selmodel}
\alias{selmodel}
\alias{selmodel.rma.uni}
\title{Selection Models}
\description{
   Function to fit selection models. \loadmathjax
}
\usage{
selmodel(x, \dots)

\method{selmodel}{rma.uni}(x, type, alternative="greater", prec,
         delta, steps, verbose=FALSE, digits, control, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"}.}
   \item{type}{character string to specify the type of selection model. Possible options are \code{"beta"}, \code{"halfnorm"}, \code{"negexp"}, \code{"logistic"}, \code{"power"}, \code{"negexppow"}, or \code{"stepfun"}. Can be abbreviated. See \sQuote{Details}.}
   \item{alternative}{character string to specify the sidedness of the hypothesis when testing the observed outcomes. Possible options are \code{"greater"} (the default), \code{"less"}, or \code{"two.sided"}. Can be abbreviated.}
   \item{prec}{optional character string to specify the measure of precision (only relevant for selection models that can incorporate this into the selection function). Possible options are \code{"sei"}, \code{"vi"}, \code{"ninv"}, or \code{"sqrtninv"}. See \sQuote{Details}.}
   \item{delta}{optional numeric vector (of the same length as the number of selection model parameters) to fix the corresponding \mjseqn{\delta} value(s). A specific \mjseqn{\delta} value can be fixed by setting the corresponding element of this argument to the desired value. A specific \mjseqn{\delta} value will be estimated if the corresponding element is set equal to \code{NA}. See \sQuote{Details}.}
   \item{steps}{numeric vector of one or more values between 0 and 1 that can or must be specified for certain selection functions. See \sQuote{Details}.}
   \item{verbose}{logical to specify whether output should be generated on the progress of the model fitting (the default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{control}{optional list of control values for the estimation algorithm. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{other arguments.}
}
\details{
   Selection models are a general class of models that attempt to model the process by which the studies included in a meta-analysis may have been influenced by some form of publication bias. If a particular selection model is an adequate approximation for the underlying selection process, then the model provides estimates of the parameters of interest (e.g., the average true outcome and the amount of heterogeneity in the true outcomes) that are \sQuote{corrected} for this selection process (i.e., they are estimates of the parameters in the population of studies before any selection has taken place). The present function fits a variety of such selection models. To do so, one should pass an object fitted with the \code{\link{rma.uni}} function to the first argument. The model that will then be fitted is of the same form as the original model combined with the specific selection model chosen (see below for possible options). For example, if the original model was a random-effects model, then a random-effects selection model will be fitted. Similarly, if the original model included moderators, then they will also be included in the selection model. Model fitting is done via maximum likelihood (ML) estimation over the fixed- and random-effects parameters (e.g., \mjseqn{\mu} and \mjseqn{\tau^2} in a random-effects model) and the selection model parameters.

   Argument \code{type} determines the specific type of selection model that should be fitted. All selection models that can be fitted are based on the idea that selection may haven taken place based on the p-values of the studies. In particular, let \mjseqn{y_i} and \mjseqn{v_i} denote the observed outcome and the corresponding sampling variance of the \mjseqn{i}th study. Then \mjseqn{z_i = y_i / \sqrt{v_i}} is the (Wald-type) test statistic for testing the null hypothesis \mjeqn{\mbox{H}_0{:}\; \theta_i = 0}{H_0: \theta_i = 0} and \mjseqn{p_i = 1 - \Phi(z_i)} (if \code{alternative="greater"}), \mjseqn{p_i = \Phi(z_i)} (if \code{alternative="less"}), or \mjseqn{p_i = 2(1 - \Phi(|z_i|))} (if \code{alternative="two.sided"}) the corresponding (one- or two-sided) p-value, where \mjseqn{\Phi()} denotes the cumulative distribution function of a standard normal distribution. Finally, let \mjseqn{w(p_i)} denote some function that specifies the relative likelihood of selection given the p-value of a study.

   If \mjseqn{w(p_i) > w(p_{i'})} when \mjseqn{p_i < p_{i'}} (i.e., \mjseqn{w(p_i)} is larger for smaller p-values), then \code{alternative="greater"} implies selection in favor of increasingly significant positive outcomes, \code{alternative="less"} implies selection in favor of increasingly significant negative outcomes, and \code{alternative="two.sided"} implies selection in favor of increasingly significant outcomes regardless of their direction.

   \subsection{Beta Selection Model}{

   When \code{type="beta"}, the function can be used to fit the \sQuote{beta selection model} by Citkowicz and Vevea (2017). For this model, the selection function is given by \mjsdeqn{w(p_i) = p_i^{\delta_1 - 1} \times (1 - p_i)^{\delta_2 - 1}} where \mjseqn{\delta_1 > 0} and \mjseqn{\delta_2 > 0}. The null hypothesis \mjeqn{\mbox{H}_0{:}\; \delta_1 = \delta_2 = 1}{H_0: \delta_1 = \delta_2 = 1} represents the case where there is no selection (at least not depending on the p-values). \ifelse{text}{}{The figure below illustrates with some examples how the relative likelihood of selection can depend on the p-value for various combinations of \mjseqn{\delta_1} and \mjseqn{\delta_2}.} Note that the model allows for a non-monotonic selection function.

   \if{html}{\figure{selmodel-beta.png}{options: width=600}}
   \if{latex}{\figure{selmodel-beta.pdf}{options: width=4in}}

   }

   \subsection{Half-Normal, Negative-Exponential, Logistic, and Power Selection Models}{

   Preston et al. (2004) suggested the first three of the following selection functions:

   \tabular{lllll}{
   \bold{name}          \tab \ics \tab \bold{\code{type}} \tab \ics \tab \bold{selection function} \cr
   half-normal          \tab \ics \tab \code{"halfnorm"}  \tab \ics \tab \mjseqn{w(p_i) = \exp(-\delta \times p_i^2)} \cr
   negative-exponential \tab \ics \tab \code{"negexp"}    \tab \ics \tab \mjseqn{w(p_i) = \exp(-\delta \times p_i)} \cr
   logistic             \tab \ics \tab \code{"logistic"}  \tab \ics \tab \mjseqn{w(p_i) = 2 \times \exp(-\delta \times p_i) / (1 + \exp(-\delta \times p_i))} \cr
   power                \tab \ics \tab \code{"power"}     \tab \ics \tab \mjseqn{w(p_i) = (1-p_i)^\delta}
   }

   The power selection model is added here as it has similar properties as the models suggested by Preston et al. (2004). For all models, assume \mjseqn{\delta \ge 0}, so that all functions imply a monotonically decreasing relationship between the p-value and the selection probability. For all functions, \mjeqn{\mbox{H}_0{:}\; \delta = 0}{H_0: \delta = 0} implies no selection. \ifelse{text}{}{The figure below shows the relative likelihood of selection as a function of the p-value for \mjseqn{\delta = 0} and for the various selection functions when \mjseqn{\delta = 6}.}

   \if{html}{\figure{selmodel-preston.png}{options: width=600}}
   \if{latex}{\figure{selmodel-preston.pdf}{options: width=4in}}

   Here, these functions are extended to allow for the possibility that \mjseqn{w(p_i) = 1} for p-values below a certain significance threshold denoted by \mjseqn{\alpha} (e.g., to model the case that the relative likelihood of selection is equally high for all significant studies but decreases monotonically for p-values above the significance threshold). To fit such a selection model, one should specify the \mjseqn{\alpha} value (with \mjseqn{0 < \alpha < 1}) via the \code{steps} argument. There must be at least one observed p-value below and above the chosen threshold to fit these models. \ifelse{text}{}{The figure below shows some examples of the relative likelihood of selection when \code{steps=.05}.}

   \if{html}{\figure{selmodel-preston-step.png}{options: width=600}}
   \if{latex}{\figure{selmodel-preston-step.pdf}{options: width=4in}}

   Preston et al. (2004) also suggested selection functions where the relatively likelihood of selection not only depends on the p-value, but also the precision (e.g., standard error) of the estimate (if two studies have similar p-values, it may be plausible to assume that the larger / more precise study has a higher probability of selection). These selection functions plus the corresponding power functions are given by:

   \tabular{lllll}{
   \bold{name}          \tab \ics \tab \bold{\code{type}} \tab \ics \tab \bold{selection function} \cr
   half-normal          \tab \ics \tab \code{"halfnorm"}  \tab \ics \tab \mjseqn{w(p_i) = \exp(-\delta \times \mathrm{prec}_i \times p_i^2)} \cr
   negative-exponential \tab \ics \tab \code{"negexp"}    \tab \ics \tab \mjseqn{w(p_i) = \exp(-\delta \times \mathrm{prec}_i \times p_i)} \cr
   logistic             \tab \ics \tab \code{"logistic"}  \tab \ics \tab \mjseqn{w(p_i) = 2 \times \exp(-\delta \times \mathrm{prec}_i \times p_i) / (1 + \exp(-\delta \times \mathrm{prec}_i \times p_i))} \cr
   power                \tab \ics \tab \code{"power"}     \tab \ics \tab \mjseqn{w(p_i) = (1-p_i)^{\delta \times \mathrm{prec}_i}}
   }

   where \mjseqn{\mathrm{prec}_i = \sqrt{v_i}} (i.e., the standard error of the \mjseqn{i}th study) according to Preston et al. (2004). Here, this idea is generalized to allow the user to specify the specific measure of precision to use (via the \code{prec} argument). Possible options are:
   \itemize{
   \item \code{prec="sei"} for the standard errors,
   \item \code{prec="vi"} for the sampling variances,
   \item \code{prec="ninv"} for the inverse of the sample sizes,
   \item \code{prec="sqrtninv"} for the inverse square root of the sample sizes.
   } Using some function of the sample sizes as a measure of precision is only possible when information about the sample sizes is actually stored within the object passed to the \code{selmodel} function. See \sQuote{Note}.

   Note that \mjseqn{\mathrm{prec}_i} is really a measure of imprecision (with higher values corresponding to lower precision). Also, regardless of the specific measure chosen, the values are actually rescaled with \mjseqn{\mathrm{prec}_i = \mathrm{prec}_i / \max(\mathrm{prec}_i)} inside of the function, such that \mjseqn{\mathrm{prec}_i = 1} for the least precise study and \mjseqn{\mathrm{prec}_i < 1} for the remaining studies (the rescaling does not actually change the fit of the model, it only helps to improve the stability of model fitting algorithm). \ifelse{text}{}{The figure below shows some examples of the relative likelihood of selection using these selection functions for two different precision values.}

   \if{html}{\figure{selmodel-preston-prec.png}{options: width=600}}
   \if{latex}{\figure{selmodel-preston-prec.pdf}{options: width=4in}}

   One can also use the \code{steps} argument as described above in combination with these selection functions (studies with p-values below the chosen threshold then have \mjseqn{w(p_i) = 1} regardless of their exact p-value or precision).

   }

   \subsection{Negative Exponential Power Selection Model}{

   As an extension of the half-normal and negative-exponential models, one can also choose \code{type="negexppow"} for a \sQuote{negative exponential power selection model}. The selection function is then given by \mjsdeqn{w(p_i) = \exp(-\delta_1 \times p_i^{1/\delta_2})} where \mjseqn{\delta_1 \ge 0} and \mjseqn{\delta_2 \ge 0} (see Begg & Mazumdar, 1994, although here a different parameterization is used, such that increasing \mjseqn{\delta_2} leads to more severe selection). \ifelse{text}{}{The figure below shows some examples of this selection function when holding \mjseqn{\delta_1} constant while increasing \mjseqn{\delta_2}.}

   \if{html}{\figure{selmodel-negexppow.png}{options: width=600}}
   \if{latex}{\figure{selmodel-negexppow.pdf}{options: width=4in}}

   This model affords greater flexibility in the shape of the selection function, but requires the estimation of the additional power parameter (the half-normal and negative-exponential models are therefore special cases when fixing \mjseqn{\delta_2} to 0.5 or 1, respectively). \mjeqn{\mbox{H}_0{:}\; \delta_1 = 0}{H_0: \delta_1 = 0} again implies no selection, but so does \mjeqn{\mbox{H}_0{:}\; \delta_2 = 0}{H_0: \delta_2 = 0}.

   One can again use the \code{steps} argument to specify a single significance threshold, \mjseqn{\alpha}, so that \mjseqn{w(p_i) = 1} for p-values below this threshold and otherwise \mjseqn{w(p_i)} follows the selection function as given above. One can also use the \code{prec} argument to specify a measure of precision in combination with this model, which leads to the selection function \mjsdeqn{w(p_i) = \exp(-\delta_1 \times \mathrm{prec}_i \times p_i^{1/\delta_2})} and hence is the logical extension of the negative exponential power selection model that also incorporates some measure of precision into the selection process.

   }

   \subsection{Step Function Selection Models}{

   When \code{type="stepfun"}, the function can be used to fit \sQuote{step function models} as described by Iyengar and Greenhouse (1988), Hedges (1992), Vevea and Hedges (1995), and Vevea and Woods (2005). For these models, one must specify one or multiple values via the \code{steps} argument, which define intervals in which the relative likelihood of selection is constant. Let \mjsdeqn{\alpha_1 < \alpha_2 < \ldots < \alpha_c} denote these cutpoints sorted in increasing order, with the constraint that \mjseqn{\alpha_c = 1} (if the highest value specified via \code{steps} is not 1, the function will automatically add this cutpoint), and define \mjseqn{\alpha_0 = 0}. The selection function is then given by \mjseqn{w(p_i) = \delta_j} if \mjseqn{\alpha_{j-1} < p_i \le \alpha_j}. To make the model identifiable, we set \mjseqn{\delta_1 = 1}. The \mjseqn{\delta_j} values therefore denote the likelihood of selection in the various intervals relative to the interval for p-values between 0 and \mjseqn{\alpha_1}. Hence, the null hypothesis \mjeqn{\mbox{H}_0{:}\; \delta_j = 1}{H_0: \delta_j = 1} for \mjseqn{j = 1, \ldots, c} implies no selection.

   For example, if \code{steps=c(.05, .10, .50, 1)}, then \mjseqn{\delta_2} is the likelihood of selection for p-values between .05 and .10, \mjseqn{\delta_3} is the likelihood of selection for p-values between .10 and .50, and \mjseqn{\delta_4} is the likelihood of selection for p-values between .50 and 1 relative to the likelihood of selection for p-values between 0 and .05. \ifelse{text}{}{The figure below shows the corresponding selection function for some arbitrarily chosen \mjseqn{\delta_j} values.}

   \if{html}{\figure{selmodel-stepfun.png}{options: width=600}}
   \if{latex}{\figure{selmodel-stepfun.pdf}{options: width=4in}}

   There must be at least one observed p-value within each interval to fit this model. If this is not the case, an error will be issued (setting \code{verbose=TRUE} provides information about the number of p-values falling into each interval).

   When specifying a single cutpoint in the context of a random-effects model, this model is sometimes called the \sQuote{three-parameter selection model} (3PSM), corresponding to the parameters \mjseqn{\mu}, \mjseqn{\tau^2}, and \mjseqn{\delta_2} (e.g., Carter et al., 2019; McShane et al., 2016; Pustejovsky & Rodgers, 2019). The same idea but in the context of a fixed-effects model was also described by Iyengar and Greenhouse (1988).

   Note that when \code{alternative="greater"} or \code{alternative="less"} (i.e., when we assume that the relative likelihood of selection is not only related to the p-values of the studies, but also the directionality of the outcomes), then it would usually make sense to divide conventional levels of significance (e.g., .05) by 2 before passing these values to the \code{steps} argument. For example, if we think that studies were selected for positive outcomes that are significant at two-tailed \mjseqn{\alpha = .05}, then we should use \code{alternative="greater"} in combination with \code{steps=c(.025, 1)}.

   One of the challenges when fitting this model with many cutpoints is the large number of parameters that need to be estimated (which is especially problematic when the number of studies is small). An alternative approach suggested by Vevea and Woods (2005) is to fix the \mjseqn{\delta_j} values to some a priori chosen values instead of estimating them. One can then conduct a sensitivity analysis by examining the results (e.g., the estimates of \mjseqn{\mu} and \mjseqn{\tau^2} in a random-effects model) for a variety of different sets of \mjseqn{\delta_j} values (reflecting more or less severe forms of selection). This can be done by specifying the \mjseqn{\delta_j} values via the \code{delta} argument. Table 1 in Vevea and Woods (2005) provides some illustrative examples of moderate and severe selection functions for one- and two-tailed selection. The code below creates a data frame that contains these functions.

   \preformatted{tab <- data.frame(
   steps = c(0.005, 0.01, 0.05, 0.10, 0.25, 0.35, 0.50, 0.65, 0.75, 0.90, 0.95, 0.99, 0.995, 1),
   delta.mod.1 = c(1, 0.99, 0.95, 0.80, 0.75, 0.65, 0.60, 0.55, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50),
   delta.sev.1 = c(1, 0.99, 0.90, 0.75, 0.60, 0.50, 0.40, 0.35, 0.30, 0.25, 0.10, 0.10, 0.10, 0.10),
   delta.mod.2 = c(1, 0.99, 0.95, 0.90, 0.80, 0.75, 0.60, 0.60, 0.75, 0.80, 0.90, 0.95, 0.99, 1.00),
   delta.sev.2 = c(1, 0.99, 0.90, 0.75, 0.60, 0.50, 0.25, 0.25, 0.50, 0.60, 0.75, 0.90, 0.99, 1.00))}

   \ifelse{text}{}{The figure below shows the corresponding selection functions.}

   \if{html}{\figure{selmodel-stepfun-fixed.png}{options: width=600}}
   \if{latex}{\figure{selmodel-stepfun-fixed.pdf}{options: width=4in}}

   These four functions are \dQuote{merely examples and should not be regarded as canonical} (Vevea & Woods, 2005).

   }

}
\value{
   An object of class \code{c("rma.uni","rma")}. The object is a list containing the same components as a regular \code{c("rma.uni","rma")} object, but the parameter estimates are based on the selection model. Most importantly, the following elements are modified based on the selection model:
   \item{beta}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{tau2}{estimated amount of (residual) heterogeneity. Always \code{0} when \code{method="FE"}.}
   \item{se.tau2}{standard error of the estimated amount of (residual) heterogeneity.}

   In addition, the object contains the following additional elements:
   \item{delta}{estimated selection model parameter(s).}
   \item{se.delta}{corresponding standard error(s).}
   \item{zval.delta}{corresponding test statistic(s).}
   \item{pval.delta}{corresponding p-value(s).}
   \item{ci.lb.delta}{lower bound of the confidence intervals for the parameter(s).}
   \item{ci.ub.delta}{upper bound of the confidence intervals for the parameter(s).}
   \item{LRT}{test statistic of the likelihood ratio test for the selection model parameter(s).}
   \item{LRTdf}{degrees of freedom for the likelihood ratio test.}
   \item{LRTp}{p-value for the likelihood ratio test.}
   \item{LRT.tau2}{test statistic of the likelihood ratio test for testing \mjeqn{\mbox{H}_0{:}\; \tau^2 = 0}{H_0: \tau^2 = 0} (\code{NA} when fitting a fixed-effects model).}
   \item{LRTp.tau2}{p-value for the likelihood ratio test.}
   \item{\dots}{some additional elements/values.}
}
\section{Methods}{
   The results of the fitted model are formatted and printed with the \code{\link{print.rma.uni}} function. The estimated selection function can be drawn with \code{\link{plot.rma.uni.selmodel}}.

   The \code{\link{profile.rma.uni.selmodel}} function can be used to obtain a plot of the log-likelihood as a function of \mjseqn{\tau^2} and/or the selection model parameter(s) of the model. Corresponding confidence intervals can be obtained with the \code{\link{confint.rma.uni.selmodel}} function.
}
\note{
   Model fitting is done via numerical optimization over the model parameters. By default, \code{\link{optim}} is used for the optimization. One can also chose a different optimizer via the \code{control} argument (e.g., \code{control=list(optimizer="nlminb")}). When using \code{\link{optim}}, one can set the particular method via the \code{optmethod} argument (e.g., \code{control=list(optimizer="optim", optmethod="BFGS")}, which is the default). Besides \code{\link{optim}} and \code{\link{nlminb}}, one can also choose one of the optimizers from the \code{minqa} package (i.e., \code{\link[minqa]{uobyqa}}, \code{\link[minqa]{newuoa}}, or \code{\link[minqa]{bobyqa}}), one of the (derivative-free) algorithms from the \code{\link[nloptr]{nloptr}} package, the Newton-type algorithm implemented in \code{\link{nlm}}, the various algorithms implemented in the \code{dfoptim} package (\code{\link[dfoptim]{hjk}} for the Hooke-Jeeves, \code{\link[dfoptim]{nmk}} for the Nelder-Mead, and \code{\link[dfoptim]{mads}} for the Mesh Adaptive Direct Searches (MADS) algorithm), the quasi-Newton type optimizer \code{\link[ucminf]{ucminf}} from the package of the same name, or the parallelized version of the L-BFGS-B algorithm implemented in \code{\link[optimParallel]{optimParallel}} from the package of the same name.

   The optimizer name must be given as a character string (i.e., in quotes). Additional control parameters can be specified via the \code{control} argument (e.g., \code{control=list(maxit=1000, reltol=1e-8)}). For \code{\link[nloptr]{nloptr}}, the default is to use the BOBYQA implementation from that package with a relative convergence criterion of \code{1e-8} on the function value (i.e., log-likelihood), but this can be changed via the \code{algorithm} and \code{ftop_rel} arguments (e.g., \code{control=list(optimizer="nloptr", algorithm="NLOPT_LN_SBPLX", ftol_rel=1e-6)}). For \code{\link[optimParallel]{optimParallel}}, the control argument \code{ncpus} can be used to specify the number of cores to use for the parallelization (e.g., \code{control=list(optimizer="optimParallel", ncpus=2)}). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.

   All selection models (except for \code{type="stepfun"}) require repeated evaluations of an integral, which is done via adaptive quadrature as implemented in the \code{\link{integrate}} function. One can adjust the arguments of the \code{integrate} function via control element \code{intCtrl}, which is a list of named arguments (e.g., \code{control = list(intCtrl = list(rel.tol=1e-4, subdivisions=100))}).

   The starting values for the fixed effects, the \mjseqn{\tau^2} value (only relevant in random/mixed-effects models), and the \mjseqn{\delta} parameter(s) are chosen automatically by the function, but one can also set the starting values manually via the \code{control} argument by specifying a vector of the appropriate length for \code{beta.init}, a single value for \code{tau2.init}, and a vector of the appropriate length for \code{delta.init}.

   By default, the \mjseqn{\delta} parameter(s) are constrained to a certain range, which improves the stability of the optimization algorithm. For all models, the maximum is set to \code{100} and the minimum to \code{0} (except for \code{type="beta"}, where the minimum for both parameters is \code{1e-05}). These defaults can be changed via the \code{control} argument by specifying a vector of the appropriate length for \code{delta.min} and/or \code{delta.max}.

   A difficulty with fitting the beta selection model (i.e., \code{type="beta"}) is the behavior of \mjseqn{w(p_i)} when \mjseqn{p_i = 0} or \mjseqn{p_i = 1}. When \mjseqn{\delta_1 < 1} or \mjseqn{\delta_2 < 1}, then this leads to selection weights equal to infinity, which causes problems when computing the likelihood function. Following Citkowicz and Vevea (2017), this problem can be avoided by censoring p-values too close to 0 or 1. The specific censoring point can be set via the \code{pval.min} element of the \code{control} argument. The default for this selection model is \code{control=list(pval.min=1e-5)}. A similar issues arises for the power selection model (i.e., \code{type="power"}) when \mjseqn{p_i = 1}. Again, \code{pval.min=1e-5} is used to circumvent this issue. For all other selection models, the default is \code{pval.min=0}.

   The variance-covariance matrix corresponding to the estimates of the fixed effects, the \mjseqn{\tau^2} value (only relevant in random/mixed-effects models), and the \mjseqn{\delta} parameter(s) is obtained by inverting the Hessian, which is numerically approximated using the \code{\link[numDeriv]{hessian}} function from the \code{numDeriv} package. This may fail, leading to \code{NA} values for the standard errors and hence test statistics, p-values, and confidence interval bounds. One can set control argument \code{hessianCtrl} to a list of named arguments to be passed on to the \code{method.args} argument of the \code{\link[numDeriv]{hessian}} function (the default is \code{control=list(hessianCtrl=list(r=6))}).

   Information on the progress of the optimization algorithm can be obtained by setting \code{verbose=TRUE} (this won't work when using parallelization). This option is useful to determine how the model fitting is progressing. One can also set \code{verbose} to an integer (\code{verbose=2} yields even more information and \code{verbose=3} also show the progress visually by drawing the selection function as the optimization proceeds).

   For selection functions where the \code{prec} argument is relevant, using (a function of) the sample sizes as the measure of precision (i.e., \code{prec="ninv"} or \code{prec="sqrtninv"}) is only possible when information about the sample sizes is actually stored within the object passed to the \code{selmodel} function. That should automatically be the case when the observed effect sizes or outcomes were computed with the \code{\link{escalc}} function or when the observed effect sizes or outcomes were computed within the model fitting function. On the other hand, this will not be the case when \code{\link{rma.uni}} was used together with the \code{yi} and \code{vi} arguments and the \code{yi} and \code{vi} values were \emph{not} computed with \code{\link{escalc}}. In that case, it is still possible to pass information about the sample sizes to the \code{\link{rma.uni}} function (e.g., use \code{rma.uni(yi, vi, ni=ni, data=dat)}, where data frame \code{dat} includes a variable called \code{ni} with the sample sizes).

   Finally, the automatic rescaling of the chosen precision measure can be switched off by setting \code{scaleprec=FALSE}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Begg, C. B., & Mazumdar, M. (1994). Operating characteristics of a rank correlation test for publication bias. \emph{Biometrics}, \bold{50}(4), 1088--1101. \verb{https://doi.org/10.2307/2533446}

   Carter, E. C., \enc{Schönbrodt}{Schoenbrodt}, F. D., Gervais, W. M., & Hilgard, J. (2019). Correcting for bias in psychology: A comparison of meta-analytic methods. \emph{Advances in Methods and Practices in Psychological Science}, \bold{2}(2), 115--144. \verb{https://doi.org/10.1177/2515245919847196}

   Citkowicz, M., & Vevea, J. L. (2017). A parsimonious weight function for modeling publication bias. \emph{Psychological Methods}, \bold{22}(1), 28--41. \verb{https://doi.org/10.1037/met0000119}

   Hedges, L. V. (1992). Modeling publication selection effects in meta-analysis. \emph{Statistical Science}, \bold{7}(2), 246--255. \verb{https://doi.org/10.1214/ss/1177011364}

   Iyengar, S., & Greenhouse, J. B. (1988). Selection models and the file drawer problem. \emph{Statistical Science}, \bold{3}(1), 109--117. \verb{https://doi.org/10.1214/ss/1177013012}

   McShane, B. B., Bockenholt, U., & Hansen, K. T. (2016). Adjusting for publication bias in meta-analysis: An evaluation of selection methods and some cautionary notes. \emph{Perspectives on Psychological Science}, \bold{11}(5), 730--749. \verb{https://doi.org/10.1177/1745691616662243}

   Preston, C., Ashby, D., & Smyth, R. (2004). Adjusting for publication bias: Modelling the selection process. \emph{Journal of Evaluation in Clinical Practice}, \bold{10}(2), 313--322. \verb{https://doi.org/10.1111/j.1365-2753.2003.00457.x}

   Pustejovsky, J. E., & Rodgers, M. A. (2019). Testing for funnel plot asymmetry of standardized mean differences. \emph{Research Synthesis Methods}, \bold{10}(1), 57--71. \verb{https://doi.org/10.1002/jrsm.1332}

   Vevea, J. L., & Hedges, L. V. (1995). A general linear model for estimating effect size in the presence of publication bias. \emph{Psychometrika}, \bold{60}(3), 419--435. \verb{https://doi.org/10.1007/BF02294384}

   Vevea, J. L., & Woods, C. M. (2005). Publication bias in research synthesis: Sensitivity analysis using a priori weight functions. \emph{Psychological Methods}, \bold{10}(4), 428--443. \verb{https://doi.org/10.1037/1082-989X.10.4.428}
}
\seealso{
   \code{\link{rma.uni}}
}
\examples{
############################################################################

### example from Citkowicz and Vevea (2017) for beta selection model

# copy data into 'dat' and examine data
dat <- dat.baskerville2012
dat

# fit random-effects model
res <- rma(smd, se^2, data=dat, method="ML", digits=3)
res

# funnel plot
funnel(res, ylim=c(0,0.6), xlab="Standardized Mean Difference")

# fit beta selection model
\dontrun{
sel <- selmodel(res, type="beta")
sel

# plot selection function
plot(sel, ylim=c(0,40))}

# fit mixed-effects meta-regression model with 'blind' dummy variable as moderator
res <- rma(smd, se^2, data=dat, mods = ~ blind, method="ML", digits=3)
res

# predicted average effect for studies that do not and that do use blinding
predict(res, newmods=c(0,1))

# fit beta selection model
\dontrun{
sel <- selmodel(res, type="beta")
sel
predict(sel, newmods=c(0,1))}

############################################################################

### example from Preston et al. (2004)

# copy data into 'dat' and examine data
dat <- dat.hahn2001
dat

### meta-analysis of (log) odds rations using the Mantel-Haenszel method
res <- rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, digits=2, slab=study)
res

# calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, drop00=TRUE)
dat

# fit fixed-effects model
res <- rma(yi, vi, data=dat, method="FE")

# predicted odds ratio (with 95\% CI)
predict(res, transf=exp, digits=2)

# funnel plot
funnel(res, atransf=exp, at=log(c(0.01,0.1,1,10,100)), ylim=c(0,2))

# fit half-normal, negative-exponential, logistic, and power selection models
\dontrun{
sel1 <- selmodel(res, type="halfnorm", alternative="less")
sel2 <- selmodel(res, type="negexp",   alternative="less")
sel3 <- selmodel(res, type="logistic", alternative="less")
sel4 <- selmodel(res, type="power",    alternative="less")

# plot selection functions
plot(sel1)
plot(sel2, add=TRUE, col="blue")
plot(sel3, add=TRUE, col="red")
plot(sel4, add=TRUE, col="green")

# show estimates of delta (and corresponding SEs)
tab <- data.frame(delta = c(sel1$delta, sel2$delta, sel3$delta, sel4$delta),
                  se    = c(sel1$se.delta, sel2$se.delta, sel3$se.delta, sel4$se.delta))
rownames(tab) <- c("Half-normal", "Negative-exponential", "Logistic", "Power")
round(tab, 2)

# predicted odds ratios (with 95\% CI)
predict(res,  transf=exp, digits=2)
predict(sel1, transf=exp, digits=2)
predict(sel2, transf=exp, digits=2)
predict(sel3, transf=exp, digits=2)
predict(sel4, transf=exp, digits=2)}

# fit selection models including standard error as precision measure (note: using
# scaleprec=FALSE here since Preston et al. (2004) did not use the rescaling)
\dontrun{
sel1 <- selmodel(res, type="halfnorm", prec="sei", alternative="less", scaleprec=FALSE)
sel2 <- selmodel(res, type="negexp",   prec="sei", alternative="less", scaleprec=FALSE)
sel3 <- selmodel(res, type="logistic", prec="sei", alternative="less", scaleprec=FALSE)
sel4 <- selmodel(res, type="power",    prec="sei", alternative="less", scaleprec=FALSE)

# show estimates of delta (and corresponding SEs)
tab <- data.frame(delta = c(sel1$delta, sel2$delta, sel3$delta, sel4$delta),
                  se    = c(sel1$se.delta, sel2$se.delta, sel3$se.delta, sel4$se.delta))
rownames(tab) <- c("Half-normal", "Negative-exponential", "Logistic", "Power")
round(tab, 2)

# predicted odds ratio (with 95\% CI)
predict(res,  transf=exp, digits=2)
predict(sel1, transf=exp, digits=2)
predict(sel2, transf=exp, digits=2)
predict(sel3, transf=exp, digits=2)
predict(sel4, transf=exp, digits=2)}

############################################################################

### meta-analysis on the effect of environmental tobacco smoke on lung cancer risk

# copy data into 'dat' and examine data
dat <- dat.hackshaw1998
dat

# fit random-effects model
res <- rma(yi, vi, data=dat, method="ML")
res

# funnel plot
funnel(res, atransf=exp, at=log(c(0.25,0.5,1,2,4,8)), ylim=c(0,0.8))

# step function selection model
\dontrun{
sel <- selmodel(res, type="stepfun", alternative="greater", steps=c(.025,.10,.50,1))
sel

# plot selection function
plot(sel)}

############################################################################

### validity of student ratings example from Vevea & Woods (2005)

# copy data into 'dat' and examine data
dat <- dat.cohen1981
dat

# calculate r-to-z transformed correlations and corresponding sampling variances
dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat[c(1,4,5)])
dat

# fit random-effects model
res <- rma(yi, vi, data=dat, method="ML", digits=3)
res

# predicted average correlation (with 95\% CI)
predict(res, transf=transf.ztor)

# funnel plot
funnel(res, ylim=c(0,0.4))

# selection functions from Vevea & Woods (2005)
tab <- data.frame(
   steps = c(0.005, 0.01, 0.05, 0.10, 0.25, 0.35, 0.50, 0.65, 0.75, 0.90, 0.95, 0.99, 0.995, 1),
   delta.mod.1 = c(1, 0.99, 0.95, 0.80, 0.75, 0.65, 0.60, 0.55, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50),
   delta.sev.1 = c(1, 0.99, 0.90, 0.75, 0.60, 0.50, 0.40, 0.35, 0.30, 0.25, 0.10, 0.10, 0.10, 0.10),
   delta.mod.2 = c(1, 0.99, 0.95, 0.90, 0.80, 0.75, 0.60, 0.60, 0.75, 0.80, 0.90, 0.95, 0.99, 1.00),
   delta.sev.2 = c(1, 0.99, 0.90, 0.75, 0.60, 0.50, 0.25, 0.25, 0.50, 0.60, 0.75, 0.90, 0.99, 1.00))

# apply step function model with a priori chosen selection weights
\dontrun{
sel <- lapply(tab[-1], function(delta) selmodel(res, type="stepfun", steps=tab$steps, delta=delta))

# estimates (transformed correlation) and tau^2 values
sav <- data.frame(estimate = round(c(res$beta, sapply(sel, function(x) x$beta)), 2),
                  varcomp  = round(c(res$tau2, sapply(sel, function(x) x$tau2)), 3))
sav}

############################################################################
}
\keyword{models}
