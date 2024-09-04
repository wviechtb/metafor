\name{rma.uni}
\alias{rma.uni}
\alias{rma}
\title{Meta-Analysis via Linear (Mixed-Effects) Models}
\description{
   Function to fit meta-analytic equal-, fixed-, and random-effects models and (mixed-effects) meta-regression models using a linear (mixed-effects) model framework. See below and the introduction to the \pkg{\link{metafor-package}} for more details on these models. \loadmathjax
}
\usage{
rma.uni(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
        m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, fi, pi, sdi, r2i, ni, mods, scale,
        measure="GEN", intercept=TRUE, data, slab, subset,
        add=1/2, to="only0", drop00=FALSE, vtype="LS",
        method="REML", weighted=TRUE, test="z",
        level=95, btt, att, tau2, verbose=FALSE, digits, control, \dots)
rma(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
        m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, fi, pi, sdi, r2i, ni, mods, scale,
        measure="GEN", intercept=TRUE, data, slab, subset,
        add=1/2, to="only0", drop00=FALSE, vtype="LS",
        method="REML", weighted=TRUE, test="z",
        level=95, btt, att, tau2, verbose=FALSE, digits, control, \dots)
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
   \item{fi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{pi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{sdi}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{r2i}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ni}{see below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \mjseqn{k} rows and as many columns as there are moderator variables. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{scale}{optional argument to include one or more predictors for the scale part in a location-scale model. See \sQuote{Details}.}
   \item{measure}{character string to specify the type of data supplied to the function. When \code{measure="GEN"} (default), the observed effect sizes or outcomes and corresponding sampling variances should be supplied to the function via the \code{yi} and \code{vi} arguments, respectively (instead of the sampling variances, one can supply the standard errors via the \code{sei} argument). Alternatively, one can set \code{measure} to one of the effect sizes or outcome measures described under the documentation for the \code{\link{escalc}} function in which case one must specify the required data via the appropriate arguments (see \code{\link{escalc}}).}
   \item{intercept}{logical to specify whether an intercept should be added to the model (the default is \code{TRUE}). Ignored when \code{mods} is a formula.}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \mjseqn{k} studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the analysis.}
   \item{add}{see the documentation of the \code{\link{escalc}} function.}
   \item{to}{see the documentation of the \code{\link{escalc}} function.}
   \item{drop00}{see the documentation of the \code{\link{escalc}} function.}
   \item{vtype}{see the documentation of the \code{\link{escalc}} function.}
   \item{method}{character string to specify whether an equal- or a random-effects model should be fitted. An equal-effects model is fitted when using \code{method="EE"}. A random-effects model is fitted by setting \code{method} equal to one of the following: \code{"DL"}, \code{"HE"}, \code{"HS"}, \code{"HSk"}, \code{"SJ"}, \code{"ML"}, \code{"REML"}, \code{"EB"}, \code{"PM"}, \code{"GENQ"}, \code{"PMM"}, or \code{"GENQM"}. The default is \code{"REML"}. See \sQuote{Details}.}
   \item{weighted}{logical to specify whether weighted (default) or unweighted estimation should be used to fit the model (the default is \code{TRUE}).}
   \item{test}{character string to specify how test statistics and confidence intervals for the fixed effects should be computed. By default (\code{test="z"}), Wald-type tests and CIs are obtained, which are based on a standard normal distribution. When \code{test="t"}, a t-distribution is used instead. When \code{test="knha"}, the method by Knapp and Hartung (2003) is used. See \sQuote{Details} and also \link[=misc-recs]{here} for some recommended practices.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95; see \link[=misc-options]{here} for details).}
   \item{btt}{optional vector of indices to specify which coefficients to include in the omnibus test of moderators. Can also be a string to \code{\link{grep}} for. See \sQuote{Details}.}
   \item{att}{optional vector of indices to specify which scale coefficients to include in the omnibus test. Only relevant for location-scale models. See \sQuote{Details}.}
   \item{tau2}{optional numeric value to specify the amount of (residual) heterogeneity in a random- or mixed-effects model (instead of estimating it). Useful for sensitivity analyses (e.g., for plotting results as a function of \mjseqn{\tau^2}). If unspecified, the value of \mjseqn{\tau^2} is estimated from the data.}
   \item{verbose}{logical to specify whether output should be generated on the progress of the model fitting (the default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is 4. See also \link[=misc-options]{here} for further details on how to control the number of digits in the output.}
   \item{control}{optional list of control values for the iterative estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{additional arguments.}
}
\details{

   \subsection{Specifying the Data}{

      The function can be used in combination with any of the usual effect sizes or outcome measures used in meta-analyses (e.g., log risk ratios, log odds ratios, risk differences, mean differences, standardized mean differences, log transformed ratios of means, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation), or, more generally, any set of estimates (with corresponding sampling variances) one would like to analyze. Simply specify the observed effect sizes or outcomes via the \code{yi} argument and the corresponding sampling variances via the \code{vi} argument. Instead of specifying \code{vi}, one can specify the standard errors (the square root of the sampling variances) via the \code{sei} argument. The \code{\link{escalc}} function can be used to compute a wide variety of effect sizes or outcome measures (and the corresponding sampling variances) based on summary statistics.

      Alternatively, the function can automatically calculate the values of a chosen effect size or outcome measure (and the corresponding sampling variances) when supplied with the necessary data. The \code{\link{escalc}} function describes which effect sizes or outcome measures are currently implemented and what data/arguments should then be specified/used. The \code{measure} argument should then be set to the desired effect size or outcome measure.

   }

   \subsection{Specifying the Model}{

      The function can be used to fit equal-, fixed-, and random-effects models, as well as (mixed-effects) meta-regression models including one or multiple moderators (the difference between the various models is described in detail on the introductory \pkg{\link{metafor-package}} help page).

      Assuming the observed effect sizes or outcomes and corresponding sampling variances are supplied via the \code{yi} and \code{vi} arguments, an \emph{equal-effects model} can be fitted with \code{rma(yi, vi, method="EE")}. Setting \code{method="FE"} fits a \emph{fixed-effects model} (see \link[=misc-models]{here} for a discussion of this model and how the interpretation of these models differ despite yielding identical results). Weighted estimation (with inverse-variance weights) is used by default. User-defined weights can be supplied via the \code{weights} argument. Unweighted estimation can be used by setting \code{weighted=FALSE} (which is the same as setting the weights equal to a constant).

      A \emph{random-effects model} can be fitted with the same code but setting the \code{method} argument to one of the various estimators for the amount of heterogeneity:

      \itemize{
      \item \code{method="DL"} = DerSimonian-Laird estimator (DerSimonian & Laird, 1986; Raudenbush, 2009),
      \item \code{method="HE"} = Hedges estimator (Hedges, 1983, 1992),
      \item \code{method="HS"} = Hunter-Schmidt estimator (Hunter & Schmidt, 1990; Viechtbauer et al., 2015),
      \item \code{method="HSk"} = Hunter-Schmidt estimator with a small sample-size correction (Brannick et al., 2019),
      \item \code{method="SJ"} = Sidik-Jonkman estimator (Sidik & Jonkman, 2005b, 2007),
      \item \code{method="ML"} = maximum likelihood estimator (Hardy & Thompson, 1996; Raudenbush, 2009),
      \item \code{method="REML"} = restricted maximum likelihood estimator (Viechtbauer, 2005; Raudenbush, 2009)
      \item \code{method="EB"} = empirical Bayes estimator (Morris, 1983; Berkey et al. 1995),
      \item \code{method="PM"} = Paule-Mandel estimator (Paule & Mandel, 1982; Viechtbauer et al., 2015),
      \item \code{method="GENQ"} = generalized Q-statistic estimator (DerSimonian & Kacker, 2007; Jackson et al., 2014),
      \item \code{method="PMM"} = median-unbiased Paule-Mandel estimator (Viechtbauer, 2021),
      \item \code{method="GENQM"} = median-unbiased generalized Q-statistic estimator (Viechtbauer, 2021).
      }

      For a description of the various estimators, see Brannick et al. (2019), DerSimonian and Kacker (2007), Raudenbush (2009), Veroniki et al. (2016), Viechtbauer (2005), and Viechtbauer et al. (2015). Note that the Hedges estimator is also called the \sQuote{variance component estimator} or \sQuote{Cochran estimator}, the Sidik-Jonkman estimator is also called the \sQuote{model error variance estimator}, the empirical Bayes estimator is actually identical to the Paule-Mandel estimator (Viechtbauer et al., 2015), and the generalized Q-statistic estimator is a general method-of-moments estimator (DerSimonian & Kacker, 2007) requiring the specification of weights (the HE and DL estimators are just special cases with equal and inverse sampling variance weights, respectively). Finally, the two median-unbiased estimators are versions of the Paule-Mandel and generalized Q-statistic estimators that equate the respective estimating equations not to their expected values, but to the medians of their theoretical distributions (Viechtbauer, 2021).

      One or more moderators can be included in a model via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate model matrix (i.e., \mjseqn{X}) with \mjseqn{k} rows and as many columns as there are moderator variables (e.g., \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for three moderator variables). The intercept is added to the model matrix by default unless \code{intercept=FALSE}.

      Alternatively, one can use standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial/spline terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept is controlled by the specified formula (e.g., \code{mods = ~ 0 + mod1 + mod2 + mod3} or \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept).

      When the observed effect sizes or outcomes and corresponding sampling variances are supplied via the \code{yi} and \code{vi} (or \code{sei}) arguments, one can also specify moderators via the \code{yi} argument (e.g., \code{rma(yi ~ mod1 + mod2 + mod3, vi)}). In that case, the \code{mods} argument is ignored and the inclusion/exclusion of the intercept is again controlled by the specified formula.

   }

   \subsection{Omnibus Test of Moderators}{

      For models including moderators, an omnibus test of all model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} (\sQuote{betas to test}) argument (i.e., to test \mjseqn{\text{H}_0{:}\; \beta_{j \in \texttt{btt}} = 0}, where \mjseqn{\beta_{j \in \texttt{btt}}} is the set of coefficients to be tested). For example, with \code{btt=c(3,4)}, only the third and fourth coefficients from the model are included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model). Instead of specifying the coefficient numbers, one can specify a string for \code{btt}. In that case, \code{\link{grep}} will be used to search for all coefficient names that match the string. The omnibus test is called the \mjseqn{Q_M}-test and follows asymptotically a chi-square distribution with \mjseqn{m} degrees of freedom (with \mjseqn{m} denoting the number of coefficients tested) under the null hypothesis (that the true value of all coefficients tested is equal to 0).

   }

   \subsection{Categorical Moderators}{

      Categorical moderator variables can be included in the model via the \code{mods} argument in the same way that appropriately (dummy) coded categorical variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to automate the coding (note that string/character variables in a model formula are automatically converted to factors). An example to illustrate these different approaches is provided below.

   }

   \subsection{Tests and Confidence Intervals}{

      By default, tests of individual coefficients in the model (and the corresponding confidence intervals) are based on a standard normal distribution, while the omnibus test is based on a chi-square distribution (see above). As an alternative, one can set \code{test="t"}, in which case tests of individual coefficients and confidence intervals are based on a t-distribution with \mjseqn{k-p} degrees of freedom, while the omnibus test then uses an F-distribution with \mjseqn{m} and \mjseqn{k-p} degrees of freedom (with \mjseqn{k} denoting the total number of estimates included in the analysis and \mjseqn{p} the total number of model coefficients including the intercept if it is present). Furthermore, when \code{test="knha"} (or equivalently, \code{test="hksj"}), the method by Hartung (1999), Sidik and Jonkman (2002), and Knapp and Hartung (2003) (the Knapp-Hartung method; also referred to as the Hartung-Knapp-Sidik-Jonkman method) is used, which applies an adjustment to the standard errors of the estimated coefficients (to account for the uncertainty in the estimate of the amount of (residual) heterogeneity) and uses t- and F-distributions as described above (see also \link[=misc-recs]{here}). Finally, one can set \code{test="adhoc"}, in which case the Knapp-Hartung method is used, but with the restriction that the adjustment to the standard errors can never result in adjusted standard errors that are smaller than the unadjusted ones (see Jackson et al., 2017, section 4.3).

   }

   \subsection{Test for (Residual) Heterogeneity}{

      A test for (residual) heterogeneity is automatically carried out by the function. Without moderators in the model, this is simply Cochran's \mjseqn{Q}-test (Cochran, 1954), which tests whether the variability in the observed effect sizes or outcomes is larger than would be expected based on sampling variability alone. A significant test suggests that the true effects/outcomes are heterogeneous. When moderators are included in the model, this is the \mjseqn{Q_E}-test for residual heterogeneity, which tests whether the variability in the observed effect sizes or outcomes not accounted for by the moderators included in the model is larger than would be expected based on sampling variability alone.

   }

   \subsection{Location-Scale Models}{

      The function can also be used to fit so-called \sQuote{location-scale models} (Viechtbauer & \enc{López-López}{Lopez-Lopez}, 2022). In such models, one can specify not only predictors for the size of the average true outcome (i.e., for their \sQuote{location}), but also predictors for the amount of heterogeneity in the outcomes (i.e., for their \sQuote{scale}). The model is given by \mjdeqn{y_i = \beta_0 + \beta_1 x_{i1} + \beta_2 x_{i2} + \ldots + \beta_{p'} x_{ip'} + u_i + \varepsilon_i,}{y_i = \beta_0 + \beta_1 x_i1 + \beta_2 x_i2 + \ldots + \beta_p' x_ip' + u_i + \epsilon_i,} \mjdeqn{u_i \sim N(0, \tau_i^2), \; \varepsilon_i \sim N(0, v_i),}{u_i ~ N(0, tau_i^2), \epsilon_i \sim N(0, v_i),} \mjdeqn{\log(\tau_i^2) = \alpha_0 + \alpha_1 z_{i1} + \alpha_2 z_{i2} + \ldots + \alpha_{q'} z_{iq'},}{log(tau^2) = \alpha_0 + \alpha z_i1 + \alpha z_i2 + \ldots + \alpha_q' z_iq',} where \mjeqn{x_{i1}, \ldots, x_{ip'}}{x_i1, \ldots, x_ip'} are the values of the \mjseqn{p'} predictor variables that may be related to the size of the average true outcome (letting \mjseqn{p = p' + 1} denote the total number of location coefficients in the model including the model intercept \mjseqn{\beta_0}) and \mjeqn{z_{i1}, \ldots, z_{iq'}}{z_i1, \ldots, z_iq'} are the values of the \mjseqn{q'} scale variables that may be related to the amount of heterogeneity in the outcomes (letting \mjseqn{q = q' + 1} denote the total number of scale coefficients in the model including the model intercept \mjseqn{\alpha_0}). Location variables can be specified via the \code{mods} argument as described above (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Scale variables can be specified via the \code{scale} argument (e.g., \code{scale = ~ var1 + var2 + var3}). A log link is used for specifying the relationship between the scale variables and the amount of heterogeneity so that \mjseqn{\tau_i^2} is guaranteed to be non-negative (one can also set (the undocumented) argument \code{link="identity"} to use an identity link, but this is more likely to lead to estimation problems). Estimates of the location and scale coefficients can be obtained either with maximum likelihood (\code{method="ML"}) or restricted maximum likelihood (\code{method="REML"}) estimation. An omnibus test of the scale coefficients is conducted as described above (where the \code{att} argument can be used to specify which scale coefficients to include in the test).

   }

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
   \item{tau2}{estimated amount of (residual) heterogeneity. Always \code{0} when \code{method="EE"}.}
   \item{se.tau2}{standard error of the estimated amount of (residual) heterogeneity.}
   \item{k}{number of studies included in the analysis.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of moderators.}
   \item{QE}{test statistic of the test for (residual) heterogeneity.}
   \item{QEp}{corresponding p-value.}
   \item{QM}{test statistic of the omnibus test of moderators.}
   \item{QMp}{corresponding p-value.}
   \item{I2}{value of \mjseqn{I^2}. See \code{\link[=print.rma.uni]{print}} for more details.}
   \item{H2}{value of \mjseqn{H^2}. See \code{\link[=print.rma.uni]{print}} for more details.}
   \item{R2}{value of \mjseqn{R^2}. See \code{\link[=print.rma.uni]{print}} for more details.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, vi, X}{the vector of outcomes, the corresponding sampling variances, and the model matrix.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values under the unrestricted and restricted likelihood.}
   \item{\dots}{some additional elements/values.}

   For location-scale models, the object is of class \code{c("rma.ls","rma.uni","rma")} and includes the following components in addition to the ones listed above:
   \item{alpha}{estimated scale coefficients of the model.}
   \item{se.alpha}{standard errors of the coefficients.}
   \item{zval.alpha}{test statistics of the coefficients.}
   \item{pval.alpha}{corresponding p-values.}
   \item{ci.lb.alpha}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub.alpha}{upper bound of the confidence intervals for the coefficients.}
   \item{va}{variance-covariance matrix of the estimated coefficients.}
   \item{tau2}{as above, but now a vector of values.}
   \item{q}{number of scale coefficients in the model (including the intercept).}
   \item{QS}{test statistic of the omnibus test of the scale coefficients.}
   \item{QSp}{corresponding p-value.}
   \item{\dots}{some additional elements/values.}
}
\section{Methods}{
   The results of the fitted model are formatted and printed with the \code{\link[=print.rma.uni]{print}} function. If fit statistics should also be given, use \code{\link[=summary.rma]{summary}} (or use the \code{\link[=fitstats.rma]{fitstats}} function to extract them). Full versus reduced model comparisons in terms of fit statistics and likelihood ratio tests can be obtained with \code{\link[=anova.rma]{anova}}. Wald-type tests for sets of model coefficients or linear combinations thereof can be obtained with the same function. Permutation tests for the model coefficient(s) can be obtained with \code{\link[=permutest.rma.uni]{permutest}}. Tests and confidence intervals based on (cluster) robust methods can be obtained with \code{\link[=robust.rma.uni]{robust}}.

   Predicted/fitted values can be obtained with \code{\link[=predict.rma]{predict}} and \code{\link[=fitted.rma]{fitted}}. For best linear unbiased predictions, see \code{\link[=blup.rma.uni]{blup}} and \code{\link[=ranef.rma.uni]{ranef}}.

   The \code{\link[=residuals.rma]{residuals}}, \code{\link[=rstandard.rma.uni]{rstandard}}, and \code{\link[=rstudent.rma.uni]{rstudent}} functions extract raw and standardized residuals. Additional model diagnostics (e.g., to determine influential studies) can be obtained with the \code{\link[=influence.rma.uni]{influence}} function. For models without moderators, leave-one-out diagnostics can also be obtained with \code{\link[=leave1out.rma.uni]{leave1out}}. For models with moderators, variance inflation factors can be obtained with \code{\link[=vif.rma]{vif}}.

   A confidence interval for the amount of (residual) heterogeneity in the random/mixed-effects model can be obtained with \code{\link[=confint.rma.uni]{confint}}. For location-scale models, \code{\link[=confint.rma.ls]{confint}} can provide confidence intervals for the scale coefficients.

   Forest, funnel, radial, \enc{L'Abbé}{L'Abbe}, and Baujat plots can be obtained with \code{\link[=forest.rma]{forest}}, \code{\link[=funnel.rma]{funnel}}, \code{\link[=radial.rma]{radial}}, \code{\link[=labbe.rma]{labbe}}, and \code{\link[=baujat.rma]{baujat}} (radial and \enc{L'Abbé}{L'Abbe} plots only for models without moderators). The \code{\link[=qqnorm.rma.uni]{qqnorm}} function provides normal QQ plots of the standardized residuals. One can also just call \code{\link[=plot.rma.uni]{plot}} on the fitted model object to obtain various plots at once. For random/mixed-effects models, the \code{\link[=profile.rma.uni]{profile}} function can be used to obtain a plot of the (restricted) log-likelihood as a function of \mjseqn{\tau^2}. For location-scale models, \code{\link[=profile.rma.ls]{profile}} draws analogous plots based on the scale coefficients. For models with moderators, \code{\link[=regplot.rma]{regplot}} draws scatter plots / bubble plots, showing the (marginal) relationship between the observed outcomes and a selected moderator from the model.

   Tests for funnel plot asymmetry (which may be indicative of publication bias) can be obtained with \code{\link{ranktest}} and \code{\link{regtest}}. For models without moderators, the \code{\link[=trimfill.rma.uni]{trimfill}} method can be used to carry out a trim and fill analysis and \code{\link[=hc.rma.uni]{hc}} provides a random-effects model analysis that is more robust to publication bias (based on the method by Henmi & Copas, 2010). The test of \sQuote{excess significance} can be carried out with the \code{\link{tes}} function. The fail-safe N (based on a file drawer analysis) can be computed using \code{\link{fsn}}. Selection models can be fitted with the \code{\link{selmodel}} function.

   For models without moderators, a cumulative meta-analysis (i.e., adding one observation at a time) can be obtained with \code{\link[=cumul.rma.uni]{cumul}}.

   Other extractor functions include \code{\link[=coef.rma]{coef}}, \code{\link[=vcov.rma]{vcov}}, \code{\link[=logLik.rma]{logLik}}, \code{\link[=deviance.rma]{deviance}}, \code{\link[=AIC.rma]{AIC}}, \code{\link[=BIC.rma]{BIC}}, \code{\link[=hatvalues.rma.uni]{hatvalues}}, and \code{\link[=weights.rma.uni]{weights}}.
}
\note{
   While the HS, HSk, HE, DL, SJ, and GENQ estimators of \mjseqn{\tau^2} are based on closed-form solutions, the ML, REML, and EB estimators must be obtained iteratively. For this, the function makes use of the Fisher scoring algorithm, which is robust to poor starting values and usually converges quickly (Harville, 1977; Jennrich & Sampson, 1976). By default, the starting value is set equal to the value of the Hedges (HE) estimator and the algorithm terminates when the change in the estimated value of \mjseqn{\tau^2} is smaller than \mjeqn{10^{-5}}{10^(-5)} from one iteration to the next. The maximum number of iterations is 100 by default (which should be sufficient in most cases). Information on the progress of the algorithm can be obtained by setting \code{verbose=TRUE}. One can also set \code{verbose} to an integer (\code{verbose=2} yields even more information and \code{verbose=3} also sets \code{option(warn=1)} temporarily).

   A different starting value, threshold, and maximum number of iterations can be specified via the \code{control} argument by setting \code{control=list(tau2.init=value, threshold=value, maxiter=value)}. The step length of the Fisher scoring algorithm can also be adjusted by a desired factor with \code{control=list(stepadj=value)} (values below 1 will reduce the step length). If using \code{verbose=TRUE} shows the estimate jumping around erratically (or cycling through a few values), decreasing the step length (and increasing the maximum number of iterations) can often help with convergence (e.g., \code{control=list(stepadj=0.5, maxiter=1000)}).

   The PM, PMM, and GENQM estimators also involve iterative algorithms, which make use of the \code{\link{uniroot}} function. By default, the desired accuracy (\code{tol}) is set equal to \code{.Machine$double.eps^0.25} and the maximum number of iterations (\code{maxiter}) to \code{100} (as above). The upper bound of the interval searched (\code{tau2.max}) is set to the larger of 100 and \code{10*mad(yi)^2} (i.e., 10 times the squared median absolute deviation of the observed effect sizes or outcomes computed with the \code{\link{mad}} function). These values can be adjusted with \code{control=list(tol=value, maxiter=value, tau2.max=value)}.

   All of the heterogeneity estimators except SJ can in principle yield negative estimates for the amount of (residual) heterogeneity. However, negative estimates of \mjseqn{\tau^2} are outside of the parameter space. For the HS, HSk, HE, DL, and GENQ estimators, negative estimates are therefore truncated to zero. For the ML, REML, and EB estimators, the Fisher scoring algorithm makes use of step halving (Jennrich & Sampson, 1976) to guarantee a non-negative estimate. Finally, for the PM, PMM, and GENQM estimators, the lower bound of the interval searched is set to zero by default. For those brave enough to step into risky territory, there is the option to set the lower bound for all these estimators to some other value besides zero (even a negative one) with \code{control=list(tau2.min=value)}, but the lowest value permitted is \code{-min(vi)} (to ensure that the marginal variances are always non-negative).

   The Hunter-Schmidt estimator for the amount of heterogeneity is defined in Hunter and Schmidt (1990) only in the context of the random-effects model when analyzing correlation coefficients. A general version of this estimator for random- and mixed-effects models not specific to any particular outcome measure is described in Viechtbauer (2005) and Viechtbauer et al. (2015) and is implemented here.

   The Sidik-Jonkman estimator starts with a crude estimate of \mjseqn{\tau^2}, which is then updated as described in Sidik and Jonkman (2005b, 2007). If, instead of the crude estimate, one wants to use a better a priori estimate, one can do so by passing this value via \code{control=list(tau2.init=value)}.

   One can also specify a vector of estimators via the \code{method} argument (e.g., \code{rma(yi, vi, method=c("REML","DL"))}). The various estimators are then applied in turn until one converges. This is mostly useful for simulation studies where an estimator (like the REML estimator) is not guaranteed to converge and one can then substitute one (like the DL estimator) that does not involve iterative methods and is guaranteed to provide an estimate.

   Outcomes with non-positive sampling variances are problematic. If a sampling variance is equal to zero, then its weight will be \mjseqn{1/0} for equal-effects models when using weighted estimation. Switching to unweighted estimation is a possible solution then. For random/mixed-effects model, some estimators of \mjseqn{\tau^2} are undefined when there is at least one sampling variance equal to zero. Other estimators may work, but it may still be necessary to switch to unweighted model fitting, especially when the estimate of \mjseqn{\tau^2} converges to zero.

   When including moderators in the model, it is possible that the model matrix is not of full rank (i.e., there is a linear relationship between the moderator variables included in the model). The function automatically tries to reduce the model matrix to full rank by removing redundant predictors, but if this fails the model cannot be fitted and an error will be issued. Deleting (redundant) moderator variables from the model as needed should solve this problem.

   Some general words of caution about the assumptions underlying the models:

   \itemize{
   \item The sampling variances (i.e., the \code{vi} values) are treated as if they are known constants, even though in practice they are usually estimates themselves. This implies that the distributions of the test statistics and corresponding confidence intervals are only exact and have nominal coverage when the within-study sample sizes are large (i.e., when the error in the sampling variance estimates is small). Certain outcome measures (e.g., the arcsine square root transformed risk difference and Fisher's r-to-z transformed correlation coefficient) are based on variance stabilizing transformations that also help to make the assumption of known sampling variances much more reasonable.
   \item When fitting a mixed/random-effects model, \mjseqn{\tau^2} is estimated and then treated as a known constant thereafter. This ignores the uncertainty in the estimate of \mjseqn{\tau^2}. As a consequence, the standard errors of the parameter estimates tend to be too small, yielding test statistics that are too large and confidence intervals that are not wide enough. The Knapp and Hartung (2003) adjustment (i.e., using \code{test="knha"}) can be used to counter this problem, yielding test statistics and confidence intervals whose properties are closer to nominal.
   \item Most effect sizes or outcome measures do not have exactly normal sampling distributions as assumed under the various models. However, the normal approximation usually becomes more accurate for most effect sizes or outcome measures as the within-study sample sizes increase. Therefore, sufficiently large within-study sample sizes are (usually) needed to be certain that the tests and confidence intervals have nominal levels/coverage. Again, certain outcome measures (e.g., Fisher's r-to-z transformed correlation coefficient) may be preferable from this perspective as well.
   }

   For location-scale models, model fitting is done via numerical optimization over the model parameters. By default, \code{\link{nlminb}} is used for the optimization. One can also chose a different optimizer from \code{\link{optim}} via the \code{control} argument (e.g., \code{control=list(optimizer="BFGS")} or \code{control=list(optimizer="Nelder-Mead")}). Besides \code{\link{nlminb}} and one of the methods from \code{\link{optim}}, one can also choose one of the optimizers from the \code{minqa} package (i.e., \code{\link[minqa]{uobyqa}}, \code{\link[minqa]{newuoa}}, or \code{\link[minqa]{bobyqa}}), one of the (derivative-free) algorithms from the \code{\link[nloptr]{nloptr}} package, the Newton-type algorithm implemented in \code{\link{nlm}}, the various algorithms implemented in the \code{dfoptim} package (\code{\link[dfoptim]{hjk}} for the Hooke-Jeeves, \code{\link[dfoptim]{nmk}} for the Nelder-Mead, and \code{\link[dfoptim]{mads}} for the Mesh Adaptive Direct Searches algorithm), the quasi-Newton type optimizers \code{\link[ucminf]{ucminf}} and \code{\link[lbfgsb3c]{lbfgsb3c}} and the subspace-searching simplex algorithm \code{\link[subplex]{subplex}} from the packages of the same name, the Barzilai-Borwein gradient decent method implemented in \code{\link[BB]{BBoptim}}, the \code{\link[optimx]{Rcgmin}} and \code{\link[optimx]{Rvmmin}} optimizers, or the parallelized version of the L-BFGS-B algorithm implemented in \code{\link[optimParallel]{optimParallel}} from the package of the same name. When using an identity link with \code{link="identity"}, constrained optimization (to ensure non-negative \mjseqn{\tau_i^2} values) as implemented in \code{\link{constrOptim}} is used by default. Alternative optimizers in this case are the \code{\link[Rsolnp]{solnp}} solver from the \code{Rsolnp} package, \code{\link[nloptr]{nloptr}}, or the augmented Lagrangian adaptive barrier minimization algorithm \code{\link[alabama]{constrOptim.nl}} from the \code{alabama} package.

   The optimizer name must be given as a character string (i.e., in quotes). Additional control parameters can be specified via the \code{control} argument (e.g., \code{control=list(iter.max=1000, rel.tol=1e-8)}). For \code{\link[nloptr]{nloptr}}, the default is to use the BOBYQA implementation from that package with a relative convergence criterion of \code{1e-8} on the function value (i.e., log-likelihood), but this can be changed via the \code{algorithm} and \code{ftop_rel} arguments (e.g., \code{control=list(optimizer="nloptr", algorithm="NLOPT_LN_SBPLX", ftol_rel=1e-6)}) (note: when using \code{optimizer="nloptr"} in combination with an identity link, the \code{"NLOPT_LN_COBYLA"} algorithm is automatically used, since it allows for inequality constraints). For \code{\link[optimParallel]{optimParallel}}, the control argument \code{ncpus} can be used to specify the number of cores to use for the parallelization (e.g., \code{control=list(optimizer="optimParallel", ncpus=2)}). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.

   Under certain circumstances (e.g., when the amount of heterogeneity is very small for certain combinations of values for the scale variables and scale coefficients), the values of the scale coefficients may try to drift towards minus or plus infinity, which can lead to problems with the optimization. One can impose constraints on the scale coefficients via \code{control=list(alpha.min=minval, alpha.max=maxval)} where \code{minval} and \code{maxval} are either scalars or vectors of the appropriate length.

   Finally, for location-scale models, the standard errors of the scale coefficients are obtained by inverting the Hessian, which is numerically approximated using the \code{\link[numDeriv]{hessian}} function from the \code{numDeriv} package. This may fail (especially when using an identity link), leading to \code{NA} values for the standard errors and hence test statistics, p-values, and confidence interval bounds. One can set control argument \code{hessianCtrl} to a list of named arguments to be passed on to the \code{method.args} argument of the \code{\link[numDeriv]{hessian}} function (the default is \code{control=list(hessianCtrl=list(r=8))}). One can also set \code{control=list(hesspack="pracma")} in which case the \code{\link[pracma]{hessian}} function from the \code{pracma} package is used instead for approximating the Hessian.

   Even if the Hessian can be approximated and inverted, the standard errors may be unreasonably large when the likelihood surface is very flat around the estimated scale coefficients. This is more likely to happen when \mjseqn{k} is small and when the amount of heterogeneity is very small under some conditions as defined by the scale coefficients/variables. Setting constraints on the scale coefficients as described above can also help to mitigate this issue.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Berkey, C. S., Hoaglin, D. C., Mosteller, F., & Colditz, G. A. (1995). A random-effects regression model for meta-analysis. \emph{Statistics in Medicine}, \bold{14}(4), 395--411. \verb{https://doi.org/10.1002/sim.4780140406}

   Brannick, M. T., Potter, S. M., Benitez, B., & Morris, S. B. (2019). Bias and precision of alternate estimators in meta-analysis: Benefits of blending Schmidt–Hunter and Hedges approaches. \emph{Organizational Research Methods}, \bold{22}(2), 490--514. \verb{https://doi.org/10.1177/1094428117741966}

   Cochran, W. G. (1954). The combination of estimates from different experiments. \emph{Biometrics}, \bold{10}(1), 101--129. \verb{https://doi.org/10.2307/3001666}

   DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials. \emph{Controlled Clinical Trials}, \bold{7}(3), 177--188. \verb{https://doi.org/10.1016/0197-2456(86)90046-2}

   DerSimonian, R., & Kacker, R. (2007). Random-effects model for meta-analysis of clinical trials: An update. \emph{Contemporary Clinical Trials}, \bold{28}(2), 105--114. \verb{https://doi.org/10.1016/j.cct.2006.04.004}

   Hardy, R. J. & Thompson, S. G. (1996). A likelihood approach to meta-analysis with random effects. \emph{Statistics in Medicine}, \bold{15}(6), 619--629. \verb{https://doi.org/10.1002/(SICI)1097-0258(19960330)15:6<619::AID-SIM188>3.0.CO;2-A}

   Hartung, J. (1999). An alternative method for meta-analysis. \emph{Biometrical Journal}, \bold{41}(8), 901--916. \verb{https://doi.org/10.1002/(SICI)1521-4036(199912)41:8<901::AID-BIMJ901>3.0.CO;2-W}

   Harville, D. A. (1977). Maximum likelihood approaches to variance component estimation and to related problems. \emph{Journal of the American Statistical Association}, \bold{72}(358), 320--338. \verb{https://doi.org/10.2307/2286796}

   Hedges, L. V. (1983). A random effects model for effect sizes. \emph{Psychological Bulletin}, \bold{93}(2), 388--395. \verb{https://doi.org/10.1037/0033-2909.93.2.388}

   Hedges, L. V. (1992). Meta-analysis. \emph{Journal of Educational Statistics}, \bold{17}(4), 279--296. \verb{https://doi.org/10.3102/10769986017004279}

   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Henmi, M., & Copas, J. B. (2010). Confidence intervals for random effects meta-analysis and robustness to publication bias. \emph{Statistics in Medicine}, \bold{29}(29), 2969--2983. \verb{https://doi.org/10.1002/sim.4029}

   Hunter, J. E., & Schmidt, F. L. (1990). \emph{Methods of meta-analysis: Correcting error and bias in research findings}. Thousand Oaks, CA: Sage.

   Jackson, D., Turner, R., Rhodes, K. & Viechtbauer, W. (2014). Methods for calculating confidence and credible intervals for the residual between-study variance in random effects meta-regression models. \emph{BMC Medical Research Methodology}, \bold{14}, 103. \verb{https://doi.org/10.1186/1471-2288-14-103}

   Jackson, D., Law, M., \enc{Rücker}{Ruecker}, G., & Schwarzer, G. (2017). The Hartung-Knapp modification for random-effects meta-analysis: A useful refinement but are there any residual concerns? \emph{Statistics in Medicine}, \bold{36}(25), 3923--3934. \verb{https://doi.org/10.1002/sim.7411}

   Jennrich, R. I., & Sampson, P. F. (1976). Newton-Raphson and related algorithms for maximum likelihood variance component estimation. \emph{Technometrics}, \bold{18}(1), 11--17. \verb{https://doi.org/10.2307/1267911}

   Knapp, G., & Hartung, J. (2003). Improved tests for a random effects meta-regression with a single covariate. \emph{Statistics in Medicine}, \bold{22}(17), 2693--2710. \verb{https://doi.org/10.1002/sim.1482}

   Morris, C. N. (1983). Parametric empirical Bayes inference: Theory and applications. \emph{Journal of the American Statistical Association}, \bold{78}(381), 47--55. \verb{https://doi.org/10.2307/2287098}

   Paule, R. C., & Mandel, J. (1982). Consensus values and weighting factors. \emph{Journal of Research of the National Bureau of Standards}, \bold{87}(5), 377--385. \verb{https://doi.org/10.6028/jres.087.022}

   Raudenbush, S. W. (2009). Analyzing effect sizes: Random effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295--315). New York: Russell Sage Foundation.

   Sidik, K. & Jonkman, J. N. (2002). A simple confidence interval for meta-analysis. \emph{Statistics in Medicine}, \bold{21}(21), 3153--3159. \verb{https://doi.org/10.1002/sim.1262}

   Sidik, K., & Jonkman, J. N. (2005a). A note on variance estimation in random effects meta-regression. \emph{Journal of Biopharmaceutical Statistics}, \bold{15}(5), 823--838. \verb{https://doi.org/10.1081/BIP-200067915}

   Sidik, K., & Jonkman, J. N. (2005b). Simple heterogeneity variance estimation for meta-analysis. \emph{Journal of the Royal Statistical Society, Series C}, \bold{54}(2), 367--384. \verb{https://doi.org/10.1111/j.1467-9876.2005.00489.x}

   Sidik, K., & Jonkman, J. N. (2007). A comparison of heterogeneity variance estimators in combining results of studies. \emph{Statistics in Medicine}, \bold{26}(9), 1964--1981. \verb{https://doi.org/10.1002/sim.2688}

   Veroniki, A. A., Jackson, D., Viechtbauer, W., Bender, R., Bowden, J., Knapp, G., Kuss, O., Higgins, J. P., Langan, D., & Salanti, G. (2016). Methods to estimate the between-study variance and its uncertainty in meta-analysis. \emph{Research Synthesis Methods}, \bold{7}(1), 55--79. \verb{https://doi.org/10.1002/jrsm.1164}

   Viechtbauer, W. (2005). Bias and efficiency of meta-analytic variance estimators in the random-effects model. \emph{Journal of Educational and Behavioral Statistics}, \bold{30}(3), 261--293. \verb{https://doi.org/10.3102/10769986030003261}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W. (2021). Median-unbiased estimators for the amount of heterogeneity in meta-analysis. \emph{European Congress of Methodology}, Valencia, Spain. \verb{https://www.wvbauer.com/lib/exe/fetch.php/talks:2021_viechtbauer_eam_median_tau2.pdf}

   Viechtbauer, W., & \enc{López-López}{Lopez-Lopez}, J. A. (2022). Location-scale models for meta-analysis. \emph{Research Synthesis Methods}. \bold{13}(6), 697--715. \verb{https://doi.org/10.1002/jrsm.1562}

   Viechtbauer, W., \enc{López-López}{Lopez-Lopez}, J. A., \enc{Sánchez-Meca}{Sanchez-Meca}, J., & \enc{Marín-Martínez}{Marin-Martinez}, F. (2015). A comparison of procedures to test for moderators in mixed-effects meta-regression models. \emph{Psychological Methods}, \bold{20}(3), 360--374. \verb{https://doi.org/10.1037/met0000023}
}
\seealso{
   \code{\link{rma.mh}}, \code{\link{rma.peto}}, \code{\link{rma.glmm}}, and \code{\link{rma.mv}} for other model fitting functions.
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit a random-effects model using the log risk ratios and sampling variances as input
### note: method="REML" is the default, so one could leave this out
rma(yi, vi, data=dat, method="REML")

### fit a random-effects model using the log risk ratios and standard errors as input
### note: the second argument of rma() is for the *sampling variances*, so we use the
### named argument 'sei' to supply the standard errors to the function
dat$sei <- sqrt(dat$vi)
rma(yi, sei=sei, data=dat)

### fit a random-effects model supplying the 2x2 table cell frequencies to the function
rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### fit a mixed-effects model with two moderators (absolute latitude and publication year)
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

### factor() is not needed as character variables are automatically converted to factors
res <- rma(yi, vi, mods = ~ alloc + year + ablat, data=dat, btt=2:3)
res

### test all pairwise differences between the 'alloc' levels
anova(res, X=pairmat(btt="alloc"))

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
res <- rma(yi, vi, mods = ~ 0 + factor(alloc), data=dat)
res

############################################################################

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

############################################################################

### an example of a location-scale model
dat <- dat.bangertdrowns2004

### fit a standard random-effects model
res <- rma(yi, vi, data=dat)
res

### fit the same model as a location-scale model
res <- rma(yi, vi, scale = ~ 1, data=dat)
res

### check that we obtain the same estimate for tau^2
predict(res, newscale=1, transf=exp)

### add the total sample size (per 100) as a location and scale predictor
dat$ni100 <- dat$ni/100
res <- rma(yi, vi, mods = ~ ni100, scale = ~ ni100, data=dat)
res

### variables in the location and scale parts can differ
res <- rma(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100 + imag, data=dat)
res
}
\keyword{models}
