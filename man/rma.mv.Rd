\name{rma.mv}
\alias{rma.mv}
\title{Meta-Analysis via Multivariate/Multilevel Linear (Mixed-Effects) Models}
\description{Function to fit meta-analytic multivariate/multilevel fixed- and random/mixed-effects models with or without moderators via linear (mixed-effects) models. See below and the documentation of the \pkg{\link{metafor-package}} for more details on these models.}
\usage{
rma.mv(yi, V, W, mods, random, struct="CS", intercept=TRUE,
       data, slab, subset, method="REML", test="z",
       level=95, digits=4, btt, R, Rscale="cor",
       sigma2, tau2, rho, gamma2, phi, sparse=FALSE, verbose=FALSE, control, \dots)
}
\arguments{
   \item{yi}{vector of length \eqn{k} with the observed effect sizes or outcomes. See \sQuote{Details}.}
   \item{V}{vector of length \eqn{k} with the corresponding sampling variances or a \eqn{k \times k}{k x k} variance-covariance matrix of the sampling errors. See \sQuote{Details}.}
   \item{W}{optional argument to specify a vector of length \eqn{k} with user-defined weights or a \eqn{k \times k}{k x k} user-defined weight matrix. See \sQuote{Details}.}
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \eqn{k} rows and as many columns as there are moderator variables. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{random}{either a single one-sided formula or list of one-sided formulas to specify the random-effects structure of the model. See \sQuote{Details}.}
   \item{struct}{character string to specify the variance structure of an \code{~ inner | outer} formula in the \code{random} argument. Either \code{"CS"} for compound symmetry, \code{"HCS"} for heteroscedastic compound symmetry, \code{"UN"} for an unstructured variance-covariance matrix, \code{"ID"} for a scaled identity matrix, \code{"DIAG"} for a diagonal matrix, \code{"AR"} for an AR(1) autoregressive structure, or \code{"HAR"} for a heteroscedastic AR(1) autoregressive structure. See \sQuote{Details}.}
   \item{intercept}{logical indicating whether an intercept term should be added to the model (the default is \code{TRUE}).}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \eqn{k} outcomes/studies.}
   \item{subset}{optional vector indicating the subset of outcomes/studies that should be used for the analysis. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the observations to include.}
   \item{method}{character string specifying whether the model should be fitted via maximum-likelihood (\code{"ML"}) or via restricted maximum-likelihood (\code{"REML"}) estimation. Default is \code{"REML"}.}
   \item{test}{character string specifying how test statistics and confidence intervals for the fixed effects should be computed. By default (\code{test="z"}), Wald-type tests and CIs are obtained, which are based on a standard normal distribution. When \code{test="t"}, a t-distribution is used instead. See \sQuote{Details}.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is 95).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is 4).}
   \item{btt}{optional vector of indices specifying which coefficients to include in the omnibus test of moderators. See \sQuote{Details}.}
   \item{R}{an optional named list of known correlation matrices corresponding to (some of) the components specified via the \code{random} argument. See \sQuote{Details}.}
   \item{Rscale}{character string, integer, or logical indicating how matrices specified via the \code{R} argument should be scaled. See \sQuote{Details}.}
   \item{sigma2}{optional numerical vector (of the same length as the number of random intercept components specified via the \code{random} argument) to fix the corresponding \if{latex}{\eqn{\sigma^2}}\if{html}{\eqn{\sigma²}}\if{text}{\eqn{\sigma^2}} value(s). A specific \if{latex}{\eqn{\sigma^2}}\if{html}{\eqn{\sigma²}}\if{text}{\eqn{\sigma^2}} value can be fixed by setting the corresponding element of this argument to the desired value. A specific \if{latex}{\eqn{\sigma^2}}\if{html}{\eqn{\sigma²}}\if{text}{\eqn{\sigma^2}} value will be estimated if the corresponding element is set equal to \code{NA}. See \sQuote{Details}.}
   \item{tau2}{optional numerical value (for \code{struct="CS"} or \code{struct="AR"}) or vector (for \code{struct="HCS"}, \code{struct="UN"}, or \code{struct="HAR"}) to fix the amount of (residual) heterogeneity in the levels of the inner factor corresponding to an \code{~ inner | outer} formula specified in the \code{random} argument. A numeric value fixes a particular \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value, while \code{NA} means that the value should be estimated. See \sQuote{Details}.}
   \item{rho}{optional numerical value (for \code{struct="CS"}, \code{struct="HCS"}, \code{struct="AR"}, or \code{struct="HAR"}) or vector (for \code{struct="UN"}) to fix the correlation between levels of the inner factor corresponding to an \code{~ inner | outer} formula specified in the \code{random} argument. A numeric value fixes a particular \eqn{\rho} value, while \code{NA} means that the value should be estimated. See \sQuote{Details}.}
   \item{gamma2}{as \code{tau2} argument, but for a second \code{~ inner | outer} formula specified in the \code{random} argument. See \sQuote{Details}.}
   \item{phi}{as \code{rho} argument, but for a second \code{~ inner | outer} formula specified in the \code{random} argument. See \sQuote{Details}.}
   \item{sparse}{logical indicating whether the function should use sparse matrix objects to the extent possible (can speed up model fitting substantially for certain models). See \sQuote{Note}.}
   \item{verbose}{logical indicating whether output should be generated on the progress of the model fitting (the default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{control}{optional list of control values for the estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{additional arguments.}
}
\details{

   \bold{Specifying the Data}

   The function can be used in conjunction with any of the usual effect size or outcome measures used in meta-analyses (e.g., log odds ratios, log relative risks, risk differences, mean differences, standardized mean differences, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation, and so on). Simply specify the observed outcomes via the \code{yi} argument and the corresponding sampling variances via the \code{V} argument. In case the sampling errors are correlated, then one can specify the entire variance-covariance matrix of the sampling errors via the \code{V} argument.

   The \code{\link{escalc}} function can be used to compute a wide variety of effect size and outcome measures (and the corresponding sampling variances) based on summary statistics. Equations for computing the covariance between sampling errors for a variety of different effect size or outcome measures can be found, for example, in Gleser and Olkin (2009). For raw and Fisher's r-to-z transformed correlations, one can find suitable equations, for example, in Steiger (1980).

   \bold{Specifying Fixed Effects}

   With \code{rma.mv(yi, V)}, a fixed-effects model is fitted to the data (note: arguments \code{struct}, \code{sigma2}, \code{tau2}, \code{rho}, \code{gamma2}, \code{phi}, \code{R}, and \code{Rscale} are not relevant then and are ignored). The model is then simply given by \if{latex}{\eqn{\mathbf{y} \sim N(\mathbf{1} \beta_0, \mathbf{V})}}\if{html}{\eqn{y ~ N(1 \beta₀, V)}}\if{text}{\eqn{y ~ N(1 \beta_0, V)}}, where \eqn{\mathbf{y}}{y} is the (column) vector with the observed effect sizes or outcomes, \eqn{\mathbf{1}}{1} is a column vector of 1's, \if{latex}{\eqn{\beta_0}}\if{html}{\eqn{\beta₀}}\if{text}{\eqn{\beta_0}} is the (average) true effect size or outcome, and \eqn{\mathbf{V}}{V} is the variance-covariance matrix of the sampling errors (if a vector of sampling variances is provided via the \code{V} argument, then \eqn{\mathbf{V}}{V} is assumed to be diagonal).

   One or more moderators can be included in the model via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \eqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate model matrix (i.e., \eqn{\mathbf{X}}{X}) with \eqn{k} rows and as many columns as there are moderator variables (e.g., \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for the three moderator variables). The intercept is added to the model matrix by default unless \code{intercept=FALSE}.

   Alternatively, one can use the standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept term is controlled by the specified formula (e.g., \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept term). One can also directly specify moderators via the \code{yi} argument (e.g., \code{rma.mv(yi ~ mod1 + mod2 + mod3, V)}). In that case, the \code{mods} argument is ignored and the inclusion/exclusion of the intercept term again is controlled by the specified formula.

   With moderators included, the model is then given by \eqn{\mathbf{y} \sim N(\mathbf{X} \mathbf{\beta}, \mathbf{V})}{y ~ N(X \beta, V)}, where \eqn{\mathbf{X}}{X} denotes the model matrix containing the moderator values (and possibly the intercept) and \eqn{\mathbf{\beta}}{\beta} is a column vector containing the corresponding model coefficients. The model coefficients (i.e., \eqn{\mathbf{\beta}}{\beta}) are then estimated with \if{latex}{\eqn{\mathbf{b} = (\mathbf{X}' \mathbf{W} \mathbf{X}')^{-1} \mathbf{X}' \mathbf{W} \mathbf{y}}}\if{html}{\eqn{b = (X'WX)⁻¹ X'Wy}}\if{text}{\eqn{b = (X'WX)^(-1) X'Wy}}, where \if{latex}{\eqn{\mathbf{W} = \mathbf{V}^{-1}}}\if{html}{\eqn{W = V⁻¹}}\if{text}{\eqn{W = V^(-1)}}. With the \code{W} argument, one can also specify user-defined weights (or a weight matrix).

   \bold{Specifying Random Effects}

   One can fit random/mixed-effects models to the data by specifying the desired random effects structure via the \code{random} argument. The \code{random} argument is either a single one-sided formula or a list of one-sided formulas. One formula type that can be specified via this argument is of the form \code{random = ~ 1 | id}. Such a formula adds random effects corresponding to the grouping variable/factor \code{id} to the model. Effects or outcomes with the same value/level of the \code{id} variable/factor receive the same random effect, while effects or outcomes with different values/levels of the \code{id} variable/factor are assumed to be independent. The variance component corresponding to such a formula is denoted by \if{latex}{\eqn{\sigma^2}}\if{html}{\eqn{\sigma²}}\if{text}{\eqn{\sigma^2}}. An arbitrary number of such formulas can be specified as a list of formulas (e.g., \code{random = list(~ 1 | id1, ~ 1 | id2)}), with variance components \if{latex}{\eqn{\sigma^2_1}}\if{html}{\eqn{\sigma²₁}}\if{text}{\eqn{\sigma^2_1}}, \if{latex}{\eqn{\sigma^2_2}}\if{html}{\eqn{\sigma²₂}}\if{text}{\eqn{\sigma^2_2}}, and so on.

   Such random effects components are useful to model clustering (i.e., correlation) induced by a multilevel structure in the data (e.g., effects derived from the same paper, lab, research group, or species may be more similar to each other than effects derived from different papers, labs, research groups, or species). See, for example, Konstantopoulos (2011) and Nakagawa and Santos (2012) for more details.

   See \code{\link{dat.konstantopoulos2011}} for an example analysis with multilevel meta-analytic data.

   In addition or alternatively to specifying one or multiple \code{~ 1 | id} terms, the \code{random} argument can also contain a formula of the form \code{~ inner | outer}. Effects or outcomes with different values/levels of the \code{outer} grouping variable/factor are assumed to be independent, while effects or outcomes with the same value/level of the \code{outer} grouping variable/factor share correlated random effects corresponding to the levels of the \code{inner} grouping variable/factor (note that the \code{inner} grouping variable must be either a factor or a character variable). The \code{struct} argument is used to specify the variance structure corresponding to the \code{inner} variable/factor. With \code{struct="CS"}, a compound symmetric structure is assumed (i.e., a single variance component \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} corresponding to all values/levels of the \code{inner} variable/factor and a single correlation coefficient \eqn{\rho} for the correlation between different values/levels). With \code{struct="HCS"}, a heteroscedastic compound symmetric structure is assumed (with variance components \if{latex}{\eqn{\tau^2_1}}\if{html}{\eqn{\tau²₁}}\if{text}{\eqn{\tau^2_1}}, \if{latex}{\eqn{\tau^2_2}}\if{html}{\eqn{\tau²₂}}\if{text}{\eqn{\tau^2_2}}, and so on, corresponding to the values/levels of the \code{inner} variable/factor and a single correlation coefficient \eqn{\rho} for the correlation between different values/levels). With \code{struct="UN"}, an unstructured variance-covariance matrix is assumed (with variance components \if{latex}{\eqn{\tau^2_1}}\if{html}{\eqn{\tau²₁}}\if{text}{\eqn{\tau^2_1}}, \if{latex}{\eqn{\tau^2_2}}\if{html}{\eqn{\tau²₂}}\if{text}{\eqn{\tau^2_2}}, and so on, corresponding to the values/levels of the \code{inner} variable/factor and correlation coefficients \if{latex}{\eqn{\rho_{12}}}\if{html}{\eqn{\rho₁₂}}\if{text}{\eqn{\rho_{12}}}, \if{latex}{\eqn{\rho_{13}}}\if{html}{\eqn{\rho₁₃}}\if{text}{\eqn{\rho_{13}}}, \if{latex}{\eqn{\rho_{23}}}\if{html}{\eqn{\rho₂₃}}\if{text}{\eqn{\rho_{23}}}, and so on, for the various combinations of the values/levels of the \code{inner} variable/factor). \if{html,latex}{For example, for an \code{inner} grouping variable/factor with four levels, the three structures correspond to variance-covariance matrices of the form:}

   \if{html}{\figure{structs1.png}{options: width=700}}
   \if{latex}{\figure{structs1.png}{options: width=5.4in}}

   Structures \code{struct="ID"} and \code{struct="DIAG"} are just like \code{struct="CS"} and \code{struct="HCS"}, respectively, except that \eqn{\rho} is automatically set to 0, so that we either get a scaled identity matrix or a diagonal matrix.

   With the \code{outer} factor corresponding to a study identification variable and the \code{inner} factor corresponding to a variable indicating the treatment type or study arm, such a random effects component could be used to estimate how strongly different true effects or outcomes within the same study are correlated and/or whether the amount of heterogeneity differs across different treatment types/arms. Network meta-analyses (also called mixed treatment comparison meta-analyses) will typically require such a random effects component (e.g., Salanti et al., 2008). The meta-analytic bivariate model (e.g., van Houwelingen, Arends, & Stijnen, 2002) can also be fitted in this manner (see the examples below). The \code{inner} factor could also correspond to a variable indicating the outcome variable, which allows for fitting multivariate models when there are multiple correlated effects per study (e.g., Berkey et al., 1998; Kalaian & Raudenbush, 1996).

   See \code{\link{dat.berkey1998}} for an example of a multivariate meta-analysis with multiple outcomes. See \code{\link{dat.hasselblad1998}} and \code{\link{dat.senn2013}} for examples of network meta-analyses.

   For meta-analyses of studies reporting outcomes at multiple time points, it may also be reasonable to assume that the true effects are correlated according to an autoregressive structure over time (Ishak et al., 2007; Trikalinos & Olkin, 2012). For this purpose, one can also choose \code{struct="AR"}, corresponding to a structure with a single variance component \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} and AR(1) autocorrelation among the random effects. The values of the \code{inner} variable/factor should then reflect the various time points, with increasing values reflecting later time points. Finally, \code{struct="HAR"} allows for fitting a heteroscedastic AR(1) structure (with variance components \if{latex}{\eqn{\tau^2_1}}\if{html}{\eqn{\tau²₁}}\if{text}{\eqn{\tau^2_1}}, \if{latex}{\eqn{\tau^2_2}}\if{html}{\eqn{\tau²₂}}\if{text}{\eqn{\tau^2_2}}, and so on). \if{html,latex}{For example, for an \code{inner} grouping variable/factor with four levels (i.e., time points), these structures correspond to variance-covariance matrices of the form:}

   \if{html}{\figure{structs2.png}{options: width=700}}
   \if{latex}{\figure{structs2.png}{options: width=5.4in}}

   See \code{\link{dat.fine1993}} and \code{\link{dat.ishak2007}} for examples involving such structures.

   The \code{random} argument can also contain a second formula of the form \code{~ inner | outer} (but no more!). A second formula of this form works exactly described as above, but its variance components are denoted by \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}} and its correlation components by \eqn{\phi}. The \code{struct} argument should then be of length 2 to specify the variance-covariance structure for the first and second component.

   When the \code{random} argument contains a formula of the form \code{~ 1 | id}, one can use the (optional) argument \code{R} to specify a corresponding known correlation matrix of the random effects (i.e., \code{R = list(id = Cor)}, where \code{Cor} is the correlation matrix). In that case, effects or outcomes with the same value/level of the \code{id} variable/factor receive the same random effect, while effects or outcomes with different values/levels of the \code{id} variable/factor receive random effects that are correlated as specified in the corresponding correlation matrix given via the \code{R} argument. The column/row names of the correlation matrix given via the \code{R} argument must therefore contain all of the values/levels of the \code{id} variable/factor. When the \code{random} argument contains multiple formulas of the form \code{~ 1 | id}, one can specify known correlation matrices for none, some, or all of those terms (e.g., with \code{random = list(~ 1 | id1, ~ 1 | id2)}, one could specify \code{R = list(id1 = Cor1)} or \code{R = list(id1 = Cor1, id2 = Cor2)}, where \code{Cor1} and \code{Cor2} are the correlation matrices corresponding to the grouping variables/factors \code{id1} and \code{id2}, respectively).

   Random effects with a known (or at least approximately known) correlation structure are useful in a variety of contexts. For example, such components can be used to account for the correlations induced by a shared phylogenetic history among organisms (e.g., plants, fungi, animals). In that case, \code{~ 1 | species} is used to specify the species and argument \code{R} is used to specify the phylogenetic correlation matrix of the species studied in the meta-analysis. The corresponding variance component then indicates how much variance/heterogeneity is attributable to the specified phylogeny. See Nakagawa and Santos (2012) for more details. As another example, in a genetic meta-analysis studying disease association for several single nucleotide polymorphisms (SNPs), linkage disequilibrium (LD) among the SNPs can induce an approximately known degree of correlation among the effects. In that case, \code{~ 1 | snp} could be used to specify the SNP and \code{R} the corresponding LD correlation map for the SNPs included in the meta-analysis.

   The \code{Rscale} argument controls how matrices specified via the \code{R} argument are scaled. With \code{Rscale="none"} (or \code{Rscale=0} or \code{Rscale=FALSE}), no scaling is used. With \code{Rscale="cor"} (or \code{Rscale=1} or \code{Rscale=TRUE}), the \code{\link{cov2cor}} function is used to ensure that the matrices are correlation matrices (assuming they were covariance matrices to begin with). With \code{Rscale="cor0"} (or \code{Rscale=2}), first \code{\link{cov2cor}} is used and then the elements of each correlation matrix are scaled with \eqn{(R - min(R)) / (1 - min(R))} (this ensures that a correlation of zero in a phylogenetic correlation matrix corresponds to the split at the root node of the tree comprising the species that are actually analyzed). Finally, \code{Rscale="cov0"} (or \code{Rscale=3}) only rescales with \eqn{(R - min(R))} (which ensures that a phylogenetic covariance matrix is rooted at the lowest split).

   Together with the variance-covariance matrix of the sampling errors (i.e., \eqn{\mathbf{V}}{V}), the random effects structure of the model implies a particular marginal variance-covariance matrix of the observed outcomes. Once estimates of the variance components (i.e., \if{latex}{\eqn{\sigma^2}}\if{html}{\eqn{\sigma²}}\if{text}{\eqn{\sigma^2}}, \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}}, \eqn{\rho}, \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}}, and/or \eqn{\phi}, values) have been obtained, the estimated marginal variance-covariance matrix can be constructed (denoted by \eqn{\mathbf{M}}{M}). The model coefficients (i.e., \eqn{\mathbf{\beta}}{\beta}) are then estimated with \if{latex}{\eqn{\mathbf{b} = (\mathbf{X}' \mathbf{W} \mathbf{X}')^{-1} \mathbf{X}' \mathbf{W} \mathbf{y}}}\if{html}{\eqn{b = (X'WX)⁻¹ X'Wy}}\if{text}{\eqn{b = (X'WX)^(-1) X'Wy}}, where \if{latex}{\eqn{\mathbf{W} = \mathbf{M}^{-1}}}\if{html}{\eqn{W = M⁻¹}}\if{text}{\eqn{W = M^(-1)}}. With the \code{W} argument, one can again specify user-defined weights (or a weight matrix).

   \bold{Fixing Variance Components and/or Correlations}

   Arguments \code{sigma2}, \code{tau2}, \code{rho}, \code{gamma2}, and \code{phi} can be used to fix particular variance components and/or correlations at a given value. This is useful for sensitivity analyses (e.g., for plotting the regular/restricted log-likelihood as a function of a particular variance component or correlation) or for imposing a desired variance-covariance structure on the data.

   For example, if \code{random = list(~ 1 | id1, ~ 1 | id2)}, then \code{sigma2} must be of length 2 (corresponding to \if{latex}{\eqn{\sigma^2_1}}\if{html}{\eqn{\sigma²₁}}\if{text}{\eqn{\sigma^2_1}} and \if{latex}{\eqn{\sigma^2_2}}\if{html}{\eqn{\sigma²₂}}\if{text}{\eqn{\sigma^2_2}}) and a fixed value can be assigned to either or both variance components. Setting a particular component to \code{NA} means that the component will be estimated by the function (e.g., \code{sigma2=c(0,NA)} would fix \if{latex}{\eqn{\sigma^2_1}}\if{html}{\eqn{\sigma²₁}}\if{text}{\eqn{\sigma^2_1}} to 0 and estimate \if{latex}{\eqn{\sigma^2_2}}\if{html}{\eqn{\sigma²₂}}\if{text}{\eqn{\sigma^2_2}}).

   Argument \code{tau2} is only relevant when the \code{random} argument contains an \code{~ inner | outer} formula. In that case, if the \code{tau2} argument is used, it must be either of length 1 (for \code{struct="CS"} and \code{struct="AR"}) or of the same length as the number of levels of the inner factor (for \code{struct="HCS"}, \code{struct="UN"}, or \code{struct="HAR"}). A numeric value in the \code{tau2} argument then fixes the corresponding variance component to that value, while \code{NA} means that the component will be estimated. Similarly, if argument \code{rho} is used, it must be either of length 1 (for \code{struct="CS"}, \code{struct="HCS"}, \code{struct="AR"}, or \code{struct="HAR"}) or of length \eqn{lvls(lvls-1)/2} (for \code{struct="UN"}), where \eqn{lvls} denotes the number of levels of the inner factor. Again, a numeric value fixes the corresponding correlation, while \code{NA} means that the correlation will be estimated. For example, with \code{struct="CS"} and \code{rho=0}, the variance-covariance matrix of the inner factor will be diagonal with \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} along the diagonal. For \code{struct="UN"}, the values specified under \code{rho} should be given in column-wise order (e.g., for an \code{inner} grouping variable/factor with four levels, the order would be \if{latex}{\eqn{\rho_{12}}}\if{html}{\eqn{\rho₁₂}}\if{text}{\eqn{\rho_{12}}}, \if{latex}{\eqn{\rho_{13}}}\if{html}{\eqn{\rho₁₃}}\if{text}{\eqn{\rho_{13}}}, \if{latex}{\eqn{\rho_{23}}}\if{html}{\eqn{\rho₂₃}}\if{text}{\eqn{\rho_{23}}}, \if{latex}{\eqn{\rho_{14}}}\if{html}{\eqn{\rho₁₄}}\if{text}{\eqn{\rho_{14}}}, \if{latex}{\eqn{\rho_{24}}}\if{html}{\eqn{\rho₂₄}}\if{text}{\eqn{\rho_{24}}}, \if{latex}{\eqn{\rho_{34}}}\if{html}{\eqn{\rho₃₄}}\if{text}{\eqn{\rho_{34}}}).

   Similarly, arguments \code{gamma2} and \code{phi} are only relevant when the \code{random} argument contains a second \code{~ inner | outer} formula. The arguments then work exactly as described above.

   \bold{Omnibus Test of Parameters}

   For models including moderators, an omnibus test of all the model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, with \code{btt=c(3,4)}, only the third and fourth coefficient from the model would be included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model).

   \bold{Categorical Moderators}

   Categorical moderator variables can be included in the model via the \code{mods} argument in the same way that appropriately (dummy) coded categorical independent variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to let \R handle the coding automatically.

   \bold{Tests and Confidence Intervals}

   By default, the test statistics of the individual coefficients in the model (and the corresponding confidence intervals) are based on a standard normal distribution, while the omnibus test is based on a chi-square distribution with \eqn{m} degrees of freedom (\eqn{m} being the number of coefficients tested). As an alternative, one can set \code{test="t"}, which slightly mimics the Knapp and Hartung (2003) method by using a t-distribution with \eqn{k-p} degrees of freedom for tests of individual coefficients and confidence intervals and an F-distribution with \eqn{m} and \eqn{k-p} degrees of freedom (\eqn{p} being the total number of model coefficients including the intercept if it is present) for the omnibus test statistic.

   \bold{Test for (Residual) Heterogeneity}

   A test for (residual) heterogeneity is automatically carried out by the function. Without moderators in the model, this test is the generalized/weighted least squares extension of Cochran's \eqn{Q}-test, which tests whether the variability in the observed effect sizes or outcomes is larger than one would expect based on sampling variability (and the given covariances among the sampling errors) alone. A significant test suggests that the true effects or outcomes are heterogeneous. When moderators are included in the model, this is the \eqn{Q_E}-test for residual heterogeneity, which tests whether the variability in the observed effect sizes or outcomes that is not accounted for by the moderators included in the model is larger than one would expect based on sampling variability (and the given covariances among the sampling errors) alone.
}
\value{
   An object of class \code{c("rma.mv","rma")}. The object is a list containing the following components:
   \item{beta}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{p-values for the test statistics.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{sigma2}{estimated \if{latex}{\eqn{\sigma^2}}\if{html}{\eqn{\sigma²}}\if{text}{\eqn{\sigma^2}} value(s).}
   \item{tau2}{estimated \if{latex}{\eqn{\tau^2}}\if{html}{\eqn{\tau²}}\if{text}{\eqn{\tau^2}} value(s).}
   \item{rho}{estimated \eqn{\rho} value(s).}
   \item{gamma2}{estimated \if{latex}{\eqn{\gamma^2}}\if{html}{\eqn{\gamma²}}\if{text}{\eqn{\gamma^2}} value(s).}
   \item{phi}{estimated \eqn{\phi} value(s).}
   \item{k}{number of studies included in the model.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of coefficients.}
   \item{QE}{test statistic for the test of (residual) heterogeneity.}
   \item{QEp}{p-value for the test of (residual) heterogeneity.}
   \item{QM}{test statistic for the omnibus test of coefficients.}
   \item{QMp}{p-value for the omnibus test of coefficients.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, V, X}{the vector of outcomes, the corresponding variance-covariance matrix of the sampling errors, and the model matrix.}
   \item{M}{the estimated marginal variance-covariance matrix of the observed outcomes.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values.}
   \item{\dots}{some additional elements/values.}

   The results of the fitted model are formated and printed with the \code{\link{print.rma.mv}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them). Full versus reduced model comparisons in terms of fit statistics and likelihoods can be obtained with \code{\link{anova.rma}}. Wald-type tests for sets of model coefficients or linear combinations thereof can be obtained with the same function. Tests and confidence intervals based on (cluster) robust methods can be obtained with \code{\link{robust.rma.mv}}.

   Predicted/fitted values can be obtained with \code{\link{predict.rma}} and \code{\link{fitted.rma}}.

   The \code{\link{residuals.rma}} and \code{\link{rstandard.rma.mv}} functions extract raw and standardized residuals.

   Confidence intervals for any variance/correlation parameters in the model can be obtained with \code{\link{confint.rma.mv}}.

   For random/mixed-effects models, the \code{\link{profile.rma.mv}} function can be used to obtain a plot of the (restricted) log-likelihood as a function of a specific variance component or correlation parameter of the model.

   Other extractor functions include \code{\link{coef.rma}}, \code{\link{vcov.rma}}, \code{\link{logLik.rma}}, \code{\link{deviance.rma}}, \code{\link{AIC.rma}}, \code{\link{BIC.rma}}, \code{\link{hatvalues.rma.mv}}, and \code{\link{weights.rma.mv}}.
}
\note{
   Argument \code{V} also accepts a list of variance-covariance matrices for the observed effects or outcomes. From the list elements, the full (block diagonal) variance-covariance matrix is then automatically constructed.

   Model fitting is done via numerical optimization over the model parameters. By default, \code{\link{nlminb}} is used for the optimization. One can also chose a different optimizer via the \code{control} argument (e.g., \code{control=list(optimizer="optim")}). Note that when using \code{\link{optim}}, one can set the particular method via the \code{optmethod} argument (e.g., \code{control=list(optimizer="optim", optmethod="BFGS")}). Besides \code{\link{nlminb}} and \code{\link{optim}}, one can also choose one of the optimizers from the \code{minqa} package (i.e., \code{\link[minqa]{uobyqa}}, \code{\link[minqa]{newuoa}}, or \code{\link[minqa]{bobyqa}}), one of the (derivative-free) algorithms from the \code{\link[nloptr]{nloptr}} package, the Newton-type algorithm implemented in \code{\link{nlm}}, the Hooke-Jeeves and Nelder-Mead algorithms as implemented in \code{\link[dfoptim]{hjk}} and \code{\link[dfoptim]{nmk}} from the \code{\link[dfoptim]{dfoptim}} package, or the quasi-Newton type optimizer \code{\link[ucminf]{ucminf}} from the package of the same name.

   Note that the optimizer name must be given as a character string (i.e., in quotes). Additional control parameters for the optimizer can be specified via the \code{control} argument (e.g., \code{control=list(iter.max=500, rel.tol=1e-8)}). For \code{\link[nloptr]{nloptr}}, the default is to use the BOBYQA implementation from that package with a relative convergence criterion of \code{1e-8} on the function value (i.e., log likelihood), but this can be changed via the \code{algorithm} and \code{ftop_rel} arguments (e.g., \code{control=list(optimizer="nloptr", algorithm="NLOPT_LN_SBPLX", ftol_rel=1e-6)}).

   At the moment, the starting values are not chosen in a terribly clever way and could be far off. As a result, the optimizer may be slow to converge or may even get stuck at a local maximum. One can set the starting values manually for the various variance components and correlations in the model via the \code{control} argument by specifying the vectors \code{sigma2.init}, \code{tau2.init}, \code{rho.init}, \code{gamma2.init}, and/or \code{phi.init} as needed. Especially for complex models, it is a good idea to try out different starting values to make sure that the same estimates are obtained.

   Information on the progress of the optimization algorithm can be obtained by setting \code{verbose=TRUE} or with \code{control=list(verbose=TRUE)}.

   Whether particular variance components and/or correlations are actually identifiable needs to be carefully examined when fitting complex models. The function does some limited checking internally to fix variances and/or correlations at zero when it is clear that insufficient information is available to estimate a particular parameter (e.g., if a particular factor has only a single level, the corresponding variance component cannot be estimated). However, it is strongly advised in general to do post model fitting checks to make sure that the likelihood surface around the ML/REML estimates is not flat for some combination of the parameter estimates (which would imply that the estimates are essentially arbitrary). For example, one can plot the (restricted) log-likelihood as a function of each variance component and correlation in the model to make sure that each profile plot shows a clear peak at the corresponding ML/REML estimates. The \code{\link{profile.rma.mv}} function can be used for this purpose.

   Finally, note that the model fitting is not done in a very efficient manner at the moment, which is partly a result of allowing for crossed random effects and correlations among the effect sizes or outcomes across the entire dataset (when using the \code{R} argument). As a result, the function works directly with the entire \eqn{k \times k}{k x k} (marginal) variance-covariance matrix of the effect sizes or outcomes (instead of working with smaller blocks in a block diagonal structure). As a result, model fitting can be slow for large \eqn{k}. However, when the variance-covariance structure is actually sparse, a lot of speed can be gained by setting \code{sparse=TRUE}, in which case sparse matrix objects are used to the extent possible (which can speed up model fitting substantially for certain models). Also, when model fitting appears to be slow, setting \code{verbose=TRUE} is useful to obtain information on how the model fitting is progressing.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Berkey, C. S., Hoaglin, D. C., Antczak-Bouckoms, A., Mosteller, F., & Colditz, G. A. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}, \bold{17}, 2537--2550.

   Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 357--376). New York: Russell Sage Foundation.

   van Houwelingen, H. C., Arends, L. R., & Stijnen, T. (2002). Advanced methods in meta-analysis: Multivariate approach and meta-regression. \emph{Statistics in Medicine}, \bold{21}, 589--624.

   Ishak, K. J., Platt, R. W., Joseph, L., Hanley, J. A., & Caro, J. J. (2007). Meta-analysis of longitudinal studies. \emph{Clinical Trials}, \bold{4}, 525--539.

   Kalaian, H. A., & Raudenbush, S. W. (1996). A multivariate mixed linear model for meta-analysis. \emph{Psychological Methods}, \bold{1}, 227-235.

   Konstantopoulos, S. (2011). Fixed effects and variance components estimation in three-level meta-analysis. \emph{Research Synthesis Methods}, \bold{2}, 61--76.

   Nakagawa, S., & Santos, E. S. A. (2012). Methodological issues and advances in biological meta-analysis. \emph{Evolutionary Ecology}, \bold{26}, 1253--1274.

   Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. \emph{Psychological Bulletin}, \bold{87}, 245--251.

   Salanti, G., Higgins, J. P. T., Ades, A. E., & Ioannidis, J. P. A. (2008). Evaluation of networks of randomized trials. \emph{Statistical Methods in Medical Research}, \bold{17}, 279--301.

   Trikalinos, T. A., & Olkin, I. (2012). Meta-analysis of effect sizes reported at multiple time points: A multivariate approach. \emph{Clinical Trials}, \bold{9}, 610--620.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.glmm}} for other model fitting functions.

   \code{\link{dat.konstantopoulos2011}}, \code{\link{dat.hasselblad1998}}, \code{\link{dat.begg1989}}, \code{\link{dat.berkey1998}}, \code{\link{dat.fine1993}}, and \code{\link{dat.ishak2007}} for further examples of the use of the \code{rma.mv} function.
}
\examples{
### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### random-effects model using rma.uni()
rma(yi, vi, data=dat)

### random-effects model using rma.mv()
### note: sigma^2 in this model is the same as tau^2 from the previous model
rma.mv(yi, vi, random = ~ 1 | trial, data=dat)

### change data into long format
dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### set levels of group variable ("exp" = experimental/vaccinated; "con" = control/non-vaccinated)
levels(dat.long$group) <- c("exp", "con")

### set "con" to reference level
dat.long$group <- relevel(dat.long$group, ref="con")

### calculate log odds and corresponding sampling variances
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)

### bivariate random-effects model using rma.mv()
res <- rma.mv(yi, vi, mods = ~ group, random = ~ group | study, struct="UN", data=dat.long)
res
}
\keyword{models}
