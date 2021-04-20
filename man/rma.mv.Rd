\name{rma.mv}
\alias{rma.mv}
\title{Meta-Analysis via Multivariate/Multilevel Linear (Mixed-Effects) Models}
\description{
   Function to fit meta-analytic multivariate/multilevel fixed- and random/mixed-effects models with or without moderators via linear (mixed-effects) models. See below and the documentation of the \pkg{\link{metafor-package}} for more details on these models. \loadmathjax
}
\usage{
rma.mv(yi, V, W, mods, random, struct="CS", intercept=TRUE,
       data, slab, subset, method="REML", test="z", dfs="residual",
       level=95, digits, btt, R, Rscale="cor",
       sigma2, tau2, rho, gamma2, phi, sparse=FALSE, verbose=FALSE, control, \dots)
}
\arguments{
   \item{yi}{vector of length \mjseqn{k} with the observed effect sizes or outcomes. See \sQuote{Details}.}
   \item{V}{vector of length \mjseqn{k} with the corresponding sampling variances or a \mjeqn{k \times k}{kxk} variance-covariance matrix of the sampling errors. See \sQuote{Details}.}
   \item{W}{optional argument to specify a vector of length \mjseqn{k} with user-defined weights or a \mjeqn{k \times k}{kxk} user-defined weight matrix. See \sQuote{Details}.}
   \item{mods}{optional argument to include one or more moderators in the model. A single moderator can be given as a vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving a matrix with \mjseqn{k} rows and as many columns as there are moderator variables. Alternatively, a model \code{\link{formula}} can be used to specify the model. See \sQuote{Details}.}
   \item{random}{either a single one-sided formula or list of one-sided formulas to specify the random-effects structure of the model. See \sQuote{Details}.}
   \item{struct}{character string to specify the variance structure of an \code{~ inner | outer} formula in the \code{random} argument. Either \code{"CS"} for compound symmetry, \code{"HCS"} for heteroscedastic compound symmetry, \code{"UN"} for an unstructured variance-covariance matrix, \code{"ID"} for a scaled identity matrix, \code{"DIAG"} for a diagonal matrix, \code{"AR"} for an AR(1) autoregressive structure, \code{"HAR"} for a heteroscedastic AR(1) autoregressive structure, \code{"CAR"} for a continuous-time autoregressive structure, or one of \code{"SPEXP"}, \code{"SPGAU"}, \code{"SPLIN"}, \code{"SPRAT"}, or \code{"SPSPH"} for one of the spatial correlation structures. See \sQuote{Details}.}
   \item{intercept}{logical to specify whether an intercept should be added to the model (the default is \code{TRUE}). Ignored when \code{mods} is a formula.}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \mjseqn{k} outcomes/studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies (or more precisely, rows of the dataset) that should be used for the analysis.}
   \item{method}{character string to specify whether the model should be fitted via maximum-likelihood (\code{"ML"}) or via restricted maximum-likelihood (\code{"REML"}) estimation. Default is \code{"REML"}.}
   \item{test}{character string to specify how test statistics and confidence intervals for the fixed effects should be computed. By default (\code{test="z"}), Wald-type tests and CIs are obtained, which are based on a standard normal distribution. When \code{test="t"}, a t-distribution is used instead. See \sQuote{Details}.}
   \item{dfs}{character string to specify how the (denominator) degrees of freedom should be calculated when \code{test="t"}. Either \code{ddf="residual"} or \code{ddf="contain"}. Can also be a numeric vector with the degrees of freedom for each model coefficient. See \sQuote{Details}.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is 4.}
   \item{btt}{optional vector of indices to specify which coefficients to include in the omnibus test of moderators. Can also be a string to grep for. See \sQuote{Details}.}
   \item{R}{an optional named list of known correlation matrices corresponding to (some of) the components specified via the \code{random} argument. See \sQuote{Details}.}
   \item{Rscale}{character string, integer, or logical to specify how matrices specified via the \code{R} argument should be scaled. See \sQuote{Details}.}
   \item{sigma2}{optional numeric vector (of the same length as the number of random intercept components specified via the \code{random} argument) to fix the corresponding \mjseqn{\sigma^2} value(s). A specific \mjseqn{\sigma^2} value can be fixed by setting the corresponding element of this argument to the desired value. A specific \mjseqn{\sigma^2} value will be estimated if the corresponding element is set equal to \code{NA}. See \sQuote{Details}.}
   \item{tau2}{optional numeric value (for \code{struct="CS"}, \code{"AR"}, \code{"CAR"}, or a spatial correlation structure) or vector (for \code{struct="HCS"}, \code{"UN"}, or \code{"HAR"}) to fix the amount of (residual) heterogeneity for the levels of the inner factor corresponding to an \code{~ inner | outer} formula specified in the \code{random} argument. A numeric value fixes a particular \mjseqn{\tau^2} value, while \code{NA} means that the value should be estimated. See \sQuote{Details}.}
   \item{rho}{optional numeric value (for \code{struct="CS"}, \code{"HCS"}, \code{"AR"}, \code{"HAR"}, \code{"CAR"}, or a spatial correlation structure) or vector (for \code{struct="UN"}) to fix the correlation between the levels of the inner factor corresponding to an \code{~ inner | outer} formula specified in the \code{random} argument. A numeric value fixes a particular \mjseqn{\rho} value, while \code{NA} means that the value should be estimated. See \sQuote{Details}.}
   \item{gamma2}{as \code{tau2} argument, but for a second \code{~ inner | outer} formula specified in the \code{random} argument. See \sQuote{Details}.}
   \item{phi}{as \code{rho} argument, but for a second \code{~ inner | outer} formula specified in the \code{random} argument. See \sQuote{Details}.}
   \item{sparse}{logical to specify whether the function should use sparse matrix objects to the extent possible (can speed up model fitting substantially for certain models). See \sQuote{Note}.}
   \item{verbose}{logical to specify whether output should be generated on the progress of the model fitting (the default is \code{FALSE}). Can also be an integer. Values > 1 generate more verbose output. See \sQuote{Note}.}
   \item{control}{optional list of control values for the estimation algorithms. If unspecified, default values are defined inside the function. See \sQuote{Note}.}
   \item{\dots}{additional arguments.}
}
\details{

   \bold{Specifying the Data}

   The function can be used in conjunction with any of the usual effect size or outcome measures used in meta-analyses (e.g., log risk ratios, log odds ratios, risk differences, mean differences, standardized mean differences, raw correlation coefficients, correlation coefficients transformed with Fisher's r-to-z transformation, and so on). Simply specify the observed effect sizes or outcomes via the \code{yi} argument and the corresponding sampling variances via the \code{V} argument. In case the sampling errors are correlated, then one can specify the entire variance-covariance matrix of the sampling errors via the \code{V} argument.

   The \code{\link{escalc}} function can be used to compute a wide variety of effect size or outcome measures (and the corresponding sampling variances) based on summary statistics. Equations for computing the covariance between the sampling errors for a variety of different effect size or outcome measures can be found, for example, in Gleser and Olkin (2009). For raw and Fisher r-to-z transformed correlations, one can find suitable equations, for example, in Steiger (1980). The latter are implemented in the \code{\link{rcalc}} function.

   \bold{Specifying Fixed Effects}

   With \code{rma.mv(yi, V)}, a fixed-effects model is fitted to the data (note: arguments \code{struct}, \code{sigma2}, \code{tau2}, \code{rho}, \code{gamma2}, \code{phi}, \code{R}, and \code{Rscale} are not relevant then and are ignored). The model is then simply given by \mjeqn{y \sim N(\theta, V)}{y ~ N(\theta, V)}, where \mjseqn{y} is a (column) vector with the observed outcomes, \mjseqn{\theta} is the (average) true outcome, and \mjseqn{V} is the variance-covariance matrix of the sampling errors (if a vector of sampling variances is provided via the \code{V} argument, then \mjseqn{V} is assumed to be diagonal).

   One or more moderators can be included in the model via the \code{mods} argument. A single moderator can be given as a (row or column) vector of length \mjseqn{k} specifying the values of the moderator. Multiple moderators are specified by giving an appropriate model matrix (i.e., \mjseqn{X}) with \mjseqn{k} rows and as many columns as there are moderator variables (e.g., \code{mods = cbind(mod1, mod2, mod3)}, where \code{mod1}, \code{mod2}, and \code{mod3} correspond to the names of the variables for the three moderator variables). The intercept is added to the model matrix by default unless \code{intercept=FALSE}.

   Alternatively, one can use the standard \code{\link{formula}} syntax to specify the model. In this case, the \code{mods} argument should be set equal to a one-sided formula of the form \code{mods = ~ model} (e.g., \code{mods = ~ mod1 + mod2 + mod3}). Interactions, polynomial terms, and factors can be easily added to the model in this manner. When specifying a model formula via the \code{mods} argument, the \code{intercept} argument is ignored. Instead, the inclusion/exclusion of the intercept is controlled by the specified formula (e.g., \code{mods = ~ mod1 + mod2 + mod3 - 1} would lead to the removal of the intercept). One can also directly specify moderators via the \code{yi} argument (e.g., \code{rma.mv(yi ~ mod1 + mod2 + mod3, V)}). In that case, the \code{mods} argument is ignored and the inclusion/exclusion of the intercept again is controlled by the specified formula.

   With moderators included, the model is then given by \mjeqn{y \sim N(X \beta, V)}{y ~ N(X \beta, V)}, where \mjseqn{X} denotes the model matrix containing the moderator values (and possibly the intercept) and \mjseqn{\beta} is a column vector containing the corresponding model coefficients. The model coefficients (i.e., \mjseqn{\beta}) are then estimated with \mjeqn{b = (X'WX')^{-1} X'Wy}{b = (X'WX)^(-1) X'Wy}, where \mjeqn{W = V^{-1}}{W = V^(-1)} is the weight matrix (without moderators, \mjseqn{X} is just a column vector of 1's). With the \code{W} argument, one can also specify user-defined weights (or a weight matrix).

   \bold{Specifying Random Effects}

   One can fit random/mixed-effects models to the data by specifying the desired random effects structure via the \code{random} argument. The \code{random} argument is either a single one-sided formula or a list of one-sided formulas. One formula type that can be specified via this argument is of the form \code{random = ~ 1 | id}. Such a formula adds a random effect corresponding to the grouping variable/factor \code{id} to the model. Outcomes with the same value/level of the \code{id} variable/factor receive the same value of the random effect, while outcomes with different values/levels of the \code{id} variable/factor are assumed to be independent. The variance component corresponding to such a formula is denoted by \mjseqn{\sigma^2}. An arbitrary number of such formulas can be specified as a list of formulas (e.g., \code{random = list(~ 1 | id1, ~ 1 | id2)}), with variance components \mjseqn{\sigma^2_1}, \mjseqn{\sigma^2_2}, and so on. Nested random effects of this form can also be added using \code{random = ~ 1 | id1/id2}, which adds a random effect corresponding to the grouping variable/factor \code{id1} and a random effect corresponding to \code{id2} within \code{id1} to the model. This can be extended to models with even more levels of nesting (e.g., \code{random = ~ 1 | id1/id2/id2}).

   Random effects of this form are useful to model clustering (and hence non-independence) induced by a multilevel structure in the data (e.g., outcomes derived from the same paper, lab, research group, or species may be more similar to each other than outcomes derived from different papers, labs, research groups, or species). See, for example, Konstantopoulos (2011) and Nakagawa and Santos (2012) for more details.

   See \code{\link{dat.konstantopoulos2011}} for an example of such a multilevel meta-analysis.

   In addition or alternatively to specifying one or multiple \code{~ 1 | id} terms, the \code{random} argument can also contain a formula of the form \code{~ inner | outer}. Outcomes with the same value/level of the \code{outer} grouping variable/factor share correlated random effects corresponding to the levels of the \code{inner} grouping variable/factor, while outcomes with different values/levels of the \code{outer} grouping variable/factor are assumed to be independent (note that the \code{inner} grouping variable must either be a factor or a character variable). The \code{struct} argument is used to specify the variance structure corresponding to the \code{inner} variable/factor. With \code{struct="CS"}, a compound symmetric structure is assumed (i.e., a single variance component \mjseqn{\tau^2} corresponding to all values/levels of the \code{inner} variable/factor and a single correlation coefficient \mjseqn{\rho} for the correlation between the different values/levels). With \code{struct="HCS"}, a heteroscedastic compound symmetric structure is assumed (with variance components \mjseqn{\tau^2_1}, \mjseqn{\tau^2_2}, and so on, corresponding to the values/levels of the \code{inner} variable/factor and a single correlation coefficient \mjseqn{\rho} for the correlation between the different values/levels). With \code{struct="UN"}, an unstructured variance-covariance matrix is assumed (with variance components \mjseqn{\tau^2_1}, \mjseqn{\tau^2_2}, and so on, corresponding to the values/levels of the \code{inner} variable/factor and correlation coefficients \mjeqn{\rho_{12}}{\rho_12}, \mjeqn{\rho_{13}}{\rho_13}, \mjeqn{\rho_{23}}{\rho_23}, and so on, for the various combinations of the values/levels of the \code{inner} variable/factor). \if{html,latex}{For example, for an \code{inner} grouping variable/factor with four levels, the three structures correspond to variance-covariance matrices of the form:}

   \if{html}{\figure{structs1.png}{options: width=700}}
   \if{latex}{\figure{structs1.png}{options: width=5.4in}}

   Structures \code{struct="ID"} and \code{struct="DIAG"} are just like \code{struct="CS"} and \code{struct="HCS"}, respectively, except that \mjseqn{\rho} is automatically set to 0, so that we either get a scaled identity matrix or a diagonal matrix.

   With the \code{outer} factor corresponding to a study identification variable and the \code{inner} factor corresponding to a variable indicating the treatment type or study arm, such a random effect could be used to estimate how strongly different treatment effects or outcomes within the same study are correlated and/or whether the amount of heterogeneity differs across different treatment types/arms. Network meta-analyses (also known as mixed treatment comparisons) will also typically require such a random effect (e.g., Salanti et al., 2008). The meta-analytic bivariate model (e.g., van Houwelingen, Arends, & Stijnen, 2002) can also be fitted in this manner (see the examples below). The \code{inner} factor could also correspond to a variable indicating different types of outcomes measured within the same study, which allows for fitting multivariate models with multiple correlated effects/outcomes per study (e.g., Berkey et al., 1998; Kalaian & Raudenbush, 1996).

   See \code{\link{dat.berkey1998}} for an example of a multivariate meta-analysis with multiple outcomes. See \code{\link{dat.hasselblad1998}} and \code{\link{dat.senn2013}} for examples of network meta-analyses.

   For meta-analyses of studies reporting outcomes at multiple time points, it may also be reasonable to assume that the true effects/outcomes are correlated over time according to an autoregressive structure (Ishak et al., 2007; Trikalinos & Olkin, 2012). For this purpose, one can also choose \code{struct="AR"}, corresponding to a structure with a single variance component \mjseqn{\tau^2} and AR(1) autocorrelation among the values of the random effect. The values of the \code{inner} variable (which does not have to be a factor here) should then reflect the various time points, with increasing values reflecting later time points. This structure assumes equally spaced time points, so the actual values of the \code{inner} variable are not relevant, only their ordering. One can also use \code{struct="HAR"}, which allows for fitting a heteroscedastic AR(1) structure (with variance components \mjseqn{\tau^2_1}, \mjseqn{\tau^2_2}, and so on). Finally, when time points are not evenly spaced, one might consider using \code{struct="CAR"} for a continuous-time autoregressive structure. \if{html,latex}{For example, for an \code{inner} grouping variable with four time points, these structures correspond to variance-covariance matrices of the form:}

   \if{html}{\figure{structs2.png}{options: width=700}}
   \if{latex}{\figure{structs2.png}{options: width=5.4in}}

   See \code{\link{dat.fine1993}} and \code{\link{dat.ishak2007}} for examples involving such structures.

   For outcomes that have a known spatial configuration, various spatial correlation structures are also available. For these structures, the formula is of the form \code{random = ~ var1 + var2 + \dots | outer}, where \code{var1}, \code{var2}, and so on are variables to indicate the spatial coordinates (e.g., longitude and latitude) based on which distances (by default Euclidean) will be computed. Let \mjseqn{d} denote the distance between two points that share the same level of the \code{outer} variable (if all true effects/outcomes are allowed to be spatially correlated, simply set \code{outer} to a variable that is a constant). Then the correlation between the true effects/outcomes corresponding to these two points is a function of \mjseqn{d} and the parameter \mjseqn{\rho}. The following table shows the types of spatial correlation structures that can be specified and the equations for the correlation. The covariance between the true effects/outcomes is then the correlation times \mjseqn{\tau^2}.

   \tabular{lllll}{
   structure          \tab \ics \tab \code{struct}  \tab \ics \tab correlation \cr
   exponential        \tab \ics \tab \code{"SPEXP"} \tab \ics \tab \mjeqn{\exp(-d/\rho)}{exp(-d/rho)} \cr
   Gaussian           \tab \ics \tab \code{"SPGAU"} \tab \ics \tab \mjeqn{\exp(-d^2/\rho^2)}{exp(-d^2/rho^2)} \cr
   linear             \tab \ics \tab \code{"SPLIN"} \tab \ics \tab \mjeqn{(1 - d/\rho) I(d < \rho)}{(1 - d/rho) I(d < rho)} \cr
   rational quadratic \tab \ics \tab \code{"SPRAT"} \tab \ics \tab \mjeqn{1 - (d/\rho)^2 / (1 + (d/\rho)^2)}{1 - (d/rho)^2 / (1 + (d/rho)^2)} \cr
   spherical          \tab \ics \tab \code{"SPSPH"} \tab \ics \tab \mjeqn{(1 - 1.5(d/\rho) + 0.5(d/\rho)^3) I(d < \rho)}{(1 - 1.5(d/rho) + 0.5(d/rho)^3) I(d < rho)}
   } Note that \mjseqn{I(d < \rho)} is equal to \mjseqn{1} if \mjseqn{d < \rho} and \mjseqn{0} otherwise. The parameterization of the various structures is based on Pinheiro and Bates (2000). Instead of Euclidean distances, one can also use other distance measures by setting (the undocumented) argument \code{dist} to either \code{"maximum"} for the maximum distance between two points (supremum norm), to \code{"manhattan"} for the absolute distance between the coordinate vectors (L1 norm), or to \code{"gcd"} for the great-circle distance (WGS84 ellipsoid method). In the latter case, only two variables, namely the longitude and latitude (in decimal degrees, with minus signs for West and South), must be specified.

   The \code{random} argument can also contain a second formula of the form \code{~ inner | outer} (but no more!). A second formula of this form works exactly described as above, but its variance components are denoted by \mjseqn{\gamma^2} and its correlation components by \mjseqn{\phi}. The \code{struct} argument should then be of length 2 to specify the variance-covariance structure for the first and second component, respectively.

   When the \code{random} argument contains a formula of the form \code{~ 1 | id}, one can use the (optional) argument \code{R} to specify a corresponding known correlation matrix for the random effect (i.e., \code{R = list(id = Cor)}, where \code{Cor} is the correlation matrix). In that case, outcomes with the same value/level of the \code{id} variable/factor receive the same value for the random effect, while outcomes with different values/levels of the \code{id} variable/factor receive values that are correlated as specified in the corresponding correlation matrix given via the \code{R} argument. The column/row names of the correlation matrix given via the \code{R} argument must therefore contain all of the values/levels of the \code{id} variable/factor. When the \code{random} argument contains multiple formulas of the form \code{~ 1 | id}, one can specify known correlation matrices for none, some, or all of those terms (e.g., with \code{random = list(~ 1 | id1, ~ 1 | id2)}, one could specify \code{R = list(id1 = Cor1)} or \code{R = list(id1 = Cor1, id2 = Cor2)}, where \code{Cor1} and \code{Cor2} are the correlation matrices corresponding to the grouping variables/factors \code{id1} and \code{id2}, respectively).

   Such a random effect with a known (or at least approximately known) correlation structure is useful in a variety of contexts. For example, such a component can be used to account for the correlations induced by the shared phylogenetic history among organisms (e.g., plants, fungi, animals). In that case, \code{~ 1 | species} is used to specify the species and argument \code{R} is used to specify the phylogenetic correlation matrix of the species studied in the meta-analysis. The corresponding variance component then indicates how much variance/heterogeneity is attributable to the specified phylogeny. See Nakagawa and Santos (2012) for more details. As another example, in a genetic meta-analysis studying disease association for several single nucleotide polymorphisms (SNPs), linkage disequilibrium (LD) among the SNPs can induce an approximately known degree of correlation among the effects/outcomes. In that case, \code{~ 1 | snp} could be used to specify the SNPs and \code{R} the corresponding LD correlation matrix for the SNPs included in the meta-analysis.

   The \code{Rscale} argument controls how matrices specified via the \code{R} argument are scaled. With \code{Rscale="none"} (or \code{Rscale=0} or \code{Rscale=FALSE}), no scaling is used. With \code{Rscale="cor"} (or \code{Rscale=1} or \code{Rscale=TRUE}), the \code{\link{cov2cor}} function is used to ensure that the matrices are correlation matrices (assuming they were covariance matrices to begin with). With \code{Rscale="cor0"} (or \code{Rscale=2}), first \code{\link{cov2cor}} is used and then the elements of each correlation matrix are scaled with \mjseqn{(R - \min(R)) / (1 - \min(R))} (this ensures that a correlation of zero in a phylogenetic correlation matrix corresponds to the split at the root node of the tree comprising the species that are actually analyzed). Finally, \code{Rscale="cov0"} (or \code{Rscale=3}) only rescales with \mjseqn{(R - \min(R))} (which ensures that a phylogenetic covariance matrix is rooted at the lowest split).

   Together with the variance-covariance matrix of the sampling errors (i.e., \mjseqn{V}), the specified random effects structure of the model implies a particular marginal variance-covariance matrix of the observed effect sizes or outcomes. Once estimates of the variance components (i.e., \mjseqn{\sigma^2}, \mjseqn{\tau^2}, \mjseqn{\rho}, \mjseqn{\gamma^2}, and/or \mjseqn{\phi}, values) have been obtained (either using maximum likelihood or restricted maximum likelihood estimation), the estimated marginal variance-covariance matrix can be constructed (denoted by \mjseqn{M}). The model coefficients (i.e., \mjseqn{\beta}) are then estimated with \mjeqn{b = (X'WX')^{-1} X'Wy}{b = (X'WX)^(-1) X'Wy}, where \mjeqn{W = M^{-1}}{W = M^(-1)} is the weight matrix. With the \code{W} argument, one can again specify user-defined weights (or a weight matrix).

   \bold{Fixing Variance Components and/or Correlations}

   Arguments \code{sigma2}, \code{tau2}, \code{rho}, \code{gamma2}, and \code{phi} can be used to fix particular variance components and/or correlations at a given value. This is useful for sensitivity analyses (e.g., for plotting the regular/restricted log-likelihood as a function of a particular variance component or correlation) or for imposing a desired variance-covariance structure on the data.

   For example, if \code{random = list(~ 1 | id1, ~ 1 | id2)}, then \code{sigma2} must be of length 2 (corresponding to \mjseqn{\sigma^2_1} and \mjseqn{\sigma^2_2}) and a fixed value can be assigned to either or both variance components. Setting a particular component to \code{NA} means that the component will be estimated by the function (e.g., \code{sigma2=c(0,NA)} would fix \mjseqn{\sigma^2_1} to 0 and estimate \mjseqn{\sigma^2_2}).

   Argument \code{tau2} is only relevant when the \code{random} argument contains an \code{~ inner | outer} formula. In that case, if the \code{tau2} argument is used, it must be either of length 1 (for \code{"CS"}, \code{"ID"}, \code{"AR"}, \code{"CAR"}, or one of the spatial correlation structures) or of the same length as the number of levels of the inner factor (for \code{"HCS"}, \code{"DIAG"}, \code{"UN"}, or \code{"HAR"}). A numeric value in the \code{tau2} argument then fixes the corresponding variance component to that value, while \code{NA} means that the component will be estimated. Similarly, if argument \code{rho} is used, it must be either of length 1 (for \code{"CS"}, \code{"HCS"}, \code{"AR"}, \code{"HAR"}, or one of the spatial correlation structures) or of length \mjseqn{lvls(lvls-1)/2} (for \code{"UN"}), where \mjseqn{lvls} denotes the number of levels of the inner factor. Again, a numeric value fixes the corresponding correlation, while \code{NA} means that the correlation will be estimated. For example, with \code{struct="CS"} and \code{rho=0}, the variance-covariance matrix of the inner factor will be diagonal with \mjseqn{\tau^2} along the diagonal. For \code{struct="UN"}, the values specified under \code{rho} should be given in column-wise order (e.g., for an \code{inner} grouping variable/factor with four levels, the order would be \mjeqn{\rho_{12}}{\rho_12}, \mjeqn{\rho_{13}}{\rho_13}, \mjeqn{\rho_{23}}{\rho_23}, \mjeqn{\rho_{14}}{\rho_14}, \mjeqn{\rho_{24}}{\rho_24}, \mjeqn{\rho_{34}}{\rho_34}.

   Similarly, arguments \code{gamma2} and \code{phi} are only relevant when the \code{random} argument contains a second \code{~ inner | outer} formula. The arguments then work exactly as described above.

   \bold{Omnibus Test of Moderators}

   For models including moderators, an omnibus test of all model coefficients is conducted that excludes the intercept (the first coefficient) if it is included in the model. If no intercept is included in the model, then the omnibus test includes all of the coefficients in the model including the first. Alternatively, one can manually specify the indices of the coefficients to test via the \code{btt} argument. For example, with \code{btt=c(3,4)}, only the third and fourth coefficient from the model would be included in the test (if an intercept is included in the model, then it corresponds to the first coefficient in the model). Instead of specifying the coefficient numbers, one can specify a string for \code{btt}. In that case, \code{\link{grep}} will be used to search for all coefficient names that match the string. The omnibus test is called the \mjseqn{Q_M}-test and follows, under the assumptions of the model, a chi-square distribution with \mjseqn{m} degrees of freedom (with \mjseqn{m} denoting the number of coefficients tested) under the null hypothesis (that the true value of all coefficients tested is equal to 0).

   \bold{Categorical Moderators}

   Categorical moderator variables can be included in the model via the \code{mods} argument in the same way that appropriately (dummy) coded categorical independent variables can be included in linear models. One can either do the dummy coding manually or use a model formula together with the \code{\link{factor}} function to let \R handle the coding automatically.

   \bold{Tests and Confidence Intervals}

   By default, tests of individual coefficients in the model (and the corresponding confidence intervals) are based on a standard normal distribution, while the omnibus test is based on a chi-square distribution (see above). As an alternative, one can set \code{test="t"}, in which case tests of individual coefficients and confidence intervals are based on a t-distribution with \mjseqn{k-p} degrees of freedom, while the omnibus test statistic then uses an F-distribution with \mjseqn{m} and \mjseqn{k-p} degrees of freedom (with \mjseqn{k} denoting the total number of estimates included in the analysis and \mjseqn{p} the total number of model coefficients including the intercept if it is present). Note that \code{test="t"} is not the same as \code{test="knha"} in \code{\link{rma.uni}}, as no adjustment to the standard errors of the estimated coefficients is made.

   The method for calculating the (denominator) degrees of freedom described above (which corresponds to \code{ddf="residual"}) is quite simplistic and may lead to tests with inflated Type I error rates and confidence intervals that are too narrow on average. As an alternative, one can set \code{ddf="contain"}, in which case the degrees of freedom for the test of a particular model coefficient, \mjseqn{b_j}, are determined by checking whether \mjseqn{x_j}, the corresponding column of the model matrix \mjseqn{X}, varies at the level corresponding to a particular random effect in the model. If such a random effect can be found, then the degrees of freedom are set to \mjseqn{l-p}, where \mjseqn{l} denotes the unique number of values/levels that the random effect can take on (i.e., for an \code{~ 1 | id} term, the unique number of values/levels of the \code{id} variable and for an \code{~ inner | outer} term, the unique number of values/levels of the \code{outer} variable). If no such random effect can be found, then \mjseqn{k-p} is used as the degrees of freedom. For the omnibus F-test, the minimum of the degrees of freedom of all coefficients involved in the test is used as the denominator degrees of freedom.

   One can also set \code{ddf} to a numeric vector with the desired values for the degrees of freedom for testing the model coefficients (e.g., if some other method for determining the degrees of freedom was used).

   \bold{Test for (Residual) Heterogeneity}

   A test for (residual) heterogeneity is automatically carried out by the function. Without moderators in the model, this test is the generalized/weighted least squares extension of Cochran's \mjseqn{Q}-test, which tests whether the variability in the observed effect sizes or outcomes is larger than one would expect based on sampling variability (and the given covariances among the sampling errors) alone. A significant test suggests that the true effect sizes or outcomes are heterogeneous. When moderators are included in the model, this is the \mjseqn{Q_E}-test for residual heterogeneity, which tests whether the variability in the observed effect sizes or outcomes that is not accounted for by the moderators included in the model is larger than one would expect based on sampling variability (and the given covariances among the sampling errors) alone.
}
\value{
   An object of class \code{c("rma.mv","rma")}. The object is a list containing the following components:
   \item{beta}{estimated coefficients of the model.}
   \item{se}{standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{variance-covariance matrix of the estimated coefficients.}
   \item{sigma2}{estimated \mjseqn{\sigma^2} value(s).}
   \item{tau2}{estimated \mjseqn{\tau^2} value(s).}
   \item{rho}{estimated \mjseqn{\rho} value(s).}
   \item{gamma2}{estimated \mjseqn{\gamma^2} value(s).}
   \item{phi}{estimated \mjseqn{\phi} value(s).}
   \item{k}{number of observed effect sizes or outcomes included in the analysis.}
   \item{p}{number of coefficients in the model (including the intercept).}
   \item{m}{number of coefficients included in the omnibus test of moderators.}
   \item{QE}{test statistic of the test for (residual) heterogeneity.}
   \item{QEp}{corresponding p-value.}
   \item{QM}{test statistic of the omnibus test of moderators.}
   \item{QMp}{corresponding p-value.}
   \item{int.only}{logical that indicates whether the model is an intercept-only model.}
   \item{yi, V, X}{the vector of outcomes, the corresponding variance-covariance matrix of the sampling errors, and the model matrix.}
   \item{M}{the estimated marginal variance-covariance matrix of the observed effect sizes or outcomes.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values.}
   \item{\dots}{some additional elements/values.}
}
\section{Methods}{
   The results of the fitted model are formatted and printed with the \code{\link{print.rma.mv}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them). Full versus reduced model comparisons in terms of fit statistics and likelihoods can be obtained with \code{\link{anova.rma}}. Wald-type tests for sets of model coefficients or linear combinations thereof can be obtained with the same function. Tests and confidence intervals based on (cluster) robust methods can be obtained with \code{\link{robust.rma.mv}}.

   Predicted/fitted values can be obtained with \code{\link{predict.rma}} and \code{\link{fitted.rma}}. For best linear unbiased predictions, see \code{\link{ranef.rma.mv}}.

   The \code{\link{residuals.rma}}, \code{\link{rstandard.rma.mv}}, and \code{\link{rstudent.rma.mv}} functions extract raw and standardized residuals. See \code{\link{influence.rma.mv}} for additional case diagnostics (e.g., to determine influential studies). For models with moderators, variance inflation factors can be obtained with \code{\link{vif.rma}}.

   Confidence intervals for any variance/correlation parameters in the model can be obtained with \code{\link{confint.rma.mv}}.

   For random/mixed-effects models, the \code{\link{profile.rma.mv}} function can be used to obtain a plot of the (restricted) log-likelihood as a function of a specific variance component or correlation parameter of the model. For models with moderators, \code{\link{regplot.rma}} draws scatter plots / bubble plots, showing the (marginal) relationship between the observed outcomes and a selected moderator from the model.

   Other extractor functions include \code{\link{coef.rma}}, \code{\link{vcov.rma}}, \code{\link{logLik.rma}}, \code{\link{deviance.rma}}, \code{\link{AIC.rma}}, \code{\link{BIC.rma}}, \code{\link{hatvalues.rma.mv}}, and \code{\link{weights.rma.mv}}.
}
\note{
   Argument \code{V} also accepts a list of variance-covariance matrices for the observed effect sizes or outcomes. From the list elements, the full (block diagonal) variance-covariance matrix is then automatically constructed.

   Model fitting is done via numerical optimization over the model parameters. By default, \code{\link{nlminb}} is used for the optimization. One can also chose a different optimizer via the \code{control} argument (e.g., \code{control=list(optimizer="optim")}). When using \code{\link{optim}}, one can set the particular method via the \code{optmethod} argument (e.g., \code{control=list(optimizer="optim", optmethod="BFGS")}). Besides \code{\link{nlminb}} and \code{\link{optim}}, one can also choose one of the optimizers from the \code{minqa} package (i.e., \code{\link[minqa]{uobyqa}}, \code{\link[minqa]{newuoa}}, or \code{\link[minqa]{bobyqa}}), one of the (derivative-free) algorithms from the \code{\link[nloptr]{nloptr}} package, the Newton-type algorithm implemented in \code{\link{nlm}}, the various algorithms implemented in the \code{dfoptim} package (\code{\link[dfoptim]{hjk}} for the Hooke-Jeeves, \code{\link[dfoptim]{nmk}} for the Nelder-Mead, and \code{\link[dfoptim]{mads}} for the Mesh Adaptive Direct Searches (MADS) algorithm), the quasi-Newton type optimizer \code{\link[ucminf]{ucminf}} from the package of the same name, or the parallelized version of the L-BFGS-B algorithm implemented in \code{\link[optimParallel]{optimParallel}} from the package of the same name.

   The optimizer name must be given as a character string (i.e., in quotes). Additional control parameters can be specified via the \code{control} argument (e.g., \code{control=list(iter.max=500, rel.tol=1e-8)}). For \code{\link[nloptr]{nloptr}}, the default is to use the BOBYQA implementation from that package with a relative convergence criterion of \code{1e-8} on the function value (i.e., log-likelihood), but this can be changed via the \code{algorithm} and \code{ftop_rel} arguments (e.g., \code{control=list(optimizer="nloptr", algorithm="NLOPT_LN_SBPLX", ftol_rel=1e-6)}). For \code{\link[optimParallel]{optimParallel}}, the control argument \code{ncpus} can be used to specify the number of cores to use for the parallelization (e.g., \code{control=list(optimizer="optimParallel", ncpus=2)}). With \code{parallel::detectCores()}, one can check on the number of available cores on the local machine.

   At the moment, the starting values are not chosen in a terribly clever way and could be far off. As a result, the optimizer may be slow to converge or may even get stuck at a local maximum. One can set the starting values manually for the various variance components and correlations in the model via the \code{control} argument by specifying the vectors \code{sigma2.init}, \code{tau2.init}, \code{rho.init}, \code{gamma2.init}, and/or \code{phi.init} as needed. Especially for complex models, it is a good idea to try out different starting values to make sure that the same estimates are obtained.

   Information on the progress of the optimization algorithm can be obtained by setting \code{verbose=TRUE} (this won't work when using parallelization). Since fitting complex models with many random effects can be computationally expensive, this option is useful to determine how the model fitting is progressing. One can also set \code{verbose} to an integer (\code{verbose=2} yields even more information and \code{verbose=3} also sets \code{option(warn=1)} temporarily).

   Whether particular variance components and/or correlations are actually identifiable needs to be carefully examined when fitting complex models. The function does some limited checking internally to fix variances and/or correlations at zero when it is clear that insufficient information is available to estimate a particular parameter (e.g., if a particular factor has only a single level, the corresponding variance component cannot be estimated). However, it is strongly advised in general to do post model fitting checks to make sure that the likelihood surface around the ML/REML estimates is not flat for some combination of the parameter estimates (which would imply that the estimates are essentially arbitrary). For example, one can plot the (restricted) log-likelihood as a function of each variance component and correlation in the model to make sure that each profile plot shows a clear peak at the corresponding ML/REML estimates. The \code{\link{profile.rma.mv}} function can be used for this purpose.

   Finally, note that the model fitting is not done in a very efficient manner at the moment, which is partly a result of allowing for crossed random effects and correlations across the entire dataset (e.g., when using the \code{R} argument). As a result, the function works directly with the entire \mjeqn{k \times k}{kxk} (marginal) variance-covariance matrix of the observed effect sizes or outcomes (instead of working with smaller blocks in a block diagonal structure). As a result, model fitting can be slow for large \mjseqn{k}. However, when the variance-covariance structure is actually sparse, a lot of speed can be gained by setting \code{sparse=TRUE}, in which case sparse matrix objects are used (via the \href{https://cran.r-project.org/package=Matrix}{Matrix} package). Also, when model fitting appears to be slow, setting \code{verbose=TRUE} is useful to obtain information on how the model fitting is progressing.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Berkey, C. S., Hoaglin, D. C., Antczak-Bouckoms, A., Mosteller, F., & Colditz, G. A. (1998). Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}, \bold{17}(22), 2537--2550. \verb{https://doi.org/10.1002/(sici)1097-0258(19981130)17:22<2537::aid-sim953>3.0.co;2-c}

   Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 357--376). New York: Russell Sage Foundation.

   van Houwelingen, H. C., Arends, L. R., & Stijnen, T. (2002). Advanced methods in meta-analysis: Multivariate approach and meta-regression. \emph{Statistics in Medicine}, \bold{21}(4), 589--624. \verb{https://doi.org/10.1002/sim.1040}

   Ishak, K. J., Platt, R. W., Joseph, L., Hanley, J. A., & Caro, J. J. (2007). Meta-analysis of longitudinal studies. \emph{Clinical Trials}, \bold{4}(5), 525--539. \verb{https://doi.org/10.1177/1740774507083567}

   Kalaian, H. A., & Raudenbush, S. W. (1996). A multivariate mixed linear model for meta-analysis. \emph{Psychological Methods}, \bold{1}(3), 227-235. \verb{https://doi.org/10.1037/1082-989X.1.3.227}

   Konstantopoulos, S. (2011). Fixed effects and variance components estimation in three-level meta-analysis. \emph{Research Synthesis Methods}, \bold{2}(1), 61--76. \verb{https://doi.org/10.1002/jrsm.35}

   Nakagawa, S., & Santos, E. S. A. (2012). Methodological issues and advances in biological meta-analysis. \emph{Evolutionary Ecology}, \bold{26}(5), 1253--1274. \verb{https://doi.org/10.1007/s10682-012-9555-5}

   Pinheiro, J. C., & Bates, D. (2000). \emph{Mixed-effects models in S and S-PLUS}. New York: Springer.

   Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. \emph{Psychological Bulletin}, \bold{87}(2), 245--251. \verb{https://doi.org/10.1037/0033-2909.87.2.245}

   Salanti, G., Higgins, J. P. T., Ades, A. E., & Ioannidis, J. P. A. (2008). Evaluation of networks of randomized trials. \emph{Statistical Methods in Medical Research}, \bold{17}(3), 279--301. \verb{https://doi.org/10.1177/0962280207080643}

   Trikalinos, T. A., & Olkin, I. (2012). Meta-analysis of effect sizes reported at multiple time points: A multivariate approach. \emph{Clinical Trials}, \bold{9}(5), 610--620. \verb{https://doi.org/10.1177/1740774512453218}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mh}}, \code{\link{rma.peto}}, and \code{\link{rma.glmm}} for other model fitting functions.

   \code{\link{dat.konstantopoulos2011}}, \code{\link{dat.hasselblad1998}}, \code{\link{dat.begg1989}}, \code{\link{dat.berkey1998}}, \code{\link{dat.fine1993}}, and \code{\link{dat.ishak2007}} for further examples of the use of the \code{rma.mv} function.
}
\examples{
### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model using rma.uni()
rma(yi, vi, data=dat)

### fit random-effects model using rma.mv()
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

### fit bivariate random-effects model using rma.mv()
res <- rma.mv(yi, vi, mods = ~ group, random = ~ group | study, struct="UN", data=dat.long)
res
}
\keyword{models}
