\name{misc-recs}
\alias{misc-recs}
\alias{misc_recs}
\title{Some Recommended Practices}
\description{
This page documents some recommended practices when working with the \pkg{metafor} package (and more generally when conducting meta-analyses). \loadmathjax
}
\details{

   \subsection{Restricted Maximum Likelihood Estimation}{

      When fitting models with the \code{\link{rma.uni}} and \code{\link{rma.mv}} functions, use of restricted maximum likelihood (REML) estimation is generally recommended. This is also the default setting (i.e., \code{method="REML"}). Various simulation studies have indicated that REML estimation tends to provide approximately unbiased estimates of the amount of heterogeneity (e.g., Langan et al., 2019; Veroniki et al., 2016; Viechtbauer, 2005), or more generally, of the variance components in more complex mixed-effects models (Harville, 1977).

      For models fitted with the \code{\link{rma.uni}} function, the empirical Bayes / Paule-Mandel estimators (i.e., \code{method="EB"} / \code{method="PM"}), which can actually be shown to be identical to each other despite their different derivations (Viechtbauer et al., 2015), also have some favorable properties. However, these estimators do not generalize in a straightforward manner to more complex models, such as those that can be fitted with the \code{\link{rma.mv}} function.

   }

   \subsection{Improved Inference Methods}{

      When fitting models with the \code{\link{rma.uni}} function, tests of individual model coefficients and the corresponding confidence intervals are by default (i.e., when \code{test="z"}) based on a standard normal distribution, while the omnibus test is based on a chi-square distribution. These inference methods may not perform nominally (i.e., the Type I error rate of tests and the coverage rate of confidence intervals may deviate from the chosen level), especially when the number of studies, \mjseqn{k}, is low. Therefore, it is highly recommended to use the method by Hartung (1999), Sidik and Jonkman (2002), and Knapp and Hartung (2003) (the Knapp-Hartung method; also referred to as the Hartung-Knapp-Sidik-Jonkman method) by setting \code{test="knha"} (or equivalently, \code{test="hksj"}). Then tests of individual coefficients and confidence intervals are based on a t-distribution with \mjseqn{k-p} degrees of freedom, while the omnibus test then uses an F-distribution with \mjseqn{m} and \mjseqn{k-p} degrees of freedom (with \mjseqn{m} denoting the number of coefficients tested and \mjseqn{p} the total number of model coefficients). Various simulation studies have shown that this method works very well in providing tests and confidence intervals with close to nominal performance (e.g., \enc{Sánchez-Meca}{Sanchez-Meca} & \enc{Marín-Martínez}{Marin-Martinez}, 2008; Viechtbauer et al., 2015).

      Alternatively, one can also conduct permutation tests using the \code{\link{permutest}} function. These also perform very well (and are, in a certain sense, \sQuote{exact} tests), but are computationally expensive.

      For models fitted with the \code{\link{rma.mv}} and \code{\link{rma.glmm}} functions, the Knapp-Hartung method and permutation tests are not available. Instead, one can set \code{test="t"} to also use t- and F-distributions for making inferences (although this does not involve the adjustment to the standard errors of the estimated model coefficients that is made as part of the Knapp-Hartung method). For \code{\link{rma.mv}}, one should also set \code{dfs="contain"}, which uses an improved method for approximating the degrees of freedom of the t- and F-distributions.

      Note that \code{test="z"} is the default for the \code{\link{rma.uni}}, \code{\link{rma.mv}}, and the \code{\link{rma.glmm}} functions. While the improved inference methods described above should ideally be the default, changing this now would break backwards compatibility.

   }

   \subsection{General Workflow for Meta-Analyses Involving Complex Dependency Structures}{

      Many meta-analyses involve observed outcomes / effect size estimates that cannot be assumed to be independent, because some estimates were computed based on the same sample of subjects (or at least a partially overlapping set). In this case, one should compute the covariances for any pair of estimates that involve (fully or partially) overlapping subjects. Doing so is difficult, but we can often construct an approximate variance-covariance matrix (say \mjseqn{V}) of such dependent estimates. This can be done with the \code{\link{vcalc}} function (and/or see the \code{\link{rcalc}} function when dealing specifically with dependent correlation coefficients). We can then fit a multivariate/multilevel model to the estimates with the \code{\link{rma.mv}} function, using \mjseqn{V} as the approximate var-cov matrix of the estimates and adding fixed and random effects to the model as deemed necessary. However, since \mjseqn{V} is often only a rough approximation (and since the random effects structure may not fully capture all dependencies in the underlying true outcomes/effects), we can then apply cluster-robust inference methods (also known as robust variance estimation) to the model. This can be done with the \code{\link{robust}} function, which also interfaces with the improved inference methods implemented in the \href{https://cran.r-project.org/package=clubSandwich}{clubSandwich} package to obtain the cluster-robust tests and confidence intervals.\mjseqn{^1} Finally, we can compute predicted outcomes (with corresponding confidence intervals) and test sets of coefficients or linear combinations thereof using the \code{\link[=predict.rma]{predict}} and \code{\link[=anova.rma]{anova}} functions. See Pustejovsky and Tipton (2022) for a paper describing such a workflow for various cases.

      To summarize, the general workflow therefore will often consist of these steps:

      \preformatted{# construct/approximate the var-cov matrix of dependent estimates
V <- vcalc(...)

# fit multivariate/multilevel model with appropriate fixed/random effects
res <- rma.mv(yi, V, mods = ~ ..., random = ~ ...)

# apply cluster-robust inference methods (robust variance estimation)
# note: use the improved methods from the clubSandwich package
sav <- robust(res, cluster = ..., clubSandwich = TRUE)
sav

# compute predicted outcomes (with corresponding CIs) as needed
predict(sav, ...)

# test sets of coefficients / linear combinations as needed
anova(sav, ...)}

      How \code{\link{vcalc}} and \code{\link{rma.mv}} should be used (and the clustering variable specified for \code{\link{robust}}) will depend on the specifics of the application.

      See \code{\link[metadat]{dat.assink2016}}, \code{\link[metadat]{dat.knapp2017}}, and \code{\link[metadat]{dat.tannersmith2016}} for some examples illustrating this workflow.

   }

   \subsection{Profile Likelihood Plots to Check Parameter Identifiability}{

      When fitting complex models, it is not guaranteed that all parameters of the model are identifiable (i.e., that there is a unique set of values for the parameters that maximizes the (restricted) likelihood function). For models fitted with the \code{\link{rma.mv}} function, this pertains especially to the variance/correlation components of the model (i.e., what is specified via the \code{random} argument). Therefore, it is strongly advised in general to do post model fitting checks to make sure that the likelihood surface around the ML/REML estimates is not flat for some combination of the parameter estimates (which would imply that the estimates are essentially arbitrary). For example, one can plot the (restricted) log-likelihood as a function of each variance/correlation component in the model to make sure that each profile plot shows a clear peak at the corresponding ML/REML estimate. The \code{\link[=profile.rma]{profile}} function can be used for this purpose. See also Raue et al. (2009) for some further discussion of parameter identifiability and the use of profile likelihoods to check for this.

      The \code{\link[=profile.rma]{profile}} function should also be used after fitting location-scale models (Viechtbauer & \enc{López-López}{Lopez-Lopez}, 2022) with the \code{\link{rma.uni}} function and after fitting selection models with the \code{\link{selmodel}} function.

   }

   ---------

   \mjseqn{^1} In small meta-analyses, the (denominator) degrees of freedom for the approximate t- and F-tests provided by the cluster-robust inference methods might be very low, in which case the tests may not be trustworthy and overly conservative (Joshi et al., 2022). Under these circumstances, one can consider the use of cluster wild bootstrapping (as implemented in the \href{https://cran.r-project.org/package=wildmeta}{wildmeta} package) as an alternative method for making inferences.

}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Hartung, J. (1999). An alternative method for meta-analysis. \emph{Biometrical Journal}, \bold{41}(8), 901--916. \verb{https://doi.org/10.1002/(SICI)1521-4036(199912)41:8<901::AID-BIMJ901>3.0.CO;2-W}

   Harville, D. A. (1977). Maximum likelihood approaches to variance component estimation and to related problems. \emph{Journal of the American Statistical Association}, \bold{72}(358), 320--338. \verb{https://doi.org/10.2307/2286796}

   Joshi, M., Pustejovsky, J. E., & Beretvas, S. N. (2022). Cluster wild bootstrapping to handle dependent effect sizes in meta-analysis with a small number of studies. \emph{Research Synthesis Methods}, \bold{13}(4), 457--477. \verb{https://doi.org/10.1002/jrsm.1554}

   Knapp, G., & Hartung, J. (2003). Improved tests for a random effects meta-regression with a single covariate. \emph{Statistics in Medicine}, \bold{22}(17), 2693--2710. \verb{https://doi.org/10.1002/sim.1482}

   Langan, D., Higgins, J. P. T., Jackson, D., Bowden, J., Veroniki, A. A., Kontopantelis, E., Viechtbauer, W. & Simmonds, M. (2019). A comparison of heterogeneity variance estimators in simulated random-effects meta-analyses. \emph{Research Synthesis Methods}, \bold{10}(1), 83--98. https://doi.org/10.1002/jrsm.1316

   Pustejovsky, J. E. & Tipton, E. (2022). Meta-analysis with robust variance estimation: Expanding the range of working models. \emph{Prevention Science}, \bold{23}, 425--438. \verb{https://doi.org/10.1007/s11121-021-01246-3}

   Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmuller, U., & Timmer, J. (2009). Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. \emph{Bioinformatics}, \bold{25}(15), 1923--1929. \verb{https://doi.org/10.1093/bioinformatics/btp358}

   \enc{Sánchez-Meca}{Sanchez-Meca}, J. & \enc{Marín-Martínez}{Marin-Martinez}, F. (2008). Confidence intervals for the overall effect size in random-effects meta-analysis. \emph{Psychological Methods}, \bold{13}(1), 31--48. \verb{https://doi.org/10.1037/1082-989x.13.1.31}

   Sidik, K. & Jonkman, J. N. (2002). A simple confidence interval for meta-analysis. \emph{Statistics in Medicine}, \bold{21}(21), 3153--3159. \verb{https://doi.org/10.1002/sim.1262}

   Veroniki, A. A., Jackson, D., Viechtbauer, W., Bender, R., Bowden, J., Knapp, G., Kuss, O., Higgins, J. P., Langan, D., & Salanti, G. (2016). Methods to estimate the between-study variance and its uncertainty in meta-analysis. \emph{Research Synthesis Methods}, \bold{7}(1), 55--79. \verb{https://doi.org/10.1002/jrsm.1164}

   Viechtbauer, W. (2005). Bias and efficiency of meta-analytic variance estimators in the random-effects model. \emph{Journal of Educational and Behavioral Statistics}, \bold{30}(3), 261--293. \verb{https://doi.org/10.3102/10769986030003261}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Viechtbauer, W., \enc{López-López}{Lopez-Lopez}, J. A., \enc{Sánchez-Meca}{Sanchez-Meca}, J., & \enc{Marín-Martínez}{Marin-Martinez}, F. (2015). A comparison of procedures to test for moderators in mixed-effects meta-regression models. \emph{Psychological Methods}, \bold{20}(3), 360--374. \verb{https://doi.org/10.1037/met0000023}

   Viechtbauer, W., & \enc{López-López}{Lopez-Lopez}, J. A. (2022). Location-scale models for meta-analysis. \emph{Research Synthesis Methods}. \bold{13}(6), 697--715. \verb{https://doi.org/10.1002/jrsm.1562}
}
\keyword{documentation}
\keyword{misc}
