\name{misc-models}
\alias{misc-models}
\title{Fixed-Effects and Random-Effects Models in Meta-Analysis \loadmathjax}
\description{
Books and articles about meta-analysis often describe and discuss the difference between the so-called \sQuote{fixed-effects model} and the \sQuote{random-effects model} (e.g., Hedges & Olkin, 1985; Cooper et al., 2009). The former term is (mostly) avoided throughout the documentation of the \pkg{metafor} package. The term \sQuote{equal-effects model} is used instead, since it more concretely describes the main assumption underlying this model (i.e., that the underlying true effects/outcomes are homogeneous, or in other words, that they are all equal to each other). The term \sQuote{common-effects model} has also sometimes been used in the literature to describe this model and is equally descriptive.

Moreover, the term \sQuote{fixed-effects model} creates a bit of a conundrum. When authors use this term, they are really typically referring to the equal-effects model. There is however another type of model, the \sQuote{real} fixed-effects model, that is different from the equal-effects model, but now we would need to invent (unnecessarily) a different term to refer to this model. Some have done so or tried to make a distinction between the \sQuote{fixed-effect model} (without the s!) and the \sQuote{fixed-effects model}, but this subtle difference in terminology is easily overlooked/missed. Using the term \sQuote{equal-effects model} avoids this confusion and is more informative.

However, the question then remains what the real fixed-effects model is all about. The purpose of this page is to describe this model and to contrast it with the well-known random-effects model.

Assume we have a set of \mjseqn{i = 1, \ldots, k} independent studies and let \mjseqn{y_i} denote the observed value of the effect size or outcome measure in the \mjeqn{i\textrm{th}}{ith} study. Let \mjseqn{\theta_i} denote the corresponding (unknown) true effect/outcome, such that \mjdeqn{y_i | \theta_i \sim N(\theta_i, v_i).}{y_i | \theta_i ~ N(\theta_i, v_i).} In other words, the observed effect sizes or outcomes are assumed to be unbiased and normally distributed estimates of the corresponding true effects/outcomes with sampling variances equal to \mjseqn{v_i}. The \mjseqn{v_i} values are assumed to be known.

The \bold{fixed-effects model} is simply given by \mjsdeqn{y_i = \theta_i + \epsilon_i,} where the \mjseqn{\theta_i} values are the (fixed) true effects/outcomes of the \mjseqn{k} studies. Therefore, the model \sQuote{conditions} on the true effects/outcomes and provides a \emph{conditional inference} about the \mjseqn{k} studies included in the meta-analysis.

When using weighted estimation (the default in \code{\link{rma.uni}} when \code{method="FE"}), this implies that the fitted model provides an estimate of \mjdeqn{\bar{\theta}_w = \frac{\sum_{i=1}^k w_i \theta_i}{\sum_{i=1}^k w_i},}{\theta_w = \sum w_i \theta_i / \sum w_i,} that is, the \emph{weighted average} of the true effects/outcomes in the \mjseqn{k} studies, with weights equal to \mjseqn{w_i = 1/v_i}.

As an example, consider the meta-analysis by Bangert-Drowns et al. (2004) on the effectiveness of writing-to-learn interventions on academic achievement. The dataset (\code{dat.bangertdrowns2004}) includes the observed standardized mean differences (variable \code{yi}) and the corresponding sampling variances (variable \code{vi}) of the 48 studies. We can fit a fixed-effects model to these data with:

\preformatted{# copy data into 'dat'
dat <- dat.bangertdrowns2004

# fit a fixed-effects model
res <- rma(yi, vi, data=dat, method="FE")
res

# Fixed-Effects Model (k = 48)
#
# I^2 (total heterogeneity / total variability):   56.12\%
# H^2 (total variability / sampling variability):  2.28
#
# Test for Heterogeneity:
# Q(df = 47) = 107.1061, p-val < .0001
#
# Model Results:
#
# estimate      se    zval    pval   ci.lb   ci.ub
#   0.1656  0.0269  6.1499  <.0001  0.1128  0.2184}

The Q-test suggests that the underlying true standardized mean differences are heterogeneous \mjeqn{(Q(\textrm{df}=47) = 107.11, p < .0001).}{(Q(df=47) = 107.11, p < .0001).} Therefore, if we believe this to be true, then the value shown under \code{estimate} is an estimate of the inverse-variance weighted average of the true standardized mean differences of these 48 studies (i.e., \mjeqn{\hat{\bar{\theta}}_w = 0.17}{\theta-bar-hat_w = 0.17}).

One can also employ an unweighted estimation method (by setting \code{weighted=FALSE} in \code{\link{rma.uni}}), which provides an estimate of the \emph{unweighted average} of the true effects/outcomes in the \mjseqn{k} studies, that is, an estimate of \mjdeqn{\bar{\theta}_u = \frac{\sum_{i=1}^k \theta_i}{k}.}{\theta_u = \sum \theta_i / k.}

Returning to the example, we then find:

\preformatted{# fit a fixed-effects model using unweighted estimation
res <- rma(yi, vi, data=dat, method="FE", weighted=FALSE)
res

# Fixed-Effects Model (k = 48)
#
# I^2 (total heterogeneity / total variability):   56.12\%
# H^2 (total variability / sampling variability):  2.28
#
# Test for Heterogeneity:
# Q(df = 47) = 107.1061, p-val < .0001
#
# Model Results:
#
# estimate      se    zval    pval   ci.lb   ci.ub
#   0.2598  0.0380  6.8366  <.0001  0.1853  0.3343}

Therefore, the \code{estimate} value is now an estimate of the average of the true standardized mean differences of these 48 studies (i.e., \mjeqn{\hat{\bar{\theta}}_u = 0.26}{\theta-bar-hat_u = 0.26}).

For weighted estimation, one could also choose to estimate \mjeqn{\bar{\theta}_w}{\theta_w}, where the \mjseqn{w_i} values are user-defined weights (via argument \code{weights} in \code{\link{rma.uni}}). Hence, using inverse-variance weights or unit weights (as in unweighted estimation) are just special cases. It is up to the user to decide to what extent \mjeqn{\bar{\theta}_w}{\theta_w} is a meaningful parameter to estimate (regardless of the weights used).

For example, we could use the sample sizes as weights:

\preformatted{# fit a fixed-effects model using the sample sizes as weights
res <- rma(yi, vi, data=dat, method="FE", weights=ni)
res

# Fixed-Effects Model (k = 48)
#
# I^2 (total heterogeneity / total variability):   56.12\%
# H^2 (total variability / sampling variability):  2.28
#
# Test for Heterogeneity:
# Q(df = 47) = 107.1061, p-val < .0001
#
# Model Results:
#
# estimate      se    zval    pval   ci.lb   ci.ub
#   0.1719  0.0269  6.3802  <.0001  0.1191  0.2248}

We therefore obtain an estimate of the sample-size weighted average of the true standardized mean differences of these 48 studies (i.e., \mjeqn{\hat{\bar{\theta}}_w = 0.17}{\theta-bar-hat_w = 0.17}). Since the sample sizes and the inverse sampling variances are highly correlated (\code{cor(dat$ni, 1/dat$vi)} yields \code{0.999}), the results are almost identical to the ones we obtained earlier using inverse-variance weighting.

The \bold{random-effects model} does not condition on the true effects/outcomes. Instead, the \mjseqn{k} studies included in the meta-analysis are assumed to be a random sample from a larger population of studies. In rare cases, the studies included in a meta-analysis are actually sampled from a larger collection of studies. More typically, the collection of studies is a hypothetical population of an essentially infinite set of studies comprising all of the studies that have been conducted, that could have been conducted, or that may be conducted in the future. We assume that \mjeqn{\theta_i \sim N(\mu, \tau^2)}{\theta_i ~ N(\mu, \tau^2)}, that is, the true effects/outcomes in the population of studies are normally distributed with \mjseqn{\mu} denoting the average true effect/outcome and \mjseqn{\tau^2} the variance of the true effects/outcomes in the population (\mjseqn{\tau^2} is therefore often referred to as the amount of \sQuote{heterogeneity} in the true effects/outcomes). The random-effects model can also be written as \mjsdeqn{y_i = \mu + u_i + \epsilon_i,} where \mjeqn{u_i \sim N(0, \tau^2)}{u_i ~ N(0, \tau^2)} and \mjeqn{\epsilon_i \sim N(0, v_i)}{\epsilon_i ~ N(0, v_i)}. The fitted model provides estimates of \mjseqn{\mu} and \mjseqn{\tau^2}. Consequently, the random-effects model provides an \emph{unconditional inference} about the average true effect/outcome in the population of studies (from which the \mjseqn{k} studies included in the meta-analysis are assumed to be a random sample).

Fitting a random-effects model to the example data yields:

\preformatted{# fit a random-effects model (note: method="REML" is the default)
res <- rma(yi, vi, data=dat)
res

# Random-Effects Model (k = 48; tau^2 estimator: REML)
#
# tau^2 (estimated amount of total heterogeneity): 0.0499 (SE = 0.0197)
# tau (square root of estimated tau^2 value):      0.2235
# I^2 (total heterogeneity / total variability):   58.37\%
# H^2 (total variability / sampling variability):  2.40
#
# Test for Heterogeneity:
# Q(df = 47) = 107.1061, p-val < .0001
#
# Model Results:
#
# estimate      se    zval    pval   ci.lb   ci.ub
#   0.2219  0.0460  4.8209  <.0001  0.1317  0.3122}

The value shown under \code{estimate} is now an estimate of the average true standardized mean difference of studies in the population of studies from which the 48 studies included in this dataset have come (i.e., \mjeqn{\hat{\mu}_w = 0.22}{\mu-hat = 0.22}).

When using weighted estimation in the context of a random-effects model, the model is fitted with weights equal to \mjseqn{w_i = 1/(\tau^2 + v_i)}, with \mjseqn{\tau^2} replaced by its estimate (the default in \code{\link{rma.uni}} when \code{method} is set to one of the possible choices for estimating \mjseqn{\tau^2}). One can also choose unweighted estimation in the context of the random-effects model (\code{weighted=FALSE}) or specify user-defined weights (via \code{weights}), although the parameter that is estimated (i.e., \mjseqn{\mu}) remains the same regardless of the estimation method and weights used (as opposed to the fixed-effect model, where the parameter estimated is different for weighted versus unweighted estimation or when using different weights than the standard inverse-variance weights). Since weighted estimation with inverse-variance weights is most efficient, it is usually to be preferred for random-effects models (while in the fixed-effect model case, we must carefully consider whether \mjeqn{\bar{\theta}_w}{\theta_w} or \mjeqn{\bar{\theta}_u}{\theta_u} is the more meaningful parameter to estimate).

Contrary to what is often stated in the literature, it is important to realize that the fixed-effects model does \emph{not} assume that the true effects/outcomes are homogeneous (i.e., that \mjseqn{\theta_i} is equal to some common value \mjseqn{\theta} in all \mjseqn{k} studies). In other words, the fixed-effects model provides perfectly valid inferences under heterogeneity, as long as one is restricting these inferences to the set of studies included in the meta-analysis and one realizes that the model does not provide an estimate of \mjseqn{\theta} or \mjseqn{\mu}, but of \mjeqn{\bar{\theta}_w}{\theta_w} or \mjeqn{\bar{\theta}_u}{\theta_u} (depending on the estimation method).

In the special case that the true effects/outcomes are actually homogeneous (the equal-effects case), the distinction between the fixed- and random-effects models disappears, since homogeneity implies that \mjeqn{\mu = \bar{\theta}_w = \bar{\theta}_u \equiv \theta}{\mu = \theta_w = \theta_u = \theta}. However, since there is no infallible method to test whether the true effects/outcomes are really homogeneous or not, a researcher should decide on the type of inference desired before examining the data and choose the model accordingly. In fact, there is nothing wrong with fitting both the fixed- and random-effects models to the same data, since these models address different questions (i.e., what was the average effect in the studies that have been conducted and included in this meta-analysis versus what is the average effect in the larger population of studies?).

Note that fitting an equal-effects model (with \code{method="EE"}) would yield the exact same results as fitting a fixed-effects model, since the equations used to fit these two models are identical. However, the interpretation of the results is different. If we fit an equal-effects model, we make the assumption that the true effects are homogeneous and, if we believe this assumption to be justified, can interpret the estimate as an estimate of \emph{the} true effect. On the other hand, if we reject the homogeneity assumption, then we should reject the model altogether. In contrast, if we fit a fixed-effects model, we do not assume homogeneity and instead interpret the estimate as an estimate of the (weighted) average true effect of the included studies.

For further discussions of the distinction between the equal-, fixed-, and random-effects models, see Laird and Mosteller (1990) and Hedges and Vevea (1998).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Cooper, H., Hedges, L. V., & Valentine, J. C. (Eds.) (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd ed.). New York: Russell Sage Foundation.

   Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. San Diego, CA: Academic Press.

   Hedges, L. V., & Vevea, J. L. (1998). Fixed- and random-effects models in meta-analysis. \emph{Psychological Methods}, \bold{3}(4), 486--504. \verb{https://doi.org/10.1037/1082-989X.3.4.486}

   Laird, N. M., & Mosteller, F. (1990). Some statistical methods for combining experimental results. \emph{International Journal of Technology Assessment in Health Care}, \bold{6}(1), 5--30. \verb{https://doi.org/10.1017/S0266462300008916}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\keyword{documentation}
\keyword{models}
\keyword{internal}
