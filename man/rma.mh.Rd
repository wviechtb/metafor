\name{rma.mh}
\alias{rma.mh}
\title{Meta-Analysis via the Mantel-Haenszel Method}
\description{
   Function to fit fixed-effects models to \mjeqn{2 \times 2}{2x2} table and person-time data via the Mantel-Haenszel method. See below and the documentation of the \pkg{\link{metafor-package}} for more details on these models. \loadmathjax
}
\usage{
rma.mh(ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
       measure="OR", data, slab, subset,
       add=1/2, to="only0", drop00=TRUE,
       correct=TRUE, level=95, digits, verbose=FALSE, \dots)
}
\arguments{
   \item{ai}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (upper left cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{bi}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (upper right cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ci}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (lower left cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{di}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (lower right cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n1i}{vector to specify the group sizes or row totals (first group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n2i}{vector to specify the group sizes or row totals (second group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{x1i}{vector to specify the number of events (first group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{x2i}{vector to specify the number of events (second group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{t1i}{vector to specify the total person-times (first group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{t2i}{vector to specify the total person-times (second group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{measure}{character string to specify the outcome measure to use for the meta-analysis. Possible options are the risk ratio (\code{"RR"}), the odds ratio (\code{"OR"}), the risk difference (\code{"RD"}), the incidence rate ratio (\code{"IRR"}), or the incidence rate difference (\code{"IRD"}).}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \mjseqn{k} studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the analysis.}
   \item{add}{non-negative number to specify the amount to add to zero cells, counts, or frequencies when calculating the observed effect sizes or outcomes of the individual studies. Can also be a vector of two numbers, where the first number is used in the calculation of the observed effect sizes or outcomes and the second number is used when applying the Mantel-Haenszel method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{character string to specify when the values under \code{add} should be added (either \code{"only0"}, \code{"all"}, \code{"if0all"}, or \code{"none"}). Can also be a character vector, where the first string again applies when calculating the observed effect sizes or outcomes and the second string when applying the Mantel-Haenszel method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{drop00}{logical to specify whether studies with no cases/events (or only cases) in both groups should be dropped when calculating the observed effect sizes or outcomes (the outcomes for such studies are set to \code{NA}). Can also be a vector of two logicals, where the first applies to the calculation of the observed effect sizes or outcomes and the second when applying the Mantel-Haenszel method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{correct}{logical to specify whether to apply a continuity correction when computing the Cochran-Mantel-Haenszel test statistic.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (if unspecified, the default is 4).}
   \item{verbose}{logical to specify whether output should be generated on the progress of the model fitting (the default is \code{FALSE}).}
   \item{\dots}{additional arguments.}
}
\details{

   \bold{Specifying the Data}

   When the outcome measure is either the risk ratio (measure=\code{"RR"}), odds ratio (\code{measure="OR"}), or risk difference (\code{measure="RD"}), the studies are assumed to provide data in terms of \mjeqn{2 \times 2}{2x2} tables of the form:
   \tabular{lccc}{
           \tab outcome 1 \tab outcome 2 \tab total      \cr
   group 1 \tab \code{ai} \tab \code{bi} \tab \code{n1i} \cr
   group 2 \tab \code{ci} \tab \code{di} \tab \code{n2i}
   } where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies and \code{n1i} and \code{n2i} the row totals. For example, in a set of randomized clinical trials (RCTs) or cohort studies, group 1 and group 2 may refer to the treatment (exposed) and placebo/control (not exposed) group, with outcome 1 denoting some event of interest (e.g., death) and outcome 2 its complement. In a set of case-control studies, group 1 and group 2 may refer to the group of cases and the group of controls, with outcome 1 denoting, for example, exposure to some risk factor and outcome 2 non-exposure. For these outcome measures, one needs to specify either \code{ai}, \code{bi}, \code{ci}, and \code{di} or alternatively \code{ai}, \code{ci}, \code{n1i}, and \code{n2i}.

   Alternatively, when the outcome measure is the incidence rate ratio (\code{measure="IRR"}) or the incidence rate difference (\code{measure="IRD"}), the studies are assumed to provide data in terms of tables of the form:
   \tabular{lcc}{
           \tab events     \tab person-time \cr
   group 1 \tab \code{x1i} \tab \code{t1i}  \cr
   group 2 \tab \code{x2i} \tab \code{t2i}
   } where \code{x1i} and \code{x2i} denote the number of events in the first and the second group, respectively, and \code{t1i} and \code{t2i} the corresponding total person-times at risk.

   \bold{Mantel-Haenszel Method}

   An approach for aggregating table data of these types was suggested by Mantel and Haenszel (1959) and later extended by various authors (see references). The Mantel-Haenszel method provides a weighted estimate under a fixed-effects model. The method is particularly advantageous when aggregating a large number of studies with small sample sizes (the so-called sparse data or increasing strata case).

   When analyzing odds ratios, the Cochran-Mantel-Haenszel (CMH) test (Cochran, 1954; Mantel & Haenszel, 1959) and Tarone's test for heterogeneity (Tarone, 1985) are also provided (by default, the CMH test statistic is computed with the continuity correction; this can be switched off with \code{correct=FALSE}). When analyzing incidence rate ratios, the Mantel-Haenszel (MH) test (Rothman et al., 2008) for person-time data is also provided (again, the \code{correct} argument controls whether the continuity correction is applied). When analyzing risk ratios, odds ratios, or incidence rate ratios, the printed results are given both in terms of the log and the raw units (for easier interpretation).

   \bold{Observed Effect Sizes or Outcomes of the Individual Studies}

   The Mantel-Haenszel method itself does not require the calculation of the observed effect sizes or outcomes of the individual studies (e.g., the observed odds or incidence rate ratios of the \mjseqn{k} studies) and directly makes use of the table/event counts. Zero cells/events are not a problem (except in extreme cases, such as when one of the two outcomes never occurs in any of the \mjeqn{2 \times 2}{2x2} tables or when there are no events for one of the two groups in any of the tables). Therefore, it is unnecessary to add some constant to the cell/event counts when there are zero cells/events.

   However, for plotting and various other functions, it is necessary to calculate the observed effect sizes or outcomes for the \mjseqn{k} studies. Here, zero cells/events can be problematic, so adding a constant value to the cell/event counts ensures that all \mjseqn{k} values can be calculated. The \code{add} and \code{to} arguments are used to specify what value should be added to the cell/event counts and under what circumstances when calculating the observed effect sizes or outcomes and when applying the Mantel-Haenszel method. Similarly, the \code{drop00} argument is used to specify how studies with no cases/events (or only cases) in both groups should be handled. The documentation of the \code{\link{escalc}} function explains how the \code{add}, \code{to}, and \code{drop00} arguments work. If only a single value for these arguments is specified (as per default), then these values are used when calculating the observed effect sizes or outcomes and no adjustment to the cell/event counts is made when applying the Mantel-Haenszel method. Alternatively, when specifying two values for these arguments, the first value applies when calculating the observed effect sizes or outcomes and the second value when applying the Mantel-Haenszel method.

   Note that \code{drop00} is set to \code{TRUE} by default. Therefore, the observed effect sizes or outcomes for studies where \code{ai=ci=0} or \code{bi=di=0} or studies where \code{x1i=x2i=0} are set to \code{NA}. When applying the Mantel-Haenszel method, such studies are not explicitly dropped (unless the second value of \code{drop00} argument is also set to \code{TRUE}), but this is practically not necessary, as they do not actually influence the results (assuming no adjustment to the cell/event counts are made when applying the Mantel-Haenszel method).
}
\value{
   An object of class \code{c("rma.mh","rma")}. The object is a list containing the following components:
   \item{beta}{aggregated log risk ratio, log odds ratio, risk difference, log rate ratio, or rate difference.}
   \item{se}{standard error of the aggregated value.}
   \item{zval}{test statistics of the aggregated value.}
   \item{pval}{corresponding p-value.}
   \item{ci.lb}{lower bound of the confidence interval.}
   \item{ci.ub}{upper bound of the confidence interval.}
   \item{QE}{test statistic of the test for heterogeneity.}
   \item{QEp}{correspinding p-value.}
   \item{MH}{Cochran-Mantel-Haenszel test statistic (\code{measure="OR"}) or Mantel-Haenszel test statistic (\code{measure="IRR"}).}
   \item{MHp}{corresponding p-value.}
   \item{TA}{test statistic of Tarone's test for heterogeneity (only when \code{measure="OR"}).}
   \item{TAp}{corresponding p-value (only when \code{measure="OR"}).}
   \item{k}{number of studies included in the analysis.}
   \item{yi, vi}{the vector of outcomes and corresponding sampling variances.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values under the unrestricted and restricted likelihood.}
   \item{\dots}{some additional elements/values.}
}
\section{Methods}{
   The results of the fitted model are formatted and printed with the \code{\link{print.rma.mh}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them).

   The \code{\link{residuals.rma}}, \code{\link{rstandard.rma.mh}}, and \code{\link{rstudent.rma.mh}} functions extract raw and standardized residuals. Leave-one-out diagnostics can be obtained with \code{\link{leave1out.rma.mh}}.

   Forest, funnel, radial, \enc{L'Abbé}{L'Abbe}, and Baujat plots can be obtained with \code{\link{forest.rma}}, \code{\link{funnel.rma}}, \code{\link{radial.rma}}, \code{\link{labbe.rma}}, and \code{\link{baujat.rma}}. The \code{\link{qqnorm.rma.mh}} function provides normal QQ plots of the standardized residuals. One can also just call \code{\link{plot.rma.mh}} on the fitted model object to obtain various plots at once.

   A cumulative meta-analysis (i.e., adding one observation at a time) can be obtained with \code{\link{cumul.rma.mh}}.

   Other extractor functions include \code{\link{coef.rma}}, \code{\link{vcov.rma}}, \code{\link{logLik.rma}}, \code{\link{deviance.rma}}, \code{\link{AIC.rma}}, and \code{\link{BIC.rma}}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Cochran, W. G. (1954). Some methods for strengthening the common \mjseqn{\chi^2} tests. \emph{Biometrics}, \bold{10}(4), 417--451. \verb{https://doi.org/10.2307/3001616}

   Greenland, S., & Robins, J. M. (1985). Estimation of a common effect parameter from sparse follow-up data. \emph{Biometrics}, \bold{41}(1), 55--68. \verb{https://doi.org/10.2307/2530643}

   Mantel, N., & Haenszel, W. (1959). Statistical aspects of the analysis of data from retrospective studies of disease. \emph{Journal of the National Cancer Institute}, \bold{22}(4), 719--748. \verb{https://doi.org/10.1093/jnci/22.4.719}

   Nurminen, M. (1981). Asymptotic efficiency of general noniterative estimators of common relative risk. \emph{Biometrika}, \bold{68}(2), 525--530. \verb{https://doi.org/10.1093/biomet/68.2.525}

   Robins, J., Breslow, N., & Greenland, S. (1986). Estimators of the Mantel-Haenszel variance consistent in both sparse data and large-strata limiting models. \emph{Biometrics}, \bold{42}(2), 311--323. \verb{https://doi.org/10.2307/2531052 }

   Rothman, K. J., Greenland, S., & Lash, T. L. (2008). \emph{Modern epidemiology} (3rd ed.). Philadelphia: Lippincott Williams & Wilkins.

   Sato, T., Greenland, S., & Robins, J. M. (1989). On the variance estimator for the Mantel-Haenszel risk difference. \emph{Biometrics}, \bold{45}(4), 1323--1324. \verb{https://www.jstor.org/stable/2531784}

   Tarone, R. E. (1981). On summary estimators of relative risk. \emph{Journal of Chronic Diseases}, \bold{34}(9-10), 463--468. \verb{https://doi.org/10.1016/0021-9681(81)90006-0}

   Tarone, R. E. (1985). On heterogeneity tests based on efficient scores. \emph{Biometrika}, \bold{72}(1), 91--95. \verb{https://doi.org/10.1093/biomet/72.1.91}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, \code{\link{rma.peto}}, and \code{\link{rma.mv}} for other model fitting functions.
}
\examples{
### meta-analysis of the (log) odds ratios using the Mantel-Haenszel method
rma.mh(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### meta-analysis of the (log) risk ratios using the Mantel-Haenszel method
rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
}
\keyword{models}
