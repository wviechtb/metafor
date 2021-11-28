\name{rma.peto}
\alias{rma.peto}
\title{Meta-Analysis via Peto's Method}
\description{
   Function to fit equal-effects models to \mjeqn{2 \times 2}{2x2} table data via Peto's method. See below and the introduction to the \pkg{\link{metafor-package}} for more details on these models. \loadmathjax
}
\usage{
rma.peto(ai, bi, ci, di, n1i, n2i,
         data, slab, subset,
         add=1/2, to="only0", drop00=TRUE,
         level=95, digits, verbose=FALSE, \dots)
}
\arguments{
   \item{ai}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (upper left cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{bi}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (upper right cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ci}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (lower left cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{di}{vector to specify the \mjeqn{2 \times 2}{2x2} table frequencies (lower right cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n1i}{vector to specify the group sizes or row totals (first group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n2i}{vector to specify the group sizes or row totals (second group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \mjseqn{k} studies.}
   \item{subset}{optional (logical or numeric) vector to specify the subset of studies that should be used for the analysis.}
   \item{add}{non-negative number to specify the amount to add to zero cells, counts, or frequencies when calculating the observed effect sizes or outcomes of the individual studies. Can also be a vector of two numbers, where the first number is used in the calculation of the observed effect sizes outcomes and the second number is used when applying Peto's method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{character string to specify when the values under \code{add} should be added (either \code{"only0"}, \code{"all"}, \code{"if0all"}, or \code{"none"}). Can also be a character vector, where the first string again applies when calculating the observed effect sizes or outcomes and the second string when applying Peto's method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{drop00}{logical to specify whether studies with no cases (or only cases) in both groups should be dropped when calculating the observed effect sizes or outcomes (the outcomes for such studies are set to \code{NA}). Can also be a vector of two logicals, where the first applies to the calculation of the observed effect sizes or outcomes and the second when applying Peto's method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (the default is 95).}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is 4. See also \link[=misc-options]{here} for further details on how to control the number of digits in the output.}
   \item{verbose}{logical to specify whether output should be generated on the progress of the model fitting (the default is \code{FALSE}).}
   \item{\dots}{additional arguments.}
}
\details{

   \subsection{Specifying the Data}{

      The studies are assumed to provide data in terms of \mjeqn{2 \times 2}{2x2} tables of the form:
      \tabular{lccc}{
              \tab outcome 1 \tab outcome 2 \tab total      \cr
      group 1 \tab \code{ai} \tab \code{bi} \tab \code{n1i} \cr
      group 2 \tab \code{ci} \tab \code{di} \tab \code{n2i}
      } where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies and \code{n1i} and \code{n2i} the row totals. For example, in a set of randomized clinical trials (RCTs) or cohort studies, group 1 and group 2 may refer to the treatment (exposed) and placebo/control (not exposed) group, with outcome 1 denoting some event of interest (e.g., death) and outcome 2 its complement. In a set of case-control studies, group 1 and group 2 may refer to the group of cases and the group of controls, with outcome 1 denoting, for example, exposure to some risk factor and outcome 2 non-exposure.

   }

   \subsection{Peto's Method}{

      An approach for aggregating data of this type was suggested by Peto (see Yusuf et al., 1985). The method provides a weighted estimate of the (log) odds ratio under an equal-effects model. The method is particularly advantageous when the event of interest is rare, but it should only be used when the group sizes within the individual studies are not too dissimilar and the effect sizes are generally small (Greenland & Salvan, 1990; Sweeting et al., 2004; Bradburn et al., 2007). Note that the printed results are given both in terms of the log and the raw units (for easier interpretation).

   }

   \subsection{Observed Effect Sizes or Outcomes of the Individual Studies}{

      Peto's method itself does not require the calculation of the observed (log) odds ratios of the individual studies and directly makes use of the \mjeqn{2 \times 2}{2x2} table counts. Zero cells are not a problem (except in extreme cases, such as when one of the two outcomes never occurs in any of the tables). Therefore, it is unnecessary to add some constant to the cell counts when there are zero cells.

      However, for plotting and various other functions, it is necessary to calculate the observed (log) odds ratios for the \mjseqn{k} studies. Here, zero cells can be problematic, so adding a constant value to the cell counts ensures that all \mjseqn{k} values can be calculated. The \code{add} and \code{to} arguments are used to specify what value should be added to the cell frequencies and under what circumstances when calculating the observed (log) odds ratios and when applying Peto's method. Similarly, the \code{drop00} argument is used to specify how studies with no cases (or only cases) in both groups should be handled. The documentation of the \code{\link{escalc}} function explains how the \code{add}, \code{to}, and \code{drop00} arguments work. If only a single value for these arguments is specified (as per default), then these values are used when calculating the observed (log) odds ratios and no adjustment to the cell counts is made when applying Peto's method. Alternatively, when specifying two values for these arguments, the first value applies when calculating the observed (log) odds ratios and the second value when applying Peto's method.

      Note that \code{drop00} is set to \code{TRUE} by default. Therefore, the observed (log) odds ratios for studies where \code{ai=ci=0} or \code{bi=di=0} are set to \code{NA}. When applying Peto's method, such studies are not explicitly dropped (unless the second value of \code{drop00} argument is also set to \code{TRUE}), but this is practically not necessary, as they do not actually influence the results (assuming no adjustment to the cell/event counts are made when applying Peto's method).

   }

}
\value{
   An object of class \code{c("rma.peto","rma")}. The object is a list containing the following components:
   \item{beta}{aggregated log odds ratio.}
   \item{se}{standard error of the aggregated value.}
   \item{zval}{test statistics of the aggregated value.}
   \item{pval}{corresponding p-value.}
   \item{ci.lb}{lower bound of the confidence interval.}
   \item{ci.ub}{upper bound of the confidence interval.}
   \item{QE}{test statistic of the test for heterogeneity.}
   \item{QEp}{corresponding p-value.}
   \item{k}{number of studies included in the analysis.}
   \item{yi, vi}{the vector of individual log odds ratios and corresponding sampling variances.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values under the unrestricted and restricted likelihood.}
   \item{\dots}{some additional elements/values.}
}
\section{Methods}{
   The results of the fitted model are formatted and printed with the \code{\link[=print.rma.peto]{print}} function. If fit statistics should also be given, use \code{\link[=summary.rma]{summary}} (or use the \code{\link[=fitstats.rma]{fitstats}} function to extract them).

   The \code{\link[=residuals.rma]{residuals}}, \code{\link[=rstandard.rma.peto]{rstandard}}, and \code{\link[=rstudent.rma.peto]{rstudent}} functions extract raw and standardized residuals. Leave-one-out diagnostics can be obtained with \code{\link[=leave1out.rma.peto]{leave1out}}.

   Forest, funnel, radial, \enc{L'Abbé}{L'Abbe}, and Baujat plots can be obtained with \code{\link[=forest.rma]{forest}}, \code{\link[=funnel.rma]{funnel}}, \code{\link[=radial.rma]{radial}}, \code{\link[=labbe.rma]{labbe}}, and \code{\link[=baujat.rma]{baujat}}. The \code{\link[=qqnorm.rma.peto]{qqnorm}} function provides normal QQ plots of the standardized residuals. One can also just call \code{\link[=plot.rma.peto]{plot}} on the fitted model object to obtain various plots at once.

   A cumulative meta-analysis (i.e., adding one observation at a time) can be obtained with \code{\link[=cumul.rma.peto]{cumul}}.

   Other extractor functions include \code{\link[=coef.rma]{coef}}, \code{\link[=vcov.rma]{vcov}}, \code{\link[=logLik.rma]{logLik}}, \code{\link[=deviance.rma]{deviance}}, \code{\link[=AIC.rma]{AIC}}, and \code{\link[=BIC.rma]{BIC}}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Bradburn, M. J., Deeks, J. J., Berlin, J. A., & Localio, A. R. (2007). Much ado about nothing: A comparison of the performance of meta-analytical methods with rare events. \emph{Statistics in Medicine}, \bold{26}(1), 53--77. \verb{https://doi.org/10.1002/sim.2528}

   Greenland, S., & Salvan, A. (1990). Bias in the one-step method for pooling study results. \emph{Statistics in Medicine}, \bold{9}(3), 247--252. \verb{https://doi.org/10.1002/sim.4780090307}

   Sweeting, M. J., Sutton, A. J., & Lambert, P. C. (2004). What to add to nothing? Use and avoidance of continuity corrections in meta-analysis of sparse data. \emph{Statistics in Medicine}, \bold{23}(9), 1351--1375. \verb{https://doi.org/10.1002/sim.1761}

   Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}(5), 335--371. \verb{https://doi.org/10.1016/s0033-0620(85)80003-7}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, \code{\link{rma.mh}}, and \code{\link{rma.mv}} for other model fitting functions.

   \code{\link[metadat]{dat.collins1985a}}, \code{\link[metadat]{dat.collins1985b}}, and \code{\link[metadat]{dat.yusuf1985}} for further examples of the use of the \code{rma.peto} function.
}
\examples{
### meta-analysis of the (log) odds ratios using Peto's method
rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
}
\keyword{models}
