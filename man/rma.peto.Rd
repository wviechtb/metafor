\name{rma.peto}
\alias{rma.peto}
\title{Meta-Analysis via Peto's Method}
\description{Function to fit fixed-effects models to \eqn{2 \times 2}{2x2} table data via Peto's method. See below and the documentation of the \pkg{\link{metafor-package}} for more details on these models.}
\usage{
rma.peto(ai, bi, ci, di, n1i, n2i,
         data, slab, subset,
         add=1/2, to="only0", drop00=TRUE,
         level=95, digits=4, verbose=FALSE, \dots)
}
\arguments{
   \item{ai}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (upper left cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{bi}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (upper right cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{ci}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (lower left cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{di}{vector to specify the \eqn{2 \times 2}{2x2} table frequencies (lower right cell). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n1i}{vector to specify the group sizes or row totals (first group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{n2i}{vector to specify the group sizes or row totals (second group). See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{data}{optional data frame containing the data supplied to the function.}
   \item{slab}{optional vector with labels for the \eqn{k} studies.}
   \item{subset}{optional vector indicating the subset of tables that should be used for the analysis. This can be a logical vector of length \eqn{k} or a numeric vector indicating the indices of the tables to include.}
   \item{add}{non-negative number indicating the amount to add to zero cells, counts, or frequencies when calculating the observed outcomes of the individual studies. Can also be a vector of two numbers, where the first number is used in the calculation of the observed outcomes and the second number is used when applying Peto's method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{to}{character string indicating when the values under \code{add} should be added (either \code{"only0"}, \code{"all"}, \code{"if0all"}, or \code{"none"}). Can also be a character vector, where the first string again applies when calculating the observed outcomes and the second string when applying Peto's method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{drop00}{logical indicating whether studies with no cases (or only cases) in both groups should be dropped when calculating the observed outcomes (the outcomes for such studies are set to \code{NA}). Can also be a vector of two logicals, where the first applies to the calculation of the observed outcomes and the second when applying Peto's method. See below and the documentation of the \code{\link{escalc}} function for more details.}
   \item{level}{numerical value between 0 and 100 specifying the confidence interval level (the default is 95).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (the default is 4).}
   \item{verbose}{logical indicating whether output should be generated on the progress of the model fitting (the default is \code{FALSE}).}
   \item{\dots}{additional arguments.}
}
\details{

   \bold{Specifying the Data}

   The studies are assumed to provide data in terms of \eqn{2 \times 2}{2x2} tables of the form:
   \tabular{lccc}{
           \tab outcome 1 \tab outcome 2 \tab total      \cr
   group 1 \tab \code{ai} \tab \code{bi} \tab \code{n1i} \cr
   group 2 \tab \code{ci} \tab \code{di} \tab \code{n2i}
   } where \code{ai}, \code{bi}, \code{ci}, and \code{di} denote the cell frequencies and \code{n1i} and \code{n2i} the row totals. For example, in a set of randomized clinical trials (RCTs) or cohort studies, group 1 and group 2 may refer to the treatment (exposed) and placebo/control (not exposed) group, with outcome 1 denoting some event of interest (e.g., death) and outcome 2 its complement. In a set of case-control studies, group 1 and group 2 may refer to the group of cases and the group of controls, with outcome 1 denoting, for example, exposure to some risk factor and outcome 2 non-exposure.

   \bold{Peto's Method}

   An approach for aggregating \eqn{2 \times 2}{2x2} table data of this type was suggested by Peto (see Yusuf et al., 1985). The method provides a weighted estimate of the (log) odds ratio under a fixed-effects model. The method is particularly advantageous when the event of interest is rare, but it should only be used when the group sizes within the individual studies are not too dissimilar and effect sizes are generally small (Greenland & Salvan, 1990; Sweeting et al., 2004; Bradburn et al., 2007). Note that the printed results are given both in terms of the log and the raw units (for easier interpretation).

   \bold{Observed Outcomes of the Individual Studies}

   Peto's method itself does not require the calculation of the observed (log) odds ratios of the individual studies and directly makes use of the \eqn{2 \times 2}{2x2} table counts. Zero cells are not a problem (except in extreme cases, such as when one of the two outcomes never occurs in any of the tables). Therefore, it is unnecessary to add some constant to the cell counts when there are zero cells.

   However, for plotting and various other functions, it is necessary to calculate the observed (log) odds ratios for the \eqn{k} studies. Here, zero cells can be problematic, so adding a constant value to the cell counts ensures that all \eqn{k} values can be calculated. The \code{add} and \code{to} arguments are used to specify what value should be added to the cell frequencies and under what circumstances when calculating the observed (log) odds ratios and when applying Peto's method. Similarly, the \code{drop00} argument is used to specify how studies with no cases (or only cases) in both groups should be handled. The documentation of the \code{\link{escalc}} function explains how the \code{add}, \code{to}, and \code{drop00} arguments work. If only a single value for these arguments is specified (as per default), then these values are used when calculating the observed (log) odds ratios and no adjustment to the cell counts is made when applying Peto's method. Alternatively, when specifying two values for these arguments, the first value applies when calculating the observed (log) odds ratios and the second value when applying Peto's method.
   
   Note that \code{drop00} is set to \code{TRUE} by default. Therefore, the observed (log) odds ratios for studies where \code{ai=ci=0} or \code{bi=di=0} are set to \code{NA}. When applying Peto's method, such studies are not explicitly dropped (unless the second value of \code{drop00} argument is also set to \code{TRUE}), but this is practically not necessary, as they do not actually influence the results (assuming no adjustment to the cell/event counts are made when applying Peto's method).
}
\value{
   An object of class \code{c("rma.peto","rma")}. The object is a list containing the following components:
   \item{beta}{aggregated log odds ratio.}
   \item{se}{standard error of the aggregated value.}
   \item{zval}{test statistics of the aggregated value.}
   \item{pval}{p-value for the test statistic.}
   \item{ci.lb}{lower bound of the confidence interval.}
   \item{ci.ub}{upper bound of the confidence interval.}
   \item{QE}{test statistic for the test of heterogeneity.}
   \item{QEp}{p-value for the test of heterogeneity.}
   \item{k}{number of tables included in the analysis.}
   \item{yi, vi}{the vector of individual log odds ratios and corresponding sampling variances.}
   \item{fit.stats}{a list with the log-likelihood, deviance, AIC, BIC, and AICc values under the unrestricted and restricted likelihood.}
   \item{\dots}{some additional elements/values.}

   The results of the fitted model are formatted and printed with the \code{\link{print.rma.peto}} function. If fit statistics should also be given, use \code{\link{summary.rma}} (or use the \code{\link{fitstats.rma}} function to extract them).

   The \code{\link{residuals.rma}}, \code{\link{rstandard.rma.peto}}, and \code{\link{rstudent.rma.peto}} functions extract raw and standardized residuals. Leave-one-out diagnostics can be obtained with \code{\link{leave1out.rma.peto}}.

   Forest, funnel, radial, L'abbé, and Baujat plots can be obtained with \code{\link{forest.rma}}, \code{\link{funnel.rma}}, \code{\link{radial.rma}}, \code{\link{labbe.rma}}, and \code{\link{baujat}}. The \code{\link{qqnorm.rma.peto}} function provides normal QQ plots of the standardized residuals. One can also just call \code{\link{plot.rma.peto}} on the fitted model object to obtain various plots at once.

   A cumulative meta-analysis (i.e., adding one observation at a time) can be obtained with \code{\link{cumul.rma.peto}}.

   Other extractor functions include \code{\link{coef.rma}}, \code{\link{vcov.rma}}, \code{\link{logLik.rma}}, \code{\link{deviance.rma}}, \code{\link{AIC.rma}}, and \code{\link{BIC.rma}}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{

   Bradburn, M. J., Deeks, J. J., Berlin, J. A., & Localio, A. R. (2007). Much ado about nothing: A comparison of the performance of meta-analytical methods with rare events. \emph{Statistics in Medicine}, \bold{26}, 53--77.

   Greenland, S., & Salvan, A. (1990). Bias in the one-step method for pooling study results. \emph{Statistics in Medicine}, \bold{9}, 247--252.

   Sweeting, M. J., Sutton, A. J., & Lambert, P. C. (2004). What to add to nothing? Use and avoidance of continuity corrections in meta-analysis of sparse data. \emph{Statistics in Medicine}, \bold{23}, 1351--1375.

   Yusuf, S., Peto, R., Lewis, J., Collins, R., & Sleight, P. (1985). Beta blockade during and after myocardial infarction: An overview of the randomized trials. \emph{Progress in Cardiovascular Disease}, \bold{27}, 335--371.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.glmm}}, \code{\link{rma.mh}}, and \code{\link{rma.mv}} for other model fitting functions.

   \code{\link{dat.collins1985a}}, \code{\link{dat.collins1985b}}, and \code{\link{dat.yusuf1985}} for further examples of the use of the \code{rma.peto} function.
}
\examples{
### meta-analysis of the (log) odds ratios using Peto's method
rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
}
\keyword{models}
