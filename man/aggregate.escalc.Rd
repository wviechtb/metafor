\name{aggregate.escalc}
\alias{aggregate}
\alias{aggregate.escalc}
\title{Aggregate Multiple Effect Sizes or Outcomes Within Studies}
\description{
   Function to aggregate multiple effect sizes or outcomes belonging to the same study (or to the same level of some other clustering variable) into a single combined effect size or outcome. \loadmathjax
}
\usage{
\method{aggregate}{escalc}(x, cluster, time, obs, V, struct="CS", rho, phi,
          weighted=TRUE, checkpd=TRUE, fun, na.rm=TRUE,
          addk=FALSE, subset, select, digits, var.names, \dots)
}
\arguments{
   \item{x}{an object of class \code{"escalc"}.}
   \item{cluster}{vector to specify the clustering variable (e.g., study).}
   \item{time}{optional vector to specify the time points (only relevant when \code{struct="CAR"}, \code{"CS+CAR"}, or \code{"CS*CAR"}).}
   \item{obs}{optional vector to distinguish different observed effect sizes or outcomes measured at the same time point (only relevant when \code{struct="CS*CAR"}).}
   \item{V}{optional argument to specify the variance-covariance matrix of the sampling errors. If unspecified, argument \code{struct} is used to specify the variance-covariance structure.}
   \item{struct}{character string to specify the variance-covariance structure of the sampling errors within the same cluster (either \code{"ID"}, \code{"CS"}, \code{"CAR"}, \code{"CS+CAR"}, or \code{"CS*CAR"}). See \sQuote{Details}.}
   \item{rho}{value of the correlation of the sampling errors within clusters (when \code{struct="CS"}, \code{"CS+CAR"}, or \code{"CS*CAR"}). Can also be a vector with the value of the correlation for each cluster.}
   \item{phi}{value of the autocorrelation of the sampling errors within clusters (when \code{struct="CAR"}, \code{"CS+CAR"}, or \code{"CS*CAR"}). Can also be a vector with the value of the autocorrelation for each cluster.}
   \item{weighted}{logical to specify whether estimates within clusters should be aggregated using inverse-variance weighting (the default is \code{TRUE}). If set to \code{FALSE}, unweighted averages are computed.}
   \item{checkpd}{logical to specify whether to check that the variance-covariance matrices of the sampling errors within clusters are positive definite (the default is \code{TRUE}).}
   \item{fun}{optional list with three functions for aggregating other variables besides the effect sizes or outcomes within clusters (for numeric/integer variables, for logicals, and for all other types, respectively).}
   \item{na.rm}{logical to specify whether \code{NA} values should be removed before aggregating values within clusters (the default is \code{TRUE}). Can also be a vector with two logicals (the first pertaining to the effect sizes or outcomes, the second to all other variables).}
   \item{addk}{logical to specify whether to add the cluster size as a new variable (called \code{ki}) to the dataset (the default is \code{FALSE}).}
   \item{subset}{optional (logical or numeric) vector to specify the subset of rows to include when aggregating the effect sizes or outcomes.}
   \item{select}{optional vector to specify the names of the variables to include in the aggregated dataset.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{var.names}{optional character vector with two elements to specify the name of the variable that contains the observed effect sizes or outcomes and the name of the variable with the corresponding sampling variances (when unspecified, the function attempts to set these automatically based on the object).}
   \item{\dots}{other arguments.}
}
\details{
   In many meta-analyses, multiple effect sizes or outcomes can be extracted from the same study. Ideally, such structures should be analyzed using an appropriate multilevel/multivariate model as can be fitted with the \code{\link{rma.mv}} function. However, there may occasionally be reasons for aggregating multiple effect sizes or outcomes belonging to the same study (or to the same level of some other clustering variable) into a single combined effect size or outcome. The present function can be used for this purpose.

   The input must be an object of class \code{"escalc"}. The error \sQuote{\code{Error in match.fun(FUN): argument "FUN" is missing, with no default}} indicates that a regular data frame was passed to the function, but this does not work. One can turn a regular data frame (containing the effect sizes or outcomes and the corresponding sampling variances) into an \code{"escalc"} object with the \code{\link{escalc}} function. See the \sQuote{Examples} below for an illustration of this.

   The \code{cluster} variable is used to specify which estimates/outcomes belong to the same study/cluster.

   In the simplest case, the estimates/outcomes within clusters (or, to be precise, their sampling errors) are assumed to be independent. This is usually a safe assumption as long as each study participant (or whatever the study units are) only contributes data to a single estimate/outcome. For example, if a study provides effect size estimates for male and female subjects separately, then the sampling errors can usually be assumed to be independent. In this case, one can set \code{struct="ID"} and multiple estimates/outcomes within the same cluster are combined using standard inverse-variance weighting (i.e., using weighted least squares) under the assumption of independence.

   In other cases, the estimates/outcomes within clusters cannot be assumed to be independent. For example, if multiple effect size estimates are computed for the same group of subjects (e.g., based on different scales to measure some construct of interest), then the estimates are likely to be correlated. If the actual correlation between the estimates is unknown, one can often still make an educated guess and set argument \code{rho} to this value, which is then assumed to be the same for all pairs of estimates within clusters when \code{struct="CS"} (for a compound symmetric structure). Multiple estimates/outcomes within the same cluster are then combined using inverse-variance weighting taking their correlation into consideration (i.e., using generalized least squares). One can also specify a different value of \code{rho} for each cluster by passing a vector (of the same length as the number of clusters) to this argument.

   If multiple effect size estimates are computed for the same group of subjects at different time points, then it may be more sensible to assume that the correlation between estimates decreases as a function of the distance between the time points. If so, one can specify \code{struct="CAR"} (for a continuous-time autoregressive structure), set \code{phi} to the autocorrelation (for two estimates one time-unit apart), and use argument \code{time} to specify the actual time points corresponding to the estimates. The correlation between two estimates, \mjeqn{y_{it}}{y_it} and \mjeqn{y_{it'}}{y_it'}, in the \mjeqn{i\textrm{th}}{ith} cluster, with time points \mjeqn{\textrm{time}_{it}}{time_it} and \mjeqn{\textrm{time}_{it'}}{time_it'}, is then given by \mjeqn{\phi^{|\textrm{time}_{it} - \textrm{time}_{it'}|}}{\phi^|time_it - time_it'|}. One can also specify a different value of \code{phi} for each cluster by passing a vector (of the same length as the number of clusters) to this argument.

   One can also combine the compound symmetric and autoregressive structures if there are multiple time points and multiple observed effect sizes or outcomes at these time points. One option is \code{struct="CS+CAR"}. In this case, one must specify the \code{time} argument and both \code{rho} and \code{phi}. The correlation between two estimates, \mjeqn{y_{it}}{y_it} and \mjeqn{y_{it'}}{y_it'}, in the \mjeqn{i\textrm{th}}{ith} cluster, with time points \mjeqn{\textrm{time}_{it}}{time_it} and \mjeqn{\textrm{time}_{it'}}{time_it'}, is then given by \mjeqn{\rho + (1 - \rho) \phi^{|\textrm{time}_{it} - \textrm{time}_{it'}|}}{\rho + (1 - \rho) * \phi^|time_it - time_it'|}.

   Alternatively, one can specify \code{struct="CS*CAR"}. In this case, one must specify both the \code{time} and \code{obs} arguments and both \code{rho} and \code{phi}. The correlation between two estimates, \mjeqn{y_{ijt}}{y_ijt} and \mjeqn{y_{ijt'}}{y_ijt'}, with the same value for \code{obs} but different values for \code{time}, is then given by \mjeqn{\phi^{|\textrm{time}_{ijt} - \textrm{time}_{ijt'}|}}{\phi^|time_ijt - time_ijt'|}, the correlation between two estimates, \mjeqn{y_{ijt}}{y_ijt} and \mjeqn{y_{ij't}}{y_ij't}, with different values for \code{obs} but the same value for \code{time}, is then given by \mjseqn{\rho}, and the correlation between two estimates, \mjeqn{y_{ijt}}{y_ijt} and \mjeqn{y_{ij't'}}{y_ij't}, with different values for \code{obs} and different values for \code{time}, is then given by \mjeqn{\rho \times \phi^{|\textrm{time}_{ijt} - \textrm{time}_{ijt'}|}}{\rho * \phi^|time_ijt - time_ijt'|}.

   Finally, if one actually knows the correlation (and hence the covariance) between each pair of estimates (or has an approximation thereof), one can also specify the entire variance-covariance matrix of the estimates (or more precisely, their sampling errors) via the \code{V} argument (in this case, arguments \code{struct}, \code{time}, \code{obs}, \code{rho}, and \code{phi} are ignored). Note that the \code{\link{vcalc}} function can be used to construct such a \code{V} matrix and provides even more flexibility for specifying various types of dependencies. See the \sQuote{Examples} below for an illustration of this.

   Instead of using inverse-variance weighting (i.e., weighted/generalized least squares) to combine the estimates within clusters, one can set \code{weighted=FALSE} in which case the estimates are averaged within clusters without any weighting (although the correlations between estimates as specified are still taken into consideration).

   Other variables (besides the estimates) will also be aggregated to the cluster level. By default, numeric/integer type variables are averaged, logicals are also averaged (yielding the proportion of \code{TRUE} values), and for all other types of variables (e.g., character variables or factors) the most frequent category/level is returned. One can also specify a list of three functions via the \code{fun} argument for aggregating variables belonging to these three types.

   Argument \code{na.rm} controls how missing values should be handled. By default, any missing estimates are first removed before aggregating the non-missing values within each cluster. The same applies when aggregating the other variables. One can also specify a vector with two logicals for the \code{na.rm} argument to control how missing values should be handled when aggregating the estimates and when aggregating all other variables.
}
\value{
   An object of class \code{c("escalc","data.frame")} that contains the (selected) variables aggregated to the cluster level.

   The object is formatted and printed with the \code{\link[=print.escalc]{print}} function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to create \code{escalc} objects.
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.konstantopoulos2011
head(dat, 11)

### aggregate estimates to the district level, assuming independent sampling
### errors for multiples studies/schools within the same district
agg <- aggregate(dat, cluster=district, struct="ID", addk=TRUE)
agg

### copy data into 'dat' and examine data
dat <- dat.assink2016
head(dat, 19)

### note: 'dat' is an 'escalc' object
class(dat)

### turn 'dat' into a regular data frame
dat <- as.data.frame(dat)
class(dat)

### turn data frame into an 'escalc' object
dat <- escalc(measure="SMD", yi=yi, vi=vi, data=dat)
class(dat)

### aggregate the estimates to the study level, assuming a CS structure for
### the sampling errors within studies with a correlation of 0.6
agg <- aggregate(dat, cluster=study, rho=0.6)
agg

### use vcalc() and then the V argument
V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)
agg <- aggregate(dat, cluster=study, V=V)
agg

### use a correlation of 0.7 for effect sizes corresponding to the same type of
### delinquent behavior and a correlation of 0.5 for effect sizes corresponding
### to different types of delinquent behavior
V <- vcalc(vi, cluster=study, type=deltype, obs=esid, data=dat, rho=c(0.7, 0.5))
agg <- aggregate(dat, cluster=study, V=V)
agg

### reshape 'dat.ishak2007' into long format
dat <- dat.ishak2007
dat <- reshape(dat.ishak2007, direction="long", idvar="study", v.names=c("yi","vi"),
               varying=list(c(2,4,6,8), c(3,5,7,9)))
dat <- dat[order(study, time),]
dat <- dat[!is.na(yi),]
rownames(dat) <- NULL
head(dat, 8)

### aggregate the estimates to the study level, assuming a CAR structure for
### the sampling errors within studies with an autocorrelation of 0.9
agg <- aggregate(dat, cluster=study, struct="CAR", time=time, phi=0.9)
head(agg, 5)
}
\keyword{models}
