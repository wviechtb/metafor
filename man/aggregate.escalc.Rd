\name{aggregate.escalc}
\alias{aggregate}
\alias{aggregate.escalc}
\title{Aggregate Multiple Effect Sizes or Outcomes Within Studies}
\description{
   The function can be used to aggregate multiple effect sizes or outcomes belonging to the same study (or to the same level of some other clustering variable) into a single combined effect size or outcome. \loadmathjax
}
\usage{
\method{aggregate}{escalc}(x, cluster, time, V, struct="CS", rho, phi,
          weighted=TRUE, fun, na.rm=TRUE, subset, select, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"escalc"}.}
   \item{cluster}{vector to specify the clustering variable (e.g., study).}
   \item{time}{optional vector to specify the time points (only relevant when \code{struct="CAR"} or \code{struct="CS+CAR"}).}
   \item{V}{optional argument to specify the variance-covariance matrix of the sampling errors. If not specified, argument \code{struct} is used to specify the variance-covariance structure.}
   \item{struct}{character string to specify the variance-covariance structure of the sampling errors within the same cluster (either \code{"ID"}, \code{"CS"}, \code{"CAR"}, or \code{"CS+CAR"}). See \sQuote{Details}.}
   \item{rho}{value of the correlation of the sampling errors within clusters (when \code{struct="CS"} or \code{struct="CS+CAR"}). Can also be a vector with the value of the correlation for each cluster.}
   \item{phi}{value of the autocorrelation of the sampling errors within clusters (when \code{struct="CAR"} or \code{struct="CS+CAR"}). Can also be a vector with the value of the autocorrelation for each cluster.}
   \item{weighted}{logical to specify whether estimates within clusters should be aggregated using inverse-variance weighting (the default is \code{TRUE}). If set to \code{FALSE}, unweighted averages are computed.}
   \item{fun}{optional list with three functions for aggregating other variables besides the effect sizes or outcomes within clusters (for numeric/integer variables, for logicals, and for all other types, respectively).}
   \item{na.rm}{logical to specify whether \code{NA} values should be removed before aggregating values within clusters. Can also be a vector with two logicals (the first pertains to the effect sizes or outcomes, the second to all other variables).}
   \item{subset}{optional (logical or numeric) vector to specify the subset of rows to include when aggregating the effect sizes or outcomes.}
   \item{select}{optional vector to specify the names of the variables to include in the aggregated dataset.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   In many meta-analyses, multiple effect size estimates or outcomes can be extracted from the same study. Ideally, such structures should be analyzed using an appropriate multilevel/multivariate model as can be fitted with the \code{\link{rma.mv}} function. However, there may occasionally be reasons for aggregating multiple effect sizes or outcomes belonging to the same study (or to the same level of some other clustering variable) into a single combined effect size or outcome. The present function can be used for this purpose.

   The input must be an object of class \code{"escalc"}. The error \sQuote{\code{Error in match.fun(FUN): argument "FUN" is missing, with no default}} indicates that a regular data frame was passed to the function, but this does not work. One can turn a regular data frame (containing the effect sizes or outcomes and the corresponding sampling variances) into an \code{"escalc"} object with the \code{\link{escalc}} function. See the \sQuote{Examples} below for an illustration of this.

   The \code{cluster} variable is used to specify which estimates/outcomes belong to the same study/cluster.

   In the simplest case, the estimates/outcomes within clusters (or, to be precise, their sampling errors) are assumed to be independent. This is usually a safe assumption as long as each study participant (or whatever the unit of analysis is) only contributes data to a single estimate/outcome. For example, if a study provides effect size estimates for male and female subjects separately, then the sampling errors can usually be assumed to be independent. In this case, one can set \code{struct="ID"} and multiple estimates/outcomes within the same cluster are combined using standard inverse-variance weighting (i.e., using weighted least squares) under the assumption of independence.

   In other cases, the estimates/outcomes within clusters cannot be assumed to be independent. For example, if multiple effect size estimates are computed for the same group of subjects (e.g., for different dependent variables), then the estimates are likely to be correlated. If the actual correlation between the estimates is unknown, one can often still make an educated guess and set argument \code{rho} to this value, which is then assumed to be the same for all pairs of estimates within clusters when \code{struct="CS"} (for a compound symmetric structure). Multiple estimates/outcomes within the same cluster are then combined using inverse-variance weighting taking their correlation into consideration (i.e., using generalized least squares). One can also specify a different value of \code{rho} for each cluster by passing a vector (of the same length as the number of clusters) to this argument.

   If multiple effect size estimates are computed for the same group of subjects at different time points, then it may be more sensible to assume that the correlation between estimates decreases as a function of the distance between the time points. If so, one can specify \code{struct="CAR"} (for a continuous-time autoregressive structure), set \code{phi} to the autocorrelation (for two estimates one time-unit apart), and use argument \code{time} to specify the actual time points corresponding to the estimates. The correlation between two estimates, \mjseqn{y_{ij}} and \mjseqn{y_{ij'}}, in the \mjseqn{i}th cluster, with time points \mjseqn{t_{ij}} and \mjseqn{t_{ij'}}, is then given by \mjseqn{\phi^{|t_{ij} - t_{ij'}|}}. One can also specify a different value of \code{phi} for each cluster by passing a vector (of the same length as the number of clusters) to this argument.

   One can also combine the compound symmetric and autoregressive structures by specifying \code{struct="CS+CAR"}. In this case, one must specify both \code{rho} and \code{phi}. The correlation between two estimates, \mjseqn{y_{ij}} and \mjseqn{y_{ij'}}, in the \mjseqn{i}th cluster, with time points \mjseqn{t_{ij}} and \mjseqn{t_{ij'}}, is then given by \mjseqn{\rho + (1 - \rho) \phi^{|t_{ij} - t_{ij'}|}}.

   Finally, if one actually knows the correlation (and hence the covariance) between each pair of estimates, one can also specify the entire variance-covariance matrix of the estimates (or more precisely, their sampling errors) via the \code{V} argument. In this case, arguments \code{struct}, \code{rho}, and  \code{phi} are ignored.

   Instead of using inverse-variance weighting (i.e., weighted/generalized least squares) to combine the estimates within clusters, one can set \code{weighted=FALSE} in which case the estimates are averaged within clusters without any weighting.

   Other variables (besides the estimates) will also be aggregated to the cluster level. By default, numeric/integer type variables are averaged, logicals are also averaged (yielding the proportion of \code{TRUE} values), and for all other types of variables (e.g., character variables or factors) the most frequent category/level is returned. One can also specify a list of three functions via the \code{fun} argument for aggregating variables belong to these three types.

   Argument \code{na.rm} controls how missing values should be handled. By default, any missing estimates are first removed before aggregating the non-missing values within each cluster. The same applies when aggregating the other variables. One can also specify a vector with two logicals for the \code{na.rm} argument to control how missing values should be handled when aggregating the estimates and when aggregating all other variables.
}
\value{
   An object of class \code{c("escalc","data.frame")} that contains the (selected) variables aggregated to the cluster level.

   The object is formatted and printed with the \code{\link{print.escalc}} function.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.konstantopoulos2011
dat

### aggregate estimates to the district level, assuming independent sampling
### errors for multiples studies/schools within the same district
agg <- aggregate(dat, cluster=district, struct="ID")
agg

### copy data into 'dat' and examine data
dat <- dat.assink2016
dat

### note: 'dat' is a regular data frame
class(dat)

### turn data frame into an 'escalc' object
dat <- escalc(yi=yi, vi=vi, data=dat)
class(dat)

### aggregate the estimates to the study level, assuming a CS structure for
### the sampling errors within studies with a correlation of 0.6
agg <- aggregate(dat, cluster=study, rho=0.6)
agg

### reshape 'dat.ishak2007' into long format
dat <- dat.ishak2007
dat <- reshape(dat.ishak2007, direction="long", idvar="study", v.names=c("yi","vi"),
                    varying=list(c(2,4,6,8), c(3,5,7,9)))
dat <- dat[order(dat$study, dat$time),]
is.miss  <- is.na(dat$yi)
dat <- dat[!is.miss,]
rownames(dat) <- NULL
dat

### aggregate the estimates to the study level, assuming a CAR structure for
### the sampling errors within studies with an autocorrelation of 0.9
agg <- aggregate(dat, cluster=study, struct="CAR", time=time, phi=0.9)
agg
}
\keyword{models}
