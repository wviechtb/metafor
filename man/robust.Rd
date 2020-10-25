\name{robust}
\alias{robust}
\alias{robust.rma.uni}
\alias{robust.rma.mv}
\title{(Cluster) Robust Tests and Confidence Intervals for 'rma' Objects}
\description{
   The function provides (cluster) robust tests and confidence intervals of the model coefficients for objects of class \code{"rma"}. \loadmathjax
}
\usage{
robust(x, cluster, \dots)

\method{robust}{rma.uni}(x, cluster, adjust=TRUE, digits, \dots)
\method{robust}{rma.mv}(x, cluster, adjust=TRUE, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"} or \code{"rma.mv"}.}
   \item{cluster}{a vector to specify a clustering variable to use for constructing the sandwich estimator of the variance-covariance matrix.}
   \item{adjust}{logical to specify whether a small-sample correction should be applied to the variance-covariance matrix.}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The function constructs a (cluster) robust estimate of the variance-covariance matrix of the model coefficients based on a sandwich-type estimator and then computes tests and confidence intervals of the model coefficients.

   Tests of individual coefficients and confidence intervals are based on a t-distribution with \mjseqn{n-p} degrees of freedom is used, while the omnibus test statistic uses an F-distribution with \mjseqn{m} and \mjseqn{n-p} degrees of freedom, where \mjseqn{n} is the number of clusters, \mjseqn{p} denotes the total number of model coefficients (including the intercept if it is present), and \mjseqn{m} denotes the number of coefficients tested (in the omnibus test).

   When \code{adjust=TRUE} (the default), the (cluster) robust estimate of the variance-covariance matrix is multiplied by the factor \mjseqn{n/(n-p)}, which serves as a small-sample adjustment that tends to improve the performance of the method when the number of clusters is small.

   For even better small-sample adjustments (Pustejovsky & Tipton, 2018), see the \pkg{clubSandwich} package, which nicely works together with the \pkg{metafor} package (see \sQuote{Examples}).
}
\value{
   An object of class \code{"robust.rma"}. The object is a list containing the following components:
   \item{beta}{estimated coefficients of the model.}
   \item{se}{robust standard errors of the coefficients.}
   \item{zval}{test statistics of the coefficients.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{lower bound of the confidence intervals for the coefficients.}
   \item{ci.ub}{upper bound of the confidence intervals for the coefficients.}
   \item{vb}{robust variance-covariance matrix of the estimated coefficients.}
   \item{QM}{test statistic of the omnibus test of moderators.}
   \item{QMp}{corresponding p-value.}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{\link{print.robust.rma}} function.
}
\note{
   The variable specified via \code{cluster} is assumed to be of the same length as the data originally passed to the \code{rma.uni} or \code{rma.mv} function. Any subsetting and removal of studies with missing values as done when fitting the original model is also automatically applied to the variable specified via \code{cluster}.

   The idea of the robust (sandwich-type) estimator for models with unspecified heteroscedasticity can be traced back to Eicker (1967), Huber (1967), and White (1980). Hence, the method in general is often referred to as the Eicker-Huber-White method. Some small-sample improvements to the method are described by MacKinnon and White (1985). The extension to the cluster robust estimator can be found in Froot (1989) and Williams (2000). Cameron and Miller (2015) provide an extensive overview of cluster robust methods. Sidik and Jonkman (2005, 2006) introduced robust methods in the meta-analytic context for standard random/mixed-effects models. The use of the cluster robust estimator for multivariate/multilevel meta-analytic models is described in Hedges, Tipton, and Johnson (2010).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide to cluster-robust inference. \emph{Journal of Human Resources}, \bold{50}(2), 317--372. \verb{https://doi.org/10.3368/jhr.50.2.317}

   Eicker, F. (1967). Limit theorems for regressions with unequal and dependent errors. In L. M. LeCam & J. Neyman (Eds.), \emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability} (pp. 59--82). Berkeley: University of California Press.

   Froot, K. A. (1989). Consistent covariance matrix estimation with cross-sectional dependence and heteroskedasticity in financial data. \emph{Journal of Financial and Quantitative Analysis}, \bold{24}(3), 333--355. \verb{https://doi.org/10.2307/2330815}

   Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance estimation in meta-regression with dependent effect size estimates. \emph{Research Synthesis Methods}, \bold{1}(1), 39--65. \verb{https://doi.org/10.1002/jrsm.5}

   Huber, P. (1967). The behavior of maximum-likelihood estimates under nonstandard conditions. In L. M. LeCam & J. Neyman (Eds.), \emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability} (pp. 221--233). Berkeley: University of California Press.

   MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. \emph{Journal of Econometrics}, \bold{29}(3), 305--325. \verb{https://doi.org/10.1016/0304-4076(85)90158-7}

   Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. \emph{Journal of Business & Economic Statistics}, \bold{36}(4), 672--683. \verb{https://doi.org/10.1080/07350015.2016.1247004}

   Sidik, K., & Jonkman, J. N. (2005). A note on variance estimation in random effects meta-regression. \emph{Journal of Biopharmaceutical Statistics}, \bold{15}(5), 823--838. \verb{https://doi.org/10.1081/BIP-200067915}

   Sidik, K., & Jonkman, J. N. (2006). Robust variance estimation for random effects meta-analysis. \emph{Computational Statistics & Data Analysis}, \bold{50}(12), 3681--3701. \verb{https://doi.org/10.1016/j.csda.2005.07.019}

   White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. \emph{Econometrica}, \bold{48}(4), 817--838. \verb{https://doi.org/10.2307/1912934}

   Williams, R. L. (2000). A note on robust variance estimation for cluster-correlated data. \emph{Biometrics}, \bold{56}(2), 645--646. \verb{https://doi.org/10.1111/j.0006-341x.2000.00645.x}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mv}}
}
\examples{
### copy data from Konstantopoulos (2011) into 'dat'
dat <- dat.konstantopoulos2011

### multilevel random-effects model
res <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat)
res

### results based on the sandwich method
robust(res, cluster=dat$district)

### illustrate use of the clubSandwich package together with metafor
\dontrun{
require(clubSandwich)
coef_test(res, vcov="CR2", cluster=dat$district)
}

### copy data from Berkey et al. (1998) into 'dat'
dat <- dat.berkey1998

### construct list with the variance-covariance matrices of the observed outcomes for the studies
V <- lapply(split(dat[c("v1i", "v2i")], dat$trial), as.matrix)

### construct block diagonal matrix
V <- bldiag(V)

### fit multivariate model
res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat)
res

### results based on sandwich method
robust(res, cluster=dat$trial)

### illustrate use of the clubSandwich package together with metafor
\dontrun{
require(clubSandwich)
coef_test(res, vcov="CR2", cluster=dat$trial)
}
}
\keyword{htest}
