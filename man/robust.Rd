\name{robust}
\alias{robust}
\alias{robust.rma.uni}
\alias{robust.rma.mv}
\title{Cluster-Robust Tests and Confidence Intervals for 'rma' Objects}
\description{
   The function provides cluster-robust tests and confidence intervals of the model coefficients for objects of class \code{"rma"}. \loadmathjax
}
\usage{
robust(x, cluster, \dots)

\method{robust}{rma.uni}(x, cluster, adjust=TRUE, clubSandwich=FALSE, digits, \dots)
\method{robust}{rma.mv}(x, cluster, adjust=TRUE, clubSandwich=FALSE, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma.uni"} or \code{"rma.mv"}.}
   \item{cluster}{vector to specify the clustering variable to use for constructing the sandwich estimator of the variance-covariance matrix.}
   \item{adjust}{logical to specify whether a small-sample correction should be applied to the variance-covariance matrix.}
   \item{clubSandwich}{logical to specify whether the \href{https://cran.r-project.org/package=clubSandwich}{clubSandwich} package should be used to obtain the cluster-robust tests and confidence intervals.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{\dots}{other arguments.}
}
\details{
   The function constructs a cluster-robust estimate of the variance-covariance matrix of the model coefficients based on a sandwich-type estimator and then computes tests and confidence intervals of the model coefficients.

   By default, tests of individual coefficients and confidence intervals are based on a t-distribution with \mjseqn{n-p} degrees of freedom, while the omnibus test statistic uses an F-distribution with \mjseqn{m} and \mjseqn{n-p} degrees of freedom, where \mjseqn{n} is the number of clusters, \mjseqn{p} denotes the total number of model coefficients (including the intercept if it is present), and \mjseqn{m} denotes the number of coefficients tested by the omnibus test. This is sometimes called the \sQuote{residual} method for approximating the (denominator) degrees of freedom.

   When \code{adjust=TRUE} (the default), the cluster-robust estimate of the variance-covariance matrix is multiplied by the factor \mjseqn{n/(n-p)}, which serves as a small-sample adjustment that tends to improve the performance of the method when the number of clusters is small. This is sometimes called the \sQuote{CR1} adjustment/estimator (in contrast to \sQuote{CR0} when \code{adjust=FALSE}).

   For an even better small-sample adjustment, one can set \code{clubSandwich=TRUE} in which case the \href{https://cran.r-project.org/package=clubSandwich}{clubSandwich} package is used to obtain the cluster-robust tests and confidence intervals. The variance-covariance matrix of the model coefficients is then estimated using the \sQuote{bias-reduced linearization} adjustment proposed by Bell and McCaffrey (2002) and further developed in Tipton (2015) and Pustejovsky and Tipton (2018). This is sometimes called the \sQuote{CR2} adjustment/estimator. The degrees of freedom of the t-tests are then estimated using a Satterthwaite approximation. F-tests are then based on an approximate Hotelling's T-squared reference distribution, with denominator degrees of freedom estimated using a method by Zhang (2012, 2013), as further described in Tipton and Pustejovky (2015).
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

   The results are formatted and printed with the \code{\link{print.rma.uni}} and \code{\link{print.rma.mv}} functions (depending on the type of model).
}
\note{
   The variable specified via \code{cluster} is assumed to be of the same length as the data originally passed to the \code{rma.uni} or \code{rma.mv} functions (and if the \code{data} argument was used in the original model fit, then the variable will be searched for within this data frame first). Any subsetting and removal of studies with missing values that was applied during the model fitting is also automatically applied to the variable specified via the \code{cluster} argument.

   The idea of the robust (sandwich-type) estimator for models with unspecified heteroscedasticity can be traced back to Eicker (1967), Huber (1967), and White (1980, 1984). Hence, the method in general is often referred to as the Eicker-Huber-White method. Some small-sample improvements to the method are described by MacKinnon and White (1985). The extension to the cluster-robust estimator can be found in Froot (1989) and Williams (2000), which is also related to the GEE approach by Liang and Zeger (1986). Cameron and Miller (2015) provide an extensive overview of cluster-robust methods. Sidik and Jonkman (2005, 2006) introduced robust methods in the meta-analytic context for standard random/mixed-effects models. The use of cluster-robust methods for multivariate/multilevel meta-analytic models was introduced by Hedges, Tipton, and Johnson (2010).
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Bell, R. M., & McCaffry, D. F. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. \emph{Survey Methodology}, \bold{28}(2), 169--181. \verb{https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X20020029058}

   Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide to cluster-robust inference. \emph{Journal of Human Resources}, \bold{50}(2), 317--372. \verb{https://doi.org/10.3368/jhr.50.2.317}

   Eicker, F. (1967). Limit theorems for regressions with unequal and dependent errors. In L. M. LeCam & J. Neyman (Eds.), \emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability} (pp. 59--82). Berkeley: University of California Press.

   Froot, K. A. (1989). Consistent covariance matrix estimation with cross-sectional dependence and heteroskedasticity in financial data. \emph{Journal of Financial and Quantitative Analysis}, \bold{24}(3), 333--355. \verb{https://doi.org/10.2307/2330815}

   Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance estimation in meta-regression with dependent effect size estimates. \emph{Research Synthesis Methods}, \bold{1}(1), 39--65. \verb{https://doi.org/10.1002/jrsm.5}

   Huber, P. (1967). The behavior of maximum-likelihood estimates under nonstandard conditions. In L. M. LeCam & J. Neyman (Eds.), \emph{Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability} (pp. 221--233). University of California Press.

   Liang, K. Y., & Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika}, \bold{73}(1), 13--22. \verb{https://doi.org/10.1093/biomet/73.1.13}

   MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. \emph{Journal of Econometrics}, \bold{29}(3), 305--325. \verb{https://doi.org/10.1016/0304-4076(85)90158-7}

   Tipton, E. (2015). Small sample adjustments for robust variance estimation with meta-regression. \emph{Psychological Methods}, \bold{20}(3), 375--393. \verb{https://doi.org/10.1037/met0000011}

   Tipton, E., & Pustejovsky, J. E. (2015). Small-sample adjustments for tests of moderators and model fit using robust variance estimation in meta-regression. \emph{Journal of Educational and Behavioral Statistics}, \bold{40}(6), 604--634. \verb{https://doi.org/10.3102/1076998615606099}

   Sidik, K., & Jonkman, J. N. (2005). A note on variance estimation in random effects meta-regression. \emph{Journal of Biopharmaceutical Statistics}, \bold{15}(5), 823--838. \verb{https://doi.org/10.1081/BIP-200067915}

   Sidik, K., & Jonkman, J. N. (2006). Robust variance estimation for random effects meta-analysis. \emph{Computational Statistics & Data Analysis}, \bold{50}(12), 3681--3701. \verb{https://doi.org/10.1016/j.csda.2005.07.019}

   White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. \emph{Econometrica}, \bold{48}(4), 817--838. \verb{https://doi.org/10.2307/1912934}

   White, H. (1984). \emph{Asymptotic theory for econometricians}. Orlando, FL: Academic Press.

   Williams, R. L. (2000). A note on robust variance estimation for cluster-correlated data. \emph{Biometrics}, \bold{56}(2), 645--646. \verb{https://doi.org/10.1111/j.0006-341x.2000.00645.x}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}

   Zhang, J.-T. (2012). An approximate Hotelling T2-test for heteroscedastic one-way MANOVA. \emph{Open Journal of Statistics}, \bold{2}(1), 1--11. \verb{https://doi.org/10.4236/ojs.2012.21001}

   Zhang, J.-T. (2013). Tests of linear hypotheses in the ANOVA under heteroscedasticity. \emph{International Journal of Advanced Statistics and Probability}, \bold{1}, 9--24. \verb{https://doi.org/10.14419/ijasp.v1i2.908}
}
\seealso{
   \code{\link{rma.uni}} and \code{\link{rma.mv}} for functions to fit models for which cluster-robust tests and confidence intervals can be obtained.
}
\examples{
############################################################################

### copy data from Bangert-Drowns et al. (2004) into 'dat'
dat <- dat.bangertdrowns2004

### fit random-effects model
res <- rma(yi, vi, data=dat)
res

### obtain results based on the sandwich method
robust(res, cluster=id)

### use methods from clubSandwich package
robust(res, cluster=id, clubSandwich=TRUE)

### fit meta-regression model
res <- rma(yi, vi, mods = ~ length, data=dat)
res

### obtain results based on the sandwich method
robust(res, cluster=id)

### use methods from clubSandwich package
robust(res, cluster=id, clubSandwich=TRUE)

############################################################################

### copy data from Konstantopoulos (2011) into 'dat'
dat <- dat.konstantopoulos2011

### fit multilevel random-effects model
res <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat)
res

### obtain results based on the sandwich method
robust(res, cluster=district)

### use methods from clubSandwich package
robust(res, cluster=district, clubSandwich=TRUE)

############################################################################

### copy data from Berkey et al. (1998) into 'dat'
dat <- dat.berkey1998

### variables v1i and v2i correspond to the 2x2 var-cov matrices of the studies;
### so use these variables to construct the V matrix (note: since v1i and v2i are
### var-cov matrices and not correlation matrices, set vi=1 for all rows)
V <- vcalc(vi=1, cluster=author, rvars=c(v1i, v2i), data=dat)

### fit multivariate model
res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat)
res

### obtain results based on sandwich method
robust(res, cluster=trial)

### use methods from clubSandwich package
robust(res, cluster=trial, clubSandwich=TRUE)

############################################################################
}
\keyword{htest}
