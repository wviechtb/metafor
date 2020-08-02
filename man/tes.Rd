\name{tes}
\alias{tes}
\alias{tes.rma}
\alias{tes.default}
\alias{print.tes}
\title{Test of Excess Significance}
\description{
   Function to conduct the test of excess significance. \loadmathjax
}
\usage{
tes(x, \dots)

\method{tes}{rma}(x, H0=0, alternative="two.sided", alpha=.05,
    test, tes.alternative="greater", progbar=TRUE, tes.alpha=.10, digits, \dots)

\method{tes}{default}(x, vi, sei, subset, H0=0, alternative="two.sided", alpha=.05, theta, tau2,
    test, tes.alternative="greater", progbar=TRUE, tes.alpha=.10, digits, \dots)

\method{print}{tes}(x, digits=x$digits, \dots)
}
\arguments{
   \emph{These arguments pertain to data input:}

   \item{x}{an object of class \code{"rma"} or a vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{sei}{vector with the corresponding standard errors (note: only one of the two, \code{vi} or \code{sei}, needs to be specified).}
   \item{subset}{optional (logical or numeric) vector indicating the subset of studies that should be included.}

   \emph{These arguments pertain to the tests of the observed effect sizes or outcomes:}

   \item{H0}{numeric value specifying the value of the effect size or outcome under the null hypothesis (the default is 0).}
   \item{alternative}{character string indicating the sidedness of the hypothesis when testing the observed effect sizes or outcomes. Possible options are \code{"two.sided"} (the default), \code{"greater"}, or \code{"less"}. Can be abbreviated.}
   \item{alpha}{alpha level for testing the observed effect sizes or outcomes (the default is .05).}

   \emph{These arguments pertain to the power of the tests:}

   \item{theta}{numeric value specifying the value of the true effect size or outcome under the alternative hypothesis. If unspecified, it will be estimated based on the data or the value is taken from the \code{"rma"} object.}
   \item{tau2}{numeric value specifying the amount of heterogeneity in the true effects or outcomes. If unspecified, the true effect sizes or outcomes are assumed to be homogeneous or the value is taken from the \code{"rma"} object.}

   \emph{These arguments pertain to the test of excess significance:}

   \item{test}{character string indicating the type of test to use for conducting the test of excess significance. Possible options are \code{"chi2"}, \code{"binom"}, or \code{"exact"}. Can be abbreviated. If unspecified, the function chooses the type of test based on the data.}
   \item{tes.alternative}{character string indicating the sidedness of the hypothesis for the test of excess significance. Possible options are \code{"greater"} (the default), \code{"two.sided"}, or \code{"less"}. Can be abbreviated.}
   \item{progbar}{logical indicating whether a progress bar should be shown (the default is \code{TRUE}). Only relevant when conducting an exact test.}
   \item{tes.alpha}{alpha level for the test of excess significance (the default is .10). Only relevant for finding the \sQuote{limit estimate}.}

   \emph{Miscellaneous arguments:}

   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is 4).}
   \item{\dots}{other arguments.}
}
\details{
   The function carries out the test of excess significance described by Ioannidis and Trikalinos (2007). The test can be used to examine whether the observed number of significant findings is greater than the number of significant findings expected given the power of the tests. An overabundance of significant tests may suggest that the collection of studies is not representative of all studies conducted on a particular topic.

   One can either pass an object of class \code{"rma"} to the function or a vector with the observed effect sizes or outcomes (via \code{x}) and the corresponding sampling variances via \code{vi} (or the standard errors via \code{sei}).

   The observed effect sizes or outcomes are tested for significance based on a standard Wald-type test, that is, by comparing \mjdeqn{z_i = \frac{y_i - \mbox{H}_0}{\sqrt{v_i}}}{z_i = (y_i - H_0) / sqrt(v_i)} against the appropriate critical value(s) of a standard normal distribution (e.g., \mjseqn{\pm 1.96} for \code{alternative="two.sided"} and \code{alpha=.05}, which are the defaults). Let \mjseqn{O} denote the observed number of significant tests.

   Given a particular value for the true effect or outcome denoted by \mjseqn{\theta} (which, if it is unspecified, is determined by computing the inverse-variance weighted average of the observed effect sizes or outcomes or the value is taken from the model object), let \mjseqn{1-\beta_i} denote the power of the \mjseqn{i}th test (where \mjseqn{\beta_i} denotes the Type II error probability). If \mjseqn{\tau^2 > 0}, let \mjseqn{1-\beta_i} denote the expected power (computed based on integrating the power over a normal distribution with mean \mjseqn{\theta} and variance \mjseqn{\tau^2}). Let \mjseqn{E = \sum_{i=1}^k (1-\beta_i)} denote the expected number of significant tests.

   The test of excess significance then tests if \mjseqn{O} is significantly greater (if \code{tes.alternative="greater"}) than \mjseqn{E}. This can be done using Pearson's chi-squared test (if \code{test="chi2"}), a binomial test (if \code{test="binomial"}), or an exact test (if \code{test="exact"}). The latter is described in Francis (2013). If argument \code{test} is unspecified, the default is to do an exact test if the number of elements in the sum that needs to be computed is less than or equal to \code{10^6} and to do a chi-square test otherwise.

   One can also iteratively find the value of \mjseqn{\theta} such that the p-value of the test of excess significance is equal to \code{tes.alpha} (which is .10 by default). The resulting value is called the \sQuote{limit estimate} and is denoted \mjeqn{\theta_{lim}}{\theta_lim} by Ioannidis and Trikalinos (2007). Note that the limit estimate is not computable if the p-value is larger than \code{tes.alpha} even if \mjeqn{\theta = \mbox{H}_0}{\theta = H_0}.
}
\value{
   An object of class \code{"tes"}. The object is a list containing the following components:
   \item{k}{the number of studies included in the analysis.}
   \item{O}{the observed number of significant tests.}
   \item{E}{the expected number of significant tests.}
   \item{OEratio}{the ratio of O over E.}
   \item{test}{the type of test conducted.}
   \item{pval}{the p-value of the test of excess significance.}
   \item{power}{the (estimated) power of the tests.}
   \item{sig}{logical vector indicating which tests were significant.}
   \item{theta}{the value of \mjseqn{\theta} used for computing the power of the tests.}
   \item{theta.lim}{the \sQuote{limit estimate} (i.e., \mjeqn{\theta_{lim}}{\theta_lim}).}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{print.tes} function.
}
\note{
   When \code{tes.alternative="greater"} (the default), then the function tests if \mjseqn{O} is significantly greater than \mjseqn{E} and hence this is indeed a test of excess significance. When \code{tes.alternative="two.sided"}, then the function tests if \mjseqn{O} differs significantly from \mjseqn{E} in either direction and hence it would be more apt to describe this as a test of (in)consistency (between \mjseqn{O} and \mjseqn{E}). Finally, one can also set \code{tes.alternative="less"}, in which case the function tests if \mjseqn{O} is significantly lower than \mjseqn{E}, which could be considered a test of excess non-significance.

   When \code{tes.alternative="two.sided"}, one can actually compute two limit estimates. The function attempts to compute both.

   The function computes the significance and power of the studies based on Wald-type tests regardless of the effect size or outcome measure used as input. This works as an adequate approximation as long as the within-study sample sizes are not too small.

   Note that the test is not a test for publication bias but a test whether the set of studies includes an unusual number of significant findings given the power of the studies. The general usefulness of the test and its usefulness under particular circumstances (e.g., when there is substantial heterogeneity in the true effects or outcomes) has been the subject of considerable debate. See Francis (2013) and the commentaries on this article in the same issue of the journal.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org}
}
\references{
   Francis, G. (2013). Replication, statistical consistency, and publication bias. \emph{Journal of Mathematical Psychology}, \bold{57}(5), 153--169. \verb{https://doi.org/10.1016/j.jmp.2013.02.003}

   Ioannidis, J. P. A., & Trikalinos, T. A. (2007). An exploratory test for an excess of significant findings. \emph{Clinical Trials}, \bold{4}(3), 245--253. \verb{https://doi.org/10.1177/1740774507079441}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{ranktest}}, \code{\link{regtest}}, \code{\link{trimfill}}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=x.a, n1i=n.a, ci=x.p, n2i=n.p, data=dat.dorn2007)

### conduct test of excess significance (using test="chi2" to speed things up)
tes(dat$yi, dat$vi, test="chi2")

### same as fitting a FE model and then passing the object to the function
res <- rma(yi, vi, data=dat, method="FE")
tes(res, test="chi2")

### illustrate limit estimate (value of theta where p-value of test is equal to tes.alpha)
thetas <- seq(0,1,length=101)
pvals <- sapply(thetas, function(theta) tes(dat$yi, dat$vi, test="chi2", theta=theta)$pval)
plot(thetas, pvals, type="o", pch=19, ylim=c(0,1))
sav <- tes(dat$yi, dat$vi, test="chi2")
abline(h=sav$tes.alpha, lty="dotted")
abline(v=sav$theta.lim, lty="dotted")

### examine significance of test as a function of alpha (to examine 'significance chasing')
alphas <- seq(.01,.99,length=101)
pvals <- sapply(alphas, function(alpha) tes(dat$yi, dat$vi, test="chi2", alpha=alpha)$pval)
plot(alphas, pvals, type="o", pch=19, ylim=c(0,1))
abline(v=.05, lty="dotted")
abline(h=.10, lty="dotted")
}
\keyword{htest}
