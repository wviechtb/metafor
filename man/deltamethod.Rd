\name{deltamethod}
\alias{deltamethod}
\title{Apply the (Multivariate) Delta Method}
\description{
   Function to apply the (multivariate) delta method to a set of estimates. \loadmathjax
}
\usage{
deltamethod(x, vcov, fun, level, H0=0, digits)
}
\arguments{
   \item{x}{either a vector of estimates or a model object from which model coefficients can be extracted via \code{coef(x)}.}
   \item{vcov}{when \code{x} is a vector of estimates, the corresponding variance-covariance matrix (ignored when \code{x} is a model object, in which case \code{vcov(x)} is used to extract the variance-covariance matrix).}
   \item{fun}{a function to apply to the estimates.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (see \link[=misc-options]{here} for details). If unspecified, this either defaults to 95 or, if possible, the corresponding value from the model object.}
   \item{H0}{numeric value to specify the value under the null hypothesis for the Wald-type test(s) (the default is 0). Can also be a vector.}
   \item{digits}{optional integer to specify the number of decimal places to which the printed results should be rounded.}
}
\details{
   Let \mjeqn{\hat{\theta}}{\theta} denote a vector of \mjseqn{p} estimates which can be specified via the \code{x} argument and let \mjseqn{\Sigma} denote the corresponding \mjeqn{p \times p}{pxp} variance-covariance matrix, which can be specified via the \code{vcov} argument. If \code{x} is not an vector with estimates, then the function assumes that \code{x} is a model object and will try to use \code{coef(x)} and \code{vcov(x)} to extract the model coefficients and the corresponding variance-covariance matrix (in this case, the \code{vcov} argument is ignored).

   Let \mjeqn{f(\cdot)}{f(.)} be a function, specified via the \code{fun} argument, with \mjseqn{p} inputs/arguments (or with a single argument that is assumed to be a vector of length \mjseqn{p}), which returns a numeric (and atomic) vector of \mjseqn{q} transformed estimates. Then the function computes \mjeqn{f(\hat{\theta})}{f(\theta)} and the corresponding variance-covariance matrix of the transformed estimates using the \href{https://en.wikipedia.org/wiki/Delta_method#Multivariate_delta_method}{multivariate delta method} (e.g., van der Vaart, 1998) with \mjdeqn{\textrm{Var}[f(\hat{\theta})] = \nabla f(\hat{\theta})' \cdot \Sigma \cdot \nabla f(\hat{\theta})}{Var[f(\theta)]  = ▽f(\theta)' \Sigma ▽f(\theta)} where \mjeqn{\nabla f(\hat{\theta})}{▽f(\theta)} denotes the gradient of \mjeqn{f(\cdot)}{f(.)} evaluated at \mjeqn{\hat{\theta}}{\theta}. The function computes the gradient numerically using the \code{\link[calculus]{derivative}} function from the \code{calculus} package.

   The function also computes Wald-type tests and confidence intervals for the \mjseqn{q} transformed estimates. The \code{level} argument can be used to control the confidence interval level.
}
\value{
   An object of class \code{"deltamethod"}. The object is a list containing the following components:
   \item{tab}{a data frame with the transformed estimates, standard errors, test statistics, p-values, and lower/upper confidence interval bounds.}
   \item{vcov}{the variance-covariance matrix of the transformed estimates.}
   \item{\dots}{some additional elements/values.}

   The results are formatted and printed with the \code{\link[=print.deltamethod]{print}} function. Extractor functions include \code{\link[=coef.deltamethod]{coef}} and \code{\link[=vcov.deltamethod]{vcov}}.
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   van der Vaart, A. W. (1998). \emph{Asymptotic statistics}. Cambridge, UK: Cambridge University Press.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{conv.delta}} for a function to apply the (univariate) delta method to observed effect sizes or outcomes and their sampling variances.
}
\examples{
############################################################################

### copy data into 'dat'
dat <- dat.craft2003

### construct dataset and var-cov matrix of the correlations
tmp <- rcalc(ri ~ var1 + var2 | study, ni=ni, data=dat)
V <- tmp$V
dat <- tmp$dat

### turn var1.var2 into a factor with the desired order of levels
dat$var1.var2 <- factor(dat$var1.var2,
   levels=c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf"))

### multivariate random-effects model
res <- rma.mv(yi, V, mods = ~ 0 + var1.var2, random = ~ var1.var2 | study, struct="UN", data=dat)
res

### restructure estimated mean correlations into a 4x4 matrix
R <- vec2mat(coef(res))
rownames(R) <- colnames(R) <- c("perf", "acog", "asom", "conf")
round(R, digits=3)

### check that order in vcov(res) corresponds to order in R
round(vcov(res), digits=4)

### fit regression model with 'perf' as outcome and 'acog', 'asom', and 'conf' as predictors
matreg(1, 2:4, R=R, V=vcov(res))

### same analysis but using the deltamethod() function
deltamethod(coef(res), vcov(res), fun=function(r1,r2,r3,r4,r5,r6) {
   R <- vec2mat(c(r1,r2,r3,r4,r5,r6))
   setNames(c(solve(R[-1,-1]) \%*\% R[2:4,1]), c("acog","asom","conf"))
})

### using a function that takes a vector as input
deltamethod(coef(res), vcov(res), fun=function(r) {
   R <- vec2mat(r)
   setNames(c(solve(R[-1,-1]) \%*\% R[2:4,1]), c("acog","asom","conf"))
})

############################################################################

### construct dataset and var-cov matrix of the r-to-z transformed correlations
dat <- dat.craft2003
tmp <- rcalc(ri ~ var1 + var2 | study, ni=ni, data=dat, rtoz=TRUE)
V <- tmp$V
dat <- tmp$dat

### turn var1.var2 into a factor with the desired order of levels
dat$var1.var2 <- factor(dat$var1.var2,
   levels=c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf"))

### multivariate random-effects model
res <- rma.mv(yi, V, mods = ~ 0 + var1.var2, random = ~ var1.var2 | study, struct="UN", data=dat)
res

### estimate the difference between r(acog,perf) and r(asom,perf)
deltamethod(res, fun=function(z1,z2,z3,z4,z5,z6) {
   transf.ztor(z1) - transf.ztor(z2)
})

### using a function that takes a vector as input
deltamethod(res, fun=function(z) {
   transf.ztor(z[1]) - transf.ztor(z[2])
})

############################################################################
}
\keyword{models}
