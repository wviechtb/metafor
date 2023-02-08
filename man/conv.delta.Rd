\name{conv.delta}
\alias{conv.delta}
\title{Transform Observed Effect Sizes or Outcomes and their Sampling Variances using the Delta Method}
\description{
   Function to transform observed effect sizes or outcomes and their sampling variances using the delta method. \loadmathjax
}
\usage{
conv.delta(yi, vi, data, include, transf, var.names, append=TRUE, replace="ifna", \dots)
}
\arguments{
   \item{yi}{vector with the observed effect sizes or outcomes.}
   \item{vi}{vector with the corresponding sampling variances.}
   \item{data}{optional data frame containing the variables given to the arguments above.}
   \item{include}{optional (logical or numeric) vector to specify the subset of studies for which the transformation should be carried out.}
   \item{transf}{a function which should be used for the transformation.}
   \item{var.names}{character vector with two elements to specify the name of the variable for the transformed effect sizes or outcomes and the name of the variable for the corresponding sampling variances (if \code{data} is an object of class \code{"escalc"}, the \code{var.names} are taken from the object; otherwise the defaults are \code{"yi"} and \code{"vi"}).}
   \item{append}{logical to specify whether the data frame provided via the \code{data} argument should be returned together with the estimated values (the default is \code{TRUE}).}
   \item{replace}{character string to specify how values in \code{var.names} should be replaced (only relevant when using the \code{data} argument and if variables \code{var.names} already exist in the data frame). See the \sQuote{Value} section for more details.}
   \item{\dots}{other arguments for the transformation function.}
}
\details{
   The \code{\link{escalc}} function can be used to compute a wide variety of effect sizes or \sQuote{outcome measures}. In some cases, it may be necessary to transform one type of measure to another. The present function provides a general method for doing so via the \href{https://en.wikipedia.org/wiki/Delta_method}{delta method}, which briefly works as follows.

   Let \mjseqn{y_i} denote the observed effect size or outcome for a particular study and \mjseqn{v_i} the corresponding sampling variance. Then \mjseqn{f(y_i)} will be the transformed effect size or outcome, where \mjeqn{f(\cdot)}{f(.)} is the function specified via the \code{transf} argument. The sampling variance of the transformed effect size or outcome is then computed with \mjseqn{v_i \times f'(y_i)^2}, where \mjseqn{f'(y_i)} denotes the derivative of \mjeqn{f(\cdot)}{f(.)} evaluated at \mjseqn{y_i}. The present function computes the derivative numerically using the \code{\link[numDeriv]{grad}} function from the \code{numDeriv} package.

   The value of the observed effect size or outcome should be the first argument of the function specified via \code{transf}. The function can have additional arguments, which can be specified via the \dots argument. However, due to the manner in which these additional arguments are evaluated, they cannot have names that match one of the arguments of the \code{\link[numDeriv]{grad}} function (an error will be issed if such a naming clash is detected).
}
\value{
   If the \code{data} argument was not specified or \code{append=FALSE}, a data frame of class \code{c("escalc","data.frame")} with two variables called \code{var.names[1]} (by default \code{"yi"}) and \code{var.names[2]} (by default \code{"vi"}) with the transformed observed effect sizes or outcomes and the corresponding sampling variances (computed as described above).

   If \code{data} was specified and \code{append=TRUE}, then the original data frame is returned. If \code{var.names[1]} is a variable in \code{data} and \code{replace="ifna"}, then missing values in this variable are replaced with the transformed observed effect sizes or outcomes (where possible) and otherwise a new variable called \code{var.names[1]} is added to the data frame. Similarly, if \code{var.names[2]} is a variable in \code{data} and \code{replace="ifna"}, then missing values in this variable are replaced with the sampling variances calculated as described above (where possible) and otherwise a new variable called \code{var.names[2]} is added to the data frame.

   If \code{replace="all"}, then all values in \code{var.names[1]} and \code{var.names[2]} are replaced, even for cases where the value in \code{var.names[1]} and \code{var.names[2]} is not missing.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{escalc}} for a function to compute various effect size measures.
}
\examples{
############################################################################

### the following examples illustrate that the use of the delta method (with numeric derivatives)
### yields essentially identical results as the analytic calculations that are done by escalc()

### compute logit transformed proportions and corresponding sampling variances for two studies
escalc(measure="PLO", xi=c(5,12), ni=c(40,80))

### compute raw proportions and corresponding sampling variances for two studies
dat <- escalc(measure="PR", xi=c(5,12), ni=c(40,80))
dat

### apply the logit transformation (note: this yields the same values as above with measure="PLO")
conv.delta(dat$yi, dat$vi, transf=transf.logit)

### using the 'data' argument
conv.delta(yi, vi, data=dat, transf=transf.logit, var.names=c("yi.t","vi.t"))

### or replace the existing 'yi' and 'vi' values
conv.delta(yi, vi, data=dat, transf=transf.logit, replace="all")

######################################

### use escalc() with measure D2ORN which transforms the standardized mean differences (computed
### from the means and standard deviations) into the corresponding log odds ratios
escalc(measure="D2ORN", m1i=m1i, sd1i=sd1i, n1i=n1i,
                        m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat.normand1999)

### use escalc() to compute standardized mean differences (without the usual bias correction) and
### then apply the same transformation to the standardized mean differences
dat <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i,
                             m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat.normand1999, correct=FALSE)
conv.delta(yi, vi, data=dat, transf=transf.d2orn, replace="all")

######################################

### an example where the transformation function takes additional arguments

### use escalc() with measure RPB which transforms the standardized mean differences (computed
### from the means and standard deviations) into the corresponding point-biserial correlations
escalc(measure="RPB", m1i=m1i, sd1i=sd1i, n1i=n1i,
                      m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat.normand1999)

### use escalc() to compute standardized mean differences (without the usual bias correction) and
### then apply the same transformation to the standardized mean differences
dat <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i,
                             m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat.normand1999, correct=FALSE)
conv.delta(yi, vi, data=dat,
           transf=function(di, n1i, n2i) di / sqrt(di^2 + (n1i+n2i-2)/n1i + (n1i+n2i-2)/n2i),
           n1i=n1i, n2i=n2i, replace="all")

############################################################################

### a more elaborate example showing how this function could be used in the data
### preparation steps for a meta-analysis of standardized mean differences (SMDs)

dat <- data.frame(study=1:4,
                  m1i=c(2.03,NA,NA,NA), sd1i=c(0.95,NA,NA,NA), n1i=c(32,95,NA,NA),
                  m2i=c(1.25,NA,NA,NA), sd2i=c(1.04,NA,NA,NA), n2i=c(30,99,NA,NA),
                  ti=c(NA,2.12,NA,NA),
                  ai=c(NA,NA,26,NA), bi=c(NA,NA,58,NA), ci=c(NA,NA,11,NA), di=c(NA,NA,74,NA),
                  or=c(NA,NA,NA,2.56), lower=c(NA,NA,NA,1.23), upper=c(NA,NA,NA,5.30))
dat

### study types:
### 1) provides means and SDs so that the SMD can be directly calculated
### 2) provides the t-statistic from an independent samples t-test
### 3) dichotomized the continuous dependent variable and reports the resulting 2x2 table
### 4) dichotomized the continuous dependent variable and reports an odds ratio (or) with 95% CI

### use escalc() to directly compute the SMD and its variance for studies 1 and 2
dat <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i,
                             m2i=m2i, sd2i=sd2i, n2i=n2i, ti=ti, data=dat)
dat

### use escalc() with measure OR2DN to compute the SMD value for study 3
dat <- escalc(measure="OR2DN", ai=ai, bi=bi, ci=ci, di=di, data=dat, replace=FALSE)
dat

### use conv.wald() to convert the OR and CI into the log odds ratio and its variance for study 4
dat <- conv.wald(out=or, ci.lb=lower, ci.ub=upper, data=dat,
                 transf=log, var.names=c("lnor","vlnor"))
dat

### use conv.delta() to transform the log odds ratio into the SMD value for study 4
dat <- conv.delta(lnor, vlnor, data=dat,
                  transf=transf.or2dn, var.names=c("yi","vi"))
dat

### fit an equal-effects model to the SMD values
rma(yi, vi, data=dat, method="EE")

############################################################################
}
\keyword{manip}
