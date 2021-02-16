\name{vif}
\alias{vif}
\alias{vif.rma}
\alias{print.vif.rma}
\title{Variance Inflation Factors for 'rma' Objects}
\description{
   Compute variance inflation factors (VIFs) for objects of class \code{"rma"}. \loadmathjax
}
\usage{
vif(x, \dots)

\method{vif}{rma}(x, btt, intercept=FALSE, table=FALSE, digits, \dots)

\method{print}{vif.rma}(x, digits=x$digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"} (for \code{vif}) or \code{"vif.rma"} (for \code{print}).}
   \item{btt}{optional vector of indices to specify a set of coefficients for which a generalized variance inflation factor (GVIF) should be computed. Can also be a string to grep for. See \sQuote{Details}.}
   \item{intercept}{logical to specify whether to include the intercept (if the model includes one) in the computation of the VIFs (the default is \code{FALSE}). See \sQuote{Note}.}
   \item{table}{logical to specify whether the VIFs should be added to the model coefficient table (the default is \code{FALSE}).}
   \item{digits}{integer to specify the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The function computes variance inflation factors (VIFs) for meta-regression models. Hence, the model specified via argument \code{x} must include moderator variables (and more than one for this to be useful, as the VIF for a model with a single moderator variable will always be equal to 1).

   Let \mjseqn{b_j} denote the estimate of the \mjseqn{j}th model coefficient of a particular meta-regression model and \mjeqn{\mbox{Var}[b_j]}{Var[b_j]} its variance (i.e., the corresponding diagonal element from the matrix obtained with the \code{\link{vcov.rma}} function). Moreover, let \mjseqn{b'_j} denote the estimate of the same model coefficient if the other moderator variables in the model had \emph{not} been included in the model and \mjeqn{\mbox{Var}[b'_j]}{Var[b'_j]} the corresponding variance. Then the VIF for the model coefficient is given by \mjdeqn{\mbox{VIF}[b_j] = \frac{\mbox{Var}[b_j]}{\mbox{Var}[b'_j]},}{VIF[b_j] = Var[b_j] / Var[b'_j],} which indicates the inflation in the variance of the estimated model coefficient due to potential collinearity of the \mjseqn{j}th moderator variable with the other moderator variables in the model. Taking the square root of a VIF gives the corresponding standard error inflation factor (SIF).

   If \code{btt} is not specified, then the VIF is computed for each individual model coefficient. However, if the model includes factors (coded in terms of multiple dummy variables) or other sets of moderator variables that belong together (e.g., for polynomials or cubic splines), one may want to examine how much the variance in all of the coefficients in the set is jointly impacted by collinearity with the other moderator variables in the model. For this, we can compute a generalized variance inflation factor (GVIF) (Fox & Monette, 1992) by setting the \code{btt} argument equal to the indices of those coefficients for which the GVIF should be computed. The square root of a GVIF indicates the inflation in the confidence ellipse/(hyper)ellipsoid for the set of coefficients corresponding to the set due to collinearity. However, to make this value more directly comparable to SIFs (based on single coefficients) or when the set includes a different number of coefficients, the function computes the generalized standard error inflation factor (GSIF) by raising the GVIF to the power of \mjseqn{1/(2m)} (where \mjseqn{m} denotes the number of coefficients in the set).
}
\value{
   If \code{btt} is not specified, either a vector (if \code{table=FALSE}) with the VIFs or a data frame (if \code{table=TRUE}) with the following elements:
   \item{estimate}{estimated model coefficients.}
   \item{se}{corresponding standard errors.}
   \item{zval}{corresponding test statistics.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{corresponding lower bound of the confidence intervals.}
   \item{ci.ub}{corresponding upper bound of the confidence intervals.}
   \item{vif}{corresponding variance inflation factors.}
   \item{sif}{corresponding standard error inflation factors.}
   If \code{btt} is specified, a list with elements \code{gvif} and \code{gsif} with the GVIF and GSIF values for the set of coefficients specified.
}
\note{
   The values of the (G)VIFs are invariant to the scaling of the predictor variables if the model includes an intercept that is removed when inverting the correlation matrix of the model coefficients to compute the (G)VIFs. This is the default behavior. See \sQuote{Examples}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). \emph{Regression diagnostics}. New York: Wiley.

   Fox, J., & Monette, G. (1992). Generalized collinearity diagnostics. \emph{Journal of the American Statistical Association}, \bold{87}(417), 178-183. \verb{https://doi.org/10.2307/2290467}

   Stewart, G. W. (1987). Collinearity and least squares regression. \emph{Statistical Science}, \bold{2}(1), 68-84. \verb{https://doi.org/10.1214/ss/1177013439}

   Wax, Y. (1992). Collinearity diagnosis for a relative risk regression-analysis: An application to assessment of diet cancer relationship in epidemiologic studies. \emph{Statistics in Medicine}, \bold{11}(10), 1273--1287. \verb{https://doi.org/10.1002/sim.4780111003}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{rma.uni}}, \code{\link{rma.mv}}, \code{\link{rma.glmm}}
}
\examples{
### copy data from Bangert-Drowns et al. (2004) into 'dat'
dat <- dat.bangertdrowns2004

### fit mixed-effects meta-regression model
res <- rma(yi, vi, mods = ~ length + wic + feedback + info + pers + imag + meta, data=dat)

### get variance inflation factors
vif(res)

### show that VIFs are not influenced by scaling of the predictors
u <- scale # to standardize the predictors
res <- rma(yi, vi, mods = ~ u(length) + u(wic) + u(feedback) + u(info) +
                            u(pers) + u(imag) + u(meta), data=dat)
vif(res)

### get full table
vif(res, table=TRUE)

### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit meta-regression model where one predictor (alloc) is a three-level factor
res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)

### get variance inflation factors for all individual coefficients
vif(res, table=TRUE)

### generalized variance inflation factor for the 'alloc' factor
vif(res, btt=3:4)

### can also specify a string to grep for
vif(res, btt="alloc")
}
\keyword{models}
