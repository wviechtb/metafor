\name{vif}
\alias{vif}
\alias{vif.rma}
\title{Variance Inflation Factors for 'rma' Objects}
\description{Compute variance inflation factors (VIFs) for objects of class \code{"rma"}.}
\usage{
vif(x, \dots)

\method{vif}{rma}(x, intercept=FALSE, table=FALSE, digits, \dots)
}
\arguments{
   \item{x}{an object of class \code{"rma"}.}
   \item{intercept}{logical indicating whether to include the intercept (if the model includes one) in the computation of the VIFs (the default is \code{FALSE}). See \sQuote{Note}.}
   \item{table}{logical indicating whether the VIFs should be added to the model coefficient table (the default is \code{FALSE}).}
   \item{digits}{integer specifying the number of decimal places to which the printed results should be rounded (if unspecified, the default is to take the value from the object).}
   \item{\dots}{other arguments.}
}
\details{
   The function computes variance inflation factors (VIFs) for meta-regression models. Hence, the model specified via \code{x} must include moderator variables. VIFs indicate the inflation in the variance of a model coefficient due to collinearity among the predictor variables.
}
\value{
   Either a vector (if \code{table=FALSE}) with the VIFs or a data frame (if \code{table=TRUE}) with the following elements:
   \item{estimate}{estimated model coefficients.}
   \item{se}{corresponding standard errors.}
   \item{zval}{corresponding test statistics.}
   \item{pval}{corresponding p-values.}
   \item{ci.lb}{corresponding lower bound of the confidence intervals.}
   \item{ci.ub}{corresponding upper bound of the confidence intervals.}
   \item{vif}{corresponding variance inflation factors.}
}
\note{
   The values of the VIFs are invariant to the scaling of the predictor variables if the model includes an intercept that is removed when inverting the correlation matrix of the model coefficients to compute the VIFs. This is the default behavior. See \sQuote{Examples}.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{http://www.metafor-project.org/}
}
\references{
   Davis, C. E., Hyde, J. E., Bangdiwala, S. I., & Nelson, J. J. (1986). An example of dependencies among variables in a conditional logistic regression. In S. H. Moolgavkar & R. L. Prentice (Eds.), \emph{Modern statistical methods in chronic disease epidemiology} (pp. 140--147). New York: Wiley.

   Wax, Y. (1992). Collinearity diagnosis for a relative risk regression-analysis: An application to assessment of diet cancer relationship in epidemiologic studies. \emph{Statistics in Medicine}, \bold{11}, 1273--1287.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
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
}
\keyword{models}
