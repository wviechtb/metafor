\name{coef.matreg}
\alias{coef.matreg}
\alias{vcov.matreg}
\alias{sigma.matreg}
\alias{confint.matreg}
\alias{print.confint.matreg}
\alias{logLik.matreg}
\alias{AIC.matreg}
\alias{BIC.matreg}
\title{Extractor Functions for 'matreg' Objects}
\description{
   Various extractor functions for objects of class \code{"matreg"}. \loadmathjax
}
\usage{
\method{coef}{matreg}(object, \dots)
\method{vcov}{matreg}(object, \dots)
\method{sigma}{matreg}(object, REML=TRUE, \dots)

\method{logLik}{matreg}(object, REML=FALSE, \dots)
\method{AIC}{matreg}(object, \dots, k=2, correct=FALSE, REML=FALSE)
\method{BIC}{matreg}(object, \dots, REML=FALSE)

\method{confint}{matreg}(object, parm, level, digits, \dots)
\method{print}{confint.matreg}(x, digits=x$digits, \dots)
}
\arguments{
   \item{object}{an object of class \code{"matreg"}.}
   \item{REML}{logical whether the returned value should be based on ML or REML estimation.}
   \item{k}{numeric value to specify the penalty per parameter. The default (\code{k=2}) is the classical AIC. See \code{\link{AIC}} for more details.}
   \item{correct}{logical to specify whether the regular (default) or corrected (i.e., AICc) should be extracted.}

   \emph{For \code{confint()}:}
   \item{parm}{this argument is here for compatibility with the generic function \code{\link{confint}}, but is (currently) ignored.}
   \item{level}{numeric value between 0 and 100 to specify the confidence interval level (see \link[=misc-options]{here} for details). If unspecified, the default is to take the value from the object.}
   \item{digits}{optional integer to specify the number of decimal places to which the results should be rounded. If unspecified, the default is to take the value from the object.}
   \item{x}{an object of class \code{"confint.matreg"}.}

   \item{\dots}{other arguments.}
}
\details{
   The \code{coef} function extracts the estimated (standardized) regression coefficients from objects of class \code{"matreg"}. The \code{vcov} function extracts the corresponding variance-covariance matrix (note: the \code{\link{se}} function can also be used to extract the standard errors). The \code{confint} function extracts the confidence intervals.

   Under the \sQuote{Regular \mjseqn{R} Matrix} case (see \code{\link{matreg}}), the \code{sigma} function extracts the square root of the estimated error variance (by default, based on the unbiased estimate of the error variance). The \code{logLik}, \code{AIC}, and \code{BIC} functions extract the corresponding values (note: for compatibility with the behavior for \code{lm} objects, these values are based by default on ML estimation).

}
\value{
   Depending on the function, either a vector, a matrix, or a scalar with the extracted value(s).
}
\author{
   Wolfgang Viechtbauer (\email{wvb@metafor-project.org}, \url{https://www.metafor-project.org}).
}
\references{
   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\seealso{
   \code{\link{matreg}} for the function to create \code{matreg} objects and \code{\link{predict.matreg}} to compute predicted values based on \code{matreg} objects.
}
\examples{
### fit a regression model with lm() to the 'mtcars' dataset
res <- lm(mpg ~ hp + wt + am, data=mtcars)
summary(res)
coef(res)
vcov(res)
se(res)
sigma(res)
confint(res)
logLik(res)
AIC(res)
BIC(res)

### covariance matrix of the dataset
S <- cov(mtcars)

### fit the same regression model using matreg()
res <- matreg(mpg ~ hp + wt + am, R=S, cov=TRUE,
              means=colMeans(mtcars), n=nrow(mtcars))
summary(res)
coef(res)
vcov(res)
se(res)
sigma(res)
confint(res)
logLik(res)
AIC(res)
BIC(res)
}
\keyword{models}
