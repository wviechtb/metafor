\name{transf}
\alias{transf}
\alias{transf.rtoz}
\alias{transf.ztor}
\alias{transf.logit}
\alias{transf.ilogit}
\alias{transf.arcsin}
\alias{transf.iarcsin}
\alias{transf.pft}
\alias{transf.ipft}
\alias{transf.ipft.hm}
\alias{transf.isqrt}
\alias{transf.irft}
\alias{transf.iirft}
\alias{transf.ahw}
\alias{transf.iahw}
\alias{transf.abt}
\alias{transf.iabt}
\alias{transf.ztor.int}
\alias{transf.exp.int}
\alias{transf.ilogit.int}
\title{Transformation Function}
\description{
   A set of transformation functions useful for meta-analyses.
}
\usage{
transf.rtoz(xi, \dots)
transf.ztor(xi, \dots)
transf.logit(xi, \dots)
transf.ilogit(xi, \dots)
transf.arcsin(xi, \dots)
transf.iarcsin(xi, \dots)
transf.pft(xi, ni, \dots)
transf.ipft(xi, ni, \dots)
transf.ipft.hm(xi, targs, \dots)
transf.isqrt(xi, \dots)
transf.irft(xi, ti, \dots)
transf.iirft(xi, ti, \dots)
transf.ahw(xi, \dots)
transf.iahw(xi, \dots)
transf.abt(xi, \dots)
transf.iabt(xi, \dots)
transf.ztor.int(xi, targs, \dots)
transf.exp.int(xi, targs, \dots)
transf.ilogit.int(xi, targs, \dots)
}
\arguments{
   \item{xi}{vector of values to be transformed.}
   \item{ni}{vector of sample sizes.}
   \item{ti}{vector of person-times at risk.}
   \item{targs}{list with additional arguments for the transformation function. See \sQuote{Details}.}
   \item{\dots}{other arguments.}
}
\details{
   The following transformation functions are currently implemented:
   \itemize{
   \item \code{transf.rtoz}: Fisher's r-to-z transformation for correlations.
   \item \code{transf.ztor}: inverse of the Fisher's r-to-z transformation.
   \item \code{transf.logit}: logit (log odds) transformation for proportions.
   \item \code{transf.ilogit}: inverse of the logit transformation.
   \item \code{transf.arcsin}: arcsine square-root transformation for proportions.
   \item \code{transf.iarcsin}: inverse of the arcsine transformation.
   \item \code{transf.pft}: Freeman-Tukey (double arcsine) transformation for proportions. See Freeman & Tukey (1950). The \code{xi} argument is used to specify the proportions and the \code{ni} argument the corresponding sample sizes.
   \item \code{transf.ipft}: inverse of the Freeman-Tukey (double arcsine) transformation for proportions. See Miller (1978).
   \item \code{transf.ipft.hm}: inverse of the Freeman-Tukey (double arcsine) transformation for proportions using the harmonic mean of the sample sizes for the back-transformation. See Miller (1978). The sample sizes are specified via the \code{targs} argument (the list element should be called \code{ni}).
   \item \code{transf.isqrt}: inverse of the square-root transformation (i.e., function to square a number).
   \item \code{transf.irft}: Freeman-Tukey transformation for incidence rates. See Freeman & Tukey (1950). The \code{xi} argument is used to specify the incidence rates and the \code{ti} argument the corresponding person-times at risk.
   \item \code{transf.iirft}: inverse of the Freeman-Tukey transformation for incidence rates.
   \item \code{transf.ahw}: Transformation of coefficient alpha as suggested by Hakstian & Whalen (1976).
   \item \code{transf.iahw}: Inverse of the transformation of coefficient alpha as suggested by Hakstian & Whalen (1976).
   \item \code{transf.abt}: Transformation of coefficient alpha as suggested by Bonett (2002).
   \item \code{transf.iabt}: Inverse of the transformation of coefficient alpha as suggested by Bonett (2002).
   \item \code{transf.ztor.int}: integral transformation method for the z-to-r transformation.
   \item \code{transf.exp.int}: integral transformation method for the exponential transformation.
   \item \code{transf.ilogit.int}: integral transformation method for the inverse of the logit transformation.
   }
   The integral transformation method for a transformation function \eqn{h(z)} integrates \eqn{h(z) f(z)} over \eqn{z} using the limits \code{targs$lower} and \code{targs$upper}, where \eqn{f(z)} is the density of a normal distribution with mean equal to \code{xi} and variance equal to \code{targs$tau2}. An example is provided below.
}
\value{
   A vector with the transformed values.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \cr
   package website: \url{http://www.metafor-project.org/} \cr
   author homepage: \url{http://www.wvbauer.com/}
}
\references{
   Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient alpha. \emph{Journal of Educational and Behavioral Statistics}, \bold{27}, 335--340.

   Fisher, R. A. (1921). On the \dQuote{probable error} of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32.

   Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular and the square root. \emph{Annals of Mathematical Statistics}, \bold{21}, 607--611.

   Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent alpha coefficients. \emph{Psychometrika}, \bold{41}, 219--231.

   Miller, J. J. (1978). The inverse of the Freeman-Tukey double arcsine transformation. \emph{American Statistician}, \bold{32}, 138.

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \url{http://www.jstatsoft.org/v36/i03/}.
}
\examples{
### meta-analysis of the log relative risks using a random-effects model
res <- rma(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### average relative risk with 95\% CI
predict(res, transf=exp)

### average relative risk with 95\% CI using the integral transformation
predict(res, transf=transf.exp.int, targs=list(tau2=res$tau2, lower=-4, upper=4))
}
\keyword{manip}
