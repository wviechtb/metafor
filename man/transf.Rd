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
\alias{transf.dtou1}
\alias{transf.dtou2}
\alias{transf.dtou3}
\alias{transf.dtocles}
\alias{transf.dtorpb}
\alias{transf.dtobesd}
\alias{transf.dtomd}
\alias{transf.logortord}
\alias{transf.logortorr}
\title{Transformation Functions}
\description{
   A set of transformation functions useful for meta-analyses. \loadmathjax
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
transf.dtou1(xi, \dots)
transf.dtou2(xi, \dots)
transf.dtou3(xi, \dots)
transf.dtorpb(xi, n1i, n2i, \dots)
transf.dtobesd(xi, \dots)
transf.dtomd(xi, targs, \dots)
transf.logortord(xi, pc, \dots)
transf.logortorr(xi, pc, \dots)
}
\arguments{
   \item{xi}{vector of values to be transformed.}
   \item{ni}{vector of sample sizes.}
   \item{n1i}{vector of sample sizes for the first group.}
   \item{n2i}{vector of sample sizes for the second group.}
   \item{ti}{vector of person-times at risk.}
   \item{pc}{control group risk (either a single value or a vector).}
   \item{targs}{list with additional arguments for the transformation function. See \sQuote{Details}.}
   \item{\dots}{other arguments.}
}
\details{
   The following transformation functions are currently implemented:
   \itemize{
   \item \code{transf.rtoz}: Fisher's r-to-z transformation for correlations.
   \item \code{transf.ztor}: inverse of Fisher's r-to-z transformation.
   \item \code{transf.logit}: logit (log odds) transformation for proportions.
   \item \code{transf.ilogit}: inverse of the logit transformation.
   \item \code{transf.arcsin}: arcsine square root transformation for proportions.
   \item \code{transf.iarcsin}: inverse of the arcsine transformation.
   \item \code{transf.pft}: Freeman-Tukey (double arcsine) transformation for proportions. See Freeman & Tukey (1950). The \code{xi} argument is used to specify the proportions and the \code{ni} argument the corresponding sample sizes.
   \item \code{transf.ipft}: inverse of the Freeman-Tukey (double arcsine) transformation for proportions. See Miller (1978).
   \item \code{transf.ipft.hm}: inverse of the Freeman-Tukey (double arcsine) transformation for proportions using the harmonic mean of the sample sizes for the back-transformation. See Miller (1978). The sample sizes are specified via the \code{targs} argument (the list element should be called \code{ni}).
   \item \code{transf.isqrt}: inverse of the square root transformation (i.e., function to square a number).
   \item \code{transf.irft}: Freeman-Tukey transformation for incidence rates. See Freeman & Tukey (1950). The \code{xi} argument is used to specify the incidence rates and the \code{ti} argument the corresponding person-times at risk.
   \item \code{transf.iirft}: inverse of the Freeman-Tukey transformation for incidence rates.
   \item \code{transf.ahw}: transformation of coefficient alpha as suggested by Hakstian & Whalen (1976).
   \item \code{transf.iahw}: inverse of the transformation of coefficient alpha as suggested by Hakstian & Whalen (1976).
   \item \code{transf.abt}: transformation of coefficient alpha as suggested by Bonett (2002).
   \item \code{transf.iabt}: inverse of the transformation of coefficient alpha as suggested by Bonett (2002).
   \item \code{transf.ztor.int}: integral transformation method for the z-to-r transformation.
   \item \code{transf.exp.int}: integral transformation method for the exponential transformation.
   \item \code{transf.ilogit.int}: integral transformation method for the inverse of the logit transformation.
   \item \code{transf.dtou1}: transformation of standardized mean differences to Cohen's U1 values (Cohen, 1988).
   \item \code{transf.dtou2}: transformation of standardized mean differences to Cohen's U2 values (Cohen, 1988).
   \item \code{transf.dtou3}: transformation of standardized mean differences to Cohen's U3 values (Cohen, 1988).
   \item \code{transf.dtocles}: transformation of standardized mean differences to common language effect size values (McGraw & Wong, 1992).
   \item \code{transf.dtorpb}: transformation of standardized mean differences to point-biserial correlations. If \code{n1i} and \code{n2i} are not specified, the function assumes \code{n1i=n2i} and uses an approximate formula. If \code{n1i} and \code{n2i} are specified, the exact transformation formula is used.
   \item \code{transf.dtobesd}: transformation of standardized mean differences to binomial effect size display values (Rosenthal & Rubin, 1982). Note that the function only provides the proportion in the first group scoring above the median (the proportion in the second group scoring above the median is simply one minus the proportion in the first group scoring above the median).
   \item \code{transf.dtoumd}: transformation of standardized mean differences to mean differences given a known standard deviation, which must be specified via the \code{targs} argument.
   \item \code{transf.logortord}: transformation of log odds ratios to risk differences, assuming a particular value for the control group risk (which needs to be specified via the \code{pc} argument).
   \item \code{transf.logortorr}: transformation of log odds ratios to risk ratios, assuming a particular value for the control group risk (which needs to be specified via the \code{pc} argument).
   }
}
\value{
   A vector with the transformed values.
}
\note{
   The integral transformation method for a transformation function \mjseqn{h(z)} integrates \mjseqn{h(z) f(z)} over \mjseqn{z} using the limits \code{targs$lower} and \code{targs$upper}, where \mjseqn{f(z)} is the density of a normal distribution with mean equal to \code{xi} and variance equal to \code{targs$tau2}. An example is provided below.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient alpha. \emph{Journal of Educational and Behavioral Statistics}, \bold{27}(4), 335--340. \verb{https://doi.org/10.3102/10769986027004335}

   Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.

   Fisher, R. A. (1921). On the \dQuote{probable error} of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32. \verb{http://hdl.handle.net/2440/15169}

   Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular and the square root. \emph{Annals of Mathematical Statistics}, \bold{21}(4), 607--611. \verb{https://doi.org/10.1214/aoms/1177729756}

   Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent alpha coefficients. \emph{Psychometrika}, \bold{41}(2), 219--231. \verb{https://doi.org/10.1007/BF02291840}

   McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. \emph{Psychological Bulletin}, \bold{111}(2), 361--365. \verb{https://doi.org/10.1037/0033-2909.111.2.361}

   Miller, J. J. (1978). The inverse of the Freeman-Tukey double arcsine transformation. \emph{American Statistician}, \bold{32}(4), 138. \verb{https://doi.org/10.1080/00031305.1978.10479283}

   Rosenthal, R., & Rubin, D. B. (1982). A simple, general purpose display of magnitude of experimental effect. \emph{Journal of Educational Psychology}, \bold{74}(2), 166--169. \verb{https://doi.org/10.1037/0022-0663.74.2.166}

   Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. \emph{Journal of Statistical Software}, \bold{36}(3), 1--48. \verb{https://doi.org/10.18637/jss.v036.i03}
}
\examples{
### calculate log risk ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects model
res <- rma(yi, vi, data=dat)

### average risk ratio with 95\% CI (but technically, this provides an
### estimate of the median risk ratio, not the mean risk ratio!)
predict(res, transf=exp)

### average risk ratio with 95\% CI using the integral transformation
predict(res, transf=transf.exp.int, targs=list(tau2=res$tau2, lower=-4, upper=4))
}
\keyword{manip}
