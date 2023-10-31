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
\alias{transf.r2toz}
\alias{transf.ztor2}
\alias{transf.ztor.int}
\alias{transf.exp.int}
\alias{transf.ilogit.int}
\alias{transf.dtou1}
\alias{transf.dtou2}
\alias{transf.dtou3}
\alias{transf.dtocles}
\alias{transf.dtobesd}
\alias{transf.dtomd}
\alias{transf.dtorpb}
\alias{transf.dtorbis}
\alias{transf.rpbtorbis}
\alias{transf.rtorpb}
\alias{transf.rtod}
\alias{transf.rpbtod}
\alias{transf.lnortord}
\alias{transf.lnortorr}
\alias{transf.lnortod.norm}
\alias{transf.lnortod.logis}
\alias{transf.dtolnor.norm}
\alias{transf.dtolnor.logis}
\alias{transf.lnortortet.pearson}
\alias{transf.lnortortet.digby}
\title{Transformation Functions}
\description{
   Functions to carry out various types of transformations that are useful for meta-analyses. \loadmathjax
}
\usage{
transf.rtoz(xi)
transf.ztor(xi)
transf.logit(xi)
transf.ilogit(xi)
transf.arcsin(xi)
transf.iarcsin(xi)
transf.pft(xi, ni)
transf.ipft(xi, ni)
transf.ipft.hm(xi, targs)
transf.isqrt(xi)
transf.irft(xi, ti)
transf.iirft(xi, ti)
transf.ahw(xi)
transf.iahw(xi)
transf.abt(xi)
transf.iabt(xi)
transf.r2toz(xi)
transf.ztor2(xi)
transf.ztor.int(xi, targs)
transf.exp.int(xi, targs)
transf.ilogit.int(xi, targs)
transf.dtou1(xi)
transf.dtou2(xi)
transf.dtou3(xi)
transf.dtobesd(xi)
transf.dtomd(xi, targs)
transf.dtorpb(xi, n1i, n2i)
transf.dtorbis(xi, n1i, n2i)
transf.rpbtorbis(xi, pi)
transf.rtorpb(xi, pi)
transf.rtod(xi, n1i, n2i)
transf.rpbtod(xi, n1i, n2i)
transf.lnortord(xi, pc)
transf.lnortorr(xi, pc)
transf.lnortod.norm(xi)
transf.lnortod.logis(xi)
transf.dtolnor.norm(xi)
transf.dtolnor.logis(xi)
transf.lnortortet.pearson(xi)
transf.lnortortet.digby(xi)
}
\arguments{
   \item{xi}{vector of values to be transformed.}
   \item{ni}{vector of sample sizes.}
   \item{n1i}{vector of sample sizes for the first group.}
   \item{n2i}{vector of sample sizes for the second group.}
   \item{ti}{vector of person-times at risk.}
   \item{pc}{control group risk (either a single value or a vector).}
   \item{pi}{proportion of individuals falling into the first of the two groups that is created by the dichotomization.}
   \item{targs}{list with additional arguments for the transformation function. See \sQuote{Details}.}
}
\details{
   The following transformation functions are currently implemented:

   \itemize{
   \item \code{transf.rtoz}: Fisher's r-to-z transformation for correlation coefficients (same as \code{atanh(x)}).
   \item \code{transf.ztor}: inverse of the former (i.e., the z-to-r transformation; same as \code{tanh(x)}).
   \item \code{transf.logit}: logit (log odds) transformation for proportions (same as \code{qlogis(x)}).
   \item \code{transf.ilogit}: inverse of the former (same as \code{plogis(x)}).
   \item \code{transf.arcsin}: arcsine square root transformation for proportions.
   \item \code{transf.iarcsin}: inverse of the former.
   \item \code{transf.pft}: Freeman-Tukey (double arcsine) transformation for proportions. See Freeman & Tukey (1950). The \code{xi} argument is used to specify the proportions and the \code{ni} argument the corresponding sample sizes.
   \item \code{transf.ipft}: inverse of the former. See Miller (1978).
   \item \code{transf.ipft.hm}: inverse of the former, using the harmonic mean of the sample sizes for the back-transformation. See Miller (1978). The sample sizes are specified via the \code{targs} argument (the list element should be called \code{ni}).
   \item \code{transf.isqrt}: inverse of the square root transformation (i.e., function to square a number).
   \item \code{transf.irft}: Freeman-Tukey transformation for incidence rates. See Freeman & Tukey (1950). The \code{xi} argument is used to specify the incidence rates and the \code{ti} argument the corresponding person-times at risk.
   \item \code{transf.iirft}: inverse of the former.
   \item \code{transf.ahw}: transformation of coefficient alpha as suggested by Hakstian & Whalen (1976), except that \mjeqn{1-(1-\alpha)^{1/3}}{1-(1-\alpha)^(1/3)} is used (so that the transformed values are a monotonically increasing function of the \mjseqn{\alpha} values).
   \item \code{transf.iahw}: inverse of the former.
   \item \code{transf.abt}: transformation of coefficient alpha as suggested by Bonett (2002), except that \mjeqn{-\ln(1-\alpha)}{-ln(1-\alpha)} is used (so that the transformed values are a monotonically increasing function of the \mjseqn{\alpha} values).
   \item \code{transf.iabt}: inverse of the former.
   \item \code{transf.r2toz}: variance stabilizing transformation for the coefficient of determination, given by \mjeqn{z_i = \frac{1}{2} \ln\mathopen{}\left(\frac{1+\sqrt{R_i^2}}{1-\sqrt{R_i^2}}\right)\mathclose{}}{z_i = 1/2 log((1+\sqrt(R_i^2))/(1-\sqrt(R_i^2)))} (see Olkin & Finn, 1995, but with the additional \mjeqn{\frac{1}{2}}{1/2} factor).
   \item \code{transf.ztor2}: inverse of the former.
   \item \code{transf.ztor.int}: integral transformation method for the z-to-r transformation.
   \item \code{transf.exp.int}: integral transformation method for the exponential transformation.
   \item \code{transf.ilogit.int}: integral transformation method for the inverse logit transformation.
   \item \code{transf.dtou1}: transformation of standardized mean differences to Cohen's \mjseqn{U_1} values (Cohen, 1988). Under the assumption that the data for those in the first (say, treated) and second (say, control) group are normally distributed with equal variances but potentially different means, Cohen's \mjseqn{U_1} indicates the proportion of non-overlap between the two distributions (i.e., when \mjseqn{d=0}, then \mjseqn{U_1} is equal to 0, which goes to 1 as \mjseqn{d} increases).
   \item \code{transf.dtou2}: transformation of standardized mean differences to Cohen's \mjseqn{U_2} values (Cohen, 1988). Under the same assumptions as above, Cohen's \mjseqn{U_2} indicates the proportion in the first group that exceeds the same proportion in the second group (i.e., when \mjseqn{d=0}, then \mjseqn{U_2} is equal to 0.5, which goes to 1 as \mjseqn{d} increases).
   \item \code{transf.dtou3}: transformation of standardized mean differences to Cohen's \mjseqn{U_3} values (Cohen, 1988). Under the same assumptions as above, Cohen's \mjseqn{U_3} indicates the proportion of individuals in the first group that have a higher value than the mean of those in the second group (i.e., when \mjseqn{d=0}, then \mjseqn{U_3} is equal to 0.5, which goes to 1 as \mjseqn{d} increases).
   \item \code{transf.dtocles}: transformation of standardized mean differences to common language effect size (CLES) values (McGraw & Wong, 1992) (also called the probability of superiority). A CLES value indicates the probability that a randomly sampled individual from the first group has a higher value than a randomly sampled individual from the second group (i.e., when \mjseqn{d=0}, then the CLES is equal to 0.5, which goes to 1 as \mjseqn{d} increases).
   \item \code{transf.dtobesd}: transformation of standardized mean differences to binomial effect size display values (Rosenthal & Rubin, 1982). Note that the function only provides the proportion in the first group scoring above the median (the proportion in the second group scoring above the median is simply one minus this value).
   \item \code{transf.dtomd}: transformation of standardized mean differences to mean differences given a known standard deviation (which needs to be specified via the \code{targs} argument).
   \item \code{transf.dtorpb}: transformation of standardized mean differences to point-biserial correlations. Arguments \code{n1i} and \code{n2i} denote the number of individuals in the first and second group, respectively. If \code{n1i} and \code{n2i} are not specified, the function assumes \code{n1i = n2i} and uses the approximate formula \mjeqn{r_{pb} = \frac{d}{\sqrt{d^2 + 4}}}{r_pb = d / \sqrt{d^2 + 4}}. If \code{n1i} and \code{n2i} are specified, the function uses the exact transformation formula \mjeqn{r_{pb} = \frac{d}{\sqrt{d^2 + h}}}{r_pb = d / \sqrt{d^2 + h}}, where \mjeqn{h = \frac{m}{n_1} + \frac{m}{n_2}}{h = m / n_1 + m / n_2} and \mjseqn{m = n_1 + n_2 - 2} (Jacobs & Viechtbauer, 2017).
   \item \code{transf.dtorbis}: transformation of standardized mean differences to biserial correlations. Like \code{transf.dtorpb}, but the point-biserial correlations are then transformed to biserial correlations with \mjeqn{r_{bis} = \frac{\sqrt{p(1-p)}}{f(z_p)} r_{pb}}{r_bis = sqrt(p*(1-p)) / f(z_p) r_pb}, where \mjseqn{p = \frac{n1}{n1+n2}} and \mjseqn{f(z_p)} denotes the density of the standard normal distribution at value \mjseqn{z_p}, which is the point for which \mjseqn{P(Z > z_p) = p}, with \mjseqn{Z} denoting a random variable following a standard normal distribution (Jacobs & Viechtbauer, 2017).
   \item \code{transf.rpbtorbis}: transformation of point-biserial correlations to biserial correlations. Argument \code{pi} denotes the proportion of individuals falling into the first of the two groups that is created by the dichotomization (hence, \code{1-pi} falls into the second group). If \code{pi} is not specified, the function assumes \code{pi=0.5}, which corresponds to dichotomization at the median. The transformation is carried out as described for \code{transf.dtorbis}.
   \item \code{transf.rtorpb}: transformation of Pearson product-moment correlations to the corresponding point-biserial correlations, when one of the two variables is dichotomized. Argument \code{pi} can be used to denote the proportion of individuals falling into the first of the two groups that is created by the dichotomization (hence, \code{1-pi} falls into the second group). If \code{pi} is not specified, the function assumes \code{pi=0.5}, which corresponds to dichotomization at the median. This function is simply the inverse of \code{transf.rpbtorbis}.
   \item \code{transf.rtod}: transformation of Pearson product-moment correlations to the corresponding standardized mean differences, when one of the two variables is dichotomized. Arguments \code{n1i} and \code{n2i} can be used to denote the number of individuals in the first and second group created by the dichotomization. If \code{n1i} and \code{n2i} are not specified, the function assumes \code{n1i = n2i}. This function is simply the inverse of \code{transf.dtorbis}.
   \item \code{transf.rpbtod}: transformation of point-biserial correlations to standardized mean differences. This is simply the inverse of \code{transf.dtorpb}.
   \item \code{transf.lnortord}: transformation of log odds ratios to risk differences, assuming a particular value for the control group risk (which needs to be specified via the \code{pc} argument).
   \item \code{transf.lnortorr}: transformation of log odds ratios to risk ratios, assuming a particular value for the control group risk (which needs to be specified via the \code{pc} argument).
   \item \code{transf.lnortod.norm}: transformation of log odds ratios to standardized mean differences (assuming normal distributions) (Cox & Snell, 1989).
   \item \code{transf.lnortod.logis}: transformation of log odds ratios to standardized mean differences (assuming logistic distributions) (Chinn, 2000).
   \item \code{transf.dtolnor.norm}: transformation of standardized mean differences to log odds ratios (assuming normal distributions) (Cox & Snell, 1989).
   \item \code{transf.dtolnor.logis}: transformation of standardized mean differences to log odds ratios (assuming logistic distributions) (Chinn, 2000).
   \item \code{transf.lnortortet.pearson}: transformation of log odds ratios to tetrachoric correlations as suggested by Pearson (1900).
   \item \code{transf.lnortortet.digby}: transformation of log odds ratios to tetrachoric correlations as suggested by Digby (1983).
   }
}
\value{
   A vector with the transformed values.
}
\note{
   The integral transformation method for a transformation function \mjseqn{h(z)} is given by \mjsdeqn{\int_{\textrm{lower}}^{\textrm{upper}} h(z) f(z) dz} using the limits \code{targs$lower} and \code{targs$upper}, where \mjseqn{f(z)} is the density of a normal distribution with mean equal to \code{xi} and variance equal to \code{targs$tau2}. An example is provided below.
}
\author{
   Wolfgang Viechtbauer \email{wvb@metafor-project.org} \url{https://www.metafor-project.org}
}
\references{
   Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient alpha. \emph{Journal of Educational and Behavioral Statistics}, \bold{27}(4), 335--340. \verb{https://doi.org/10.3102/10769986027004335}

   Chinn, S. (2000). A simple method for converting an odds ratio to effect size for use in meta-analysis. \emph{Statistics in Medicine}, \bold{19}(22), 3127--3131. \verb{https://doi.org/10.1002/1097-0258(20001130)19:22<3127::aid-sim784>3.0.co;2-m}

   Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.

   Cox, D. R., & Snell, E. J. (1989). \emph{Analysis of binary data} (2nd ed.). London: Chapman & Hall.

   Digby, P. G. N. (1983). Approximating the tetrachoric correlation coefficient. \emph{Biometrics}, \bold{39}(3), 753--757. \verb{https://doi.org/10.2307/2531104}

   Fisher, R. A. (1921). On the \dQuote{probable error} of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32. \verb{http://hdl.handle.net/2440/15169}

   Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular and the square root. \emph{Annals of Mathematical Statistics}, \bold{21}(4), 607--611. \verb{https://doi.org/10.1214/aoms/1177729756}

   Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent alpha coefficients. \emph{Psychometrika}, \bold{41}(2), 219--231. \verb{https://doi.org/10.1007/BF02291840}

   Jacobs, P., & Viechtbauer, W. (2017). Estimation of the biserial correlation and its sampling variance for use in meta-analysis. \emph{Research Synthesis Methods}, \bold{8}(2), 161--180. \verb{https://doi.org/10.1002/jrsm.1218}

   McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. \emph{Psychological Bulletin}, \bold{111}(2), 361--365. \verb{https://doi.org/10.1037/0033-2909.111.2.361}

   Miller, J. J. (1978). The inverse of the Freeman-Tukey double arcsine transformation. \emph{American Statistician}, \bold{32}(4), 138. \verb{https://doi.org/10.1080/00031305.1978.10479283}

   Olkin, I., & Finn, J. D. (1995). Correlations redux. \emph{Psychological Bulletin}, \bold{118}(1), 155--164. \verb{https://doi.org/10.1037/0033-2909.118.1.155}

   Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. \emph{Philosophical Transactions of the Royal Society of London, Series A}, \bold{195}, 1--47. \verb{https://doi.org/10.1098/rsta.1900.0022}

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

### this also works
predict(res, transf=transf.exp.int, targs=list(tau2=res$tau2))

### this as well
predict(res, transf=transf.exp.int, targs=res$tau2)
}
\keyword{manip}
