\name{dat.hackshaw1998}
\docType{data}
\alias{dat.hackshaw1998}
\title{Studies on the Risk of Lung Cancer in Women Exposed to Environmental Tobacco Smoke}
\description{Results from 37 studies on the risk of lung cancer in women exposed to environmental tobacco smoke (ETS) from their smoking spouse.}
\usage{dat.hackshaw1998}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{numeric}    \tab study number \cr
\bold{author}     \tab \code{character}  \tab first author of study \cr
\bold{year}       \tab \code{numeric}    \tab publication year \cr
\bold{country}    \tab \code{character}  \tab country where study was conducted \cr
\bold{design}     \tab \code{character}  \tab study design (either cohort or case-control) \cr
\bold{cases}      \tab \code{numeric}    \tab number of lung cancer cases \cr
\bold{or}         \tab \code{numeric}    \tab odds ratio \cr
\bold{or.lb}      \tab \code{numeric}    \tab lower bound of 95\% CI for the odds ratio \cr
\bold{or.ub}      \tab \code{numeric}    \tab upper bound of 95\% CI for the odds ratio \cr
\bold{yi}         \tab \code{numeric}    \tab log odds ratio \cr
\bold{vi}         \tab \code{numeric}    \tab corresponding sampling variance
}
}
\details{
   The dataset includes the results from 37 studies (4 cohort, 33 case-control) examining if women (who are lifelong nonsmokers) have an elevated risk for lung cancer due to exposure to environmental tobacco smoke (ETS) from their smoking spouse. Values of the log odds ratio greater than 0 indicate an increased risk of cancer in exposed women compared to women not exposed to ETS from their spouse.

   Note that the log odds ratios and corresponding sampling variances were back-calculated from the reported odds ratios and confidence interval (CI) bounds (see \sQuote{Examples}). Since the reported values were rounded to some extent, this introduces some minor inaccuracies into the back-calculations. The overall estimate reported in Hackshaw et al. (1997) and Hackshaw (1998) can be fully reproduced though.
}
\source{
   Hackshaw, A. K., Law, M. R., & Wald, N. J. (1997). The accumulated evidence on lung cancer and environmental tobacco smoke. \emph{British Medical Journal}, \bold{315}(7114), 980--988. \verb{https://doi.org/10.1136/bmj.315.7114.980}

   Hackshaw, A. K. (1998). Lung cancer and passive smoking. \emph{Statistical Methods in Medical Research}, \bold{7}(2), 119--136. \verb{https://doi.org/10.1177/096228029800700203}
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.hackshaw1998
dat

### random-effects model using the log odds ratios
res <- rma(yi, vi, data=dat, method="DL")
res

### estimated average odds ratio with CI (and prediction interval)
predict(res, transf=exp, digits=2)

### illustrate how the log odds ratios and corresponding sampling variances
### were back-calculated based on the reported odds ratios and CI bounds
dat$yi <- log(dat$or)
dat$vi <- ((log(dat$or.ub) - log(dat$or.lb)) / (2*qnorm(.975)))^2
}
\keyword{datasets}
