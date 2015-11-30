\name{dat.hackshaw1998}
\docType{data}
\alias{dat.hackshaw1998}
\title{Studies on Lung Cancer Risk from ETS Exposure}
\description{Results from 37 studies on the risk of lung cancer from environmental tobacco smoke (ETS) exposure.}
\usage{dat.hackshaw1998}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{numeric}    \tab study number \cr
\bold{author}     \tab \code{character}  \tab first author of study \cr
\bold{year}       \tab \code{numeric}    \tab publication year \cr
\bold{country}    \tab \code{character}  \tab country where study was conducted \cr
\bold{design}     \tab \code{character}  \tab either cohort or case-control study \cr
\bold{cases}      \tab \code{numeric}    \tab number of lung cancer cases \cr
\bold{yi}         \tab \code{numeric}    \tab log odds ratio \cr
\bold{vi}         \tab \code{numeric}    \tab corresponding sampling variance
}
}
\details{
   The dataset includes the results from 37 studies (4 cohort, 33 case-control) on the risk of lung cancer from environmental tobacco smoke (ETS) exposure from the spouse in women who are lifelong nonsmokers. Values of the log odds ratio greater than 0 indicate an increased risk of cancer in exposed women compared to women not exposed to ETS from their spouse.
}
\source{
   Hackshaw, A. K., Law, M. R., & Wald, N. J. (1997). The accumulated evidence on lung cancer and environmental tobacco smoke. \emph{British Medical Journal}, \bold{315}, 980--988.

   Hackshaw, A. K. (1998). Lung cancer and passive smoking. \emph{Statistical Methods in Medical Research}, \bold{7}, 119--136.
}
\examples{
### load data
dat <- get(data(dat.hackshaw1998))

### random-effects model using the log odds ratios
res <- rma(yi, vi, data=dat, method="DL")
res

### estimated average odds ratio with CI (and credibility/prediction interval)
predict(res, transf=exp, digits=2)
}
\keyword{datasets}
