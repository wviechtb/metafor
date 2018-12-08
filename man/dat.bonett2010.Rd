\name{dat.bonett2010}
\docType{data}
\alias{dat.bonett2010}
\title{Studies on the Reliability of the CES-D Scale}
\description{Results from 9 studies on the reliability of the Center for Epidemiologic Studies Depression (CES-D) Scale administered to children providing care to an elderly parent.}
\usage{dat.bonett2010}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}      \tab \code{numeric}    \tab study number \cr
\bold{source}     \tab \code{character}  \tab source of data \cr
\bold{ni}         \tab \code{numeric}    \tab sample size \cr
\bold{mi}         \tab \code{numeric}    \tab number of items in the scale \cr
\bold{ai}         \tab \code{numeric}    \tab observed value of Cronbach's alpha \cr
\bold{caregivers} \tab \code{character}  \tab gender of the children in the sample
}
}
\details{
   The Center for Epidemiologic Studies Depression (CES-D) Scale is a 20-item questionnaire assessing various symptoms of depression, with each item scored on a 4-point scale. The scale has been used in several studies to examine depressive symptoms in children providing care to an elderly parent. The dataset includes information on the reliability of the scale as measured with Cronbach's alpha in 9 such studies. Also, the gender composition of the children in each sample is indicated.
}
\source{
   Bonett, D. G. (2010). Varying coefficient meta-analytic methods for alpha reliability. \emph{Psychological Methods}, \bold{15}, 368--385.
}
\references{
   Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient alpha. \emph{Journal of Educational and Behavioral Statistics}, \bold{27}, 335--340.

   Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent alpha coefficients. \emph{Psychometrika}, \bold{41}, 219--231.
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.bonett2010
dat

### meta-analysis using the raw alpha values
res <- rma(measure="ARAW", ai=ai, mi=mi, ni=ni, data=dat)
res

### meta-analysis using transformed alpha values (using the
### transformation suggested by Hakstian & Whalen, 1976)
res <- rma(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat)
res
predict(res, transf=transf.iahw)

### meta-analysis using transformed alpha values (using the
### transformation suggested by Bonett, 2002)
res <- rma(measure="ABT", ai=ai, mi=mi, ni=ni, data=dat)
res
predict(res, transf=transf.iabt)

### examine whether female/mixed samples yield different alphas (with raw alphas)
res <- rma(measure="ARAW", ai=ai, mi=mi, ni=ni, mods = ~ caregivers, data=dat)
res
predict(res, newmods=c(0,1), digits=2)
}
\keyword{datasets}
