\name{dat.collins1985b}
\docType{data}
\alias{dat.collins1985b}
\title{Studies on the Effects of Diuretics in Pregnancy}
\description{Results from 9 studies examining the effects of diuretics in pregnancy on various outcomes.}
\usage{dat.collins1985b}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}      \tab \code{numeric}   \tab study number \cr
\bold{author}  \tab \code{character} \tab study author(s) \cr
\bold{year}    \tab \code{numeric}   \tab publication year \cr
\bold{pre.nti} \tab \code{numeric}   \tab number of women in treatment group followed up for pre-eclampsia outcome \cr
\bold{pre.nci} \tab \code{numeric}   \tab number of women in control/placebo group followed up for pre-eclampsia outcome \cr
\bold{pre.xti} \tab \code{numeric}   \tab number of women in treatment group with any form of pre-eclampsia \cr
\bold{pre.xci} \tab \code{numeric}   \tab number of women in control/placebo group with any form of pre-eclampsia \cr
\bold{oedema}  \tab \code{numeric}   \tab dummy variable indicating whether oedema was a diagnostic criterion \cr
\bold{fup.nti} \tab \code{numeric}   \tab number of women in treatment group followed up for mortality outcomes \cr
\bold{fup.nci} \tab \code{numeric}   \tab number of women in control/placebo group followed up for mortality outcomes \cr
\bold{ped.xti} \tab \code{numeric}   \tab number of perinatal deaths in treatment group \cr
\bold{ped.xci} \tab \code{numeric}   \tab number of perinatal deaths in control/placebo group \cr
\bold{stb.xti} \tab \code{numeric}   \tab number of stillbirths in treatment group \cr
\bold{stb.xci} \tab \code{numeric}   \tab number of stillbirths in control/placebo group \cr
\bold{ned.xti} \tab \code{numeric}   \tab number of neonatal deaths in treatment group \cr
\bold{ned.xci} \tab \code{numeric}   \tab number of neonatal deaths in control/placebo group
}
}
\details{
   The 9 studies in this dataset examined the effects of diuretics in pregnancy on various outcomes, including the presence of any form of pre-eclampsia, perinatal death, stillbirth, and neonatal death.
}
\source{
   Collins, R., Yusuf, S., & Peto, R. (1985). Overview of randomised trials of diuretics in pregnancy. \emph{British Medical Journal}, \bold{290}(6461), 17--23. \verb{https://doi.org/10.1136/bmj.290.6461.17}
}
\examples{
### copy data into 'dat'
dat <- dat.collins1985b

### calculate (log) odds ratio and sampling variance
dat <- escalc(measure="OR", n1i=pre.nti, n2i=pre.nci, ai=pre.xti, ci=pre.xci, data=dat)
summary(dat, digits=2, transf=exp)

### meta-analysis using Peto's method for any form of pre-eclampsia
rma.peto(n1i=pre.nti, n2i=pre.nci, ai=pre.xti, ci=pre.xci, data=dat, digits=2)

### meta-analysis including only studies where oedema was not a diagnostic criterion
rma.peto(n1i=pre.nti, n2i=pre.nci, ai=pre.xti, ci=pre.xci, data=dat, digits=2, subset=(oedema==0))

### meta-analyses of mortality outcomes (perinatal deaths, stillbirths, and neonatal deaths)
rma.peto(n1i=fup.nti, n2i=fup.nci, ai=ped.xti, ci=ped.xci, data=dat, digits=2)
rma.peto(n1i=fup.nti, n2i=fup.nci, ai=stb.xti, ci=stb.xci, data=dat, digits=2)
rma.peto(n1i=fup.nti, n2i=fup.nci, ai=ned.xti, ci=ned.xci, data=dat, digits=2)
}
\keyword{datasets}
