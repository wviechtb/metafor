\name{dat.pagliaro1992}
\docType{data}
\alias{dat.pagliaro1992}
\title{Studies on the Effectiveness of Nonsurgical Treatments in Cirrhosis}
\description{Results from 26 trials examining the effectiveness of beta-blockers and sclerotherapy for the prevention of first bleeding in patients with cirrhosis}
\usage{dat.pagliaro1992}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{study}    \tab \code{numeric}   \tab study id \cr
\bold{trt}      \tab \code{character} \tab either beta-blockers, sclerotherapy, or control \cr
\bold{xi}       \tab \code{numeric}   \tab number of patients with first bleeding \cr
\bold{ni}       \tab \code{numeric}   \tab number of patients treated
}
}
\details{
   The dataset includes the results from 26 randomized controlled trials examining the effectiveness of nonsurgical treatments for the prevention of first bleeding in patients with cirrhosis. Patients were either treated with beta-blockers, endoscopic sclerotherapy, or with a nonactive treatment (control). Two trials included all three treatment conditions, 7 trials compared beta-blockers against control, and 17 trials compared sclerotherapy against control. The dataset has been used in various papers to illustrate methods for conducting a network meta-analysis / mixed treatment comparison.
}
\source{
   Pagliaro, L., D'Amico, G., \enc{Sörensen}{Soerensen}, T. I. A., Lebrec, D., Burroughs, A. K., Morabito, A., \enc{Tiné}{Tine}, F., Politi, F., & Traina, M. (1992). Prevention of first bleeding in cirrhosis: A meta-analysis of randomized trials of nonsurgical treatment. \emph{Annals of Internal Medicine}, \bold{117}, 59--70.
}
\examples{
### copy data into 'dat' and examine data
dat <- dat.pagliaro1992
dat

### restructure dataset to a contrast-based format
dat.c <- to.wide(dat, study="study", grp="trt", grpvars=3:4)
dat.c

### Mantel-Haenszel results for beta-blockers and sclerotherapy versus control, respectively
rma.mh(measure="OR", ai=xi.1, n1i=ni.1, ci=xi.2, n2i=ni.2,
       data=dat.c, subset=(trt.1=="beta-blockers"), digits=2)
rma.mh(measure="OR", ai=xi.1, n1i=ni.1, ci=xi.2, n2i=ni.2,
       data=dat.c, subset=(trt.1=="sclerotherapy"), digits=2)

### calculate log odds for each study arm
dat <- escalc(measure="PLO", xi=xi, ni=ni, data=dat)
dat

### turn treatment variable into factor and set reference level
dat$trt <- relevel(factor(dat$trt), ref="control")

### add a space before each level (this makes the output a bit more legible)
levels(dat$trt) <- paste0(" ", levels(dat$trt))

### network meta-analysis using an arm-based random-effects model with fixed study effects
### (by setting rho=1/2, tau^2 reflects the amount of heterogeneity for all treatment comparisons)
res <- rma.mv(yi, vi, mods = ~ factor(study) + trt - 1, random = ~ trt | study, rho=1/2, data=dat)
res

### average odds ratio comparing beta-blockers and sclerotherapy versus control, respectively
predict(res, newmods=c(rep(0,26), 1, 0), transf=exp, digits=2)
predict(res, newmods=c(rep(0,26), 0, 1), transf=exp, digits=2)

### average odds ratio comparing beta-blockers versus sclerotherapy
predict(res, newmods=c(rep(0,26), 1, -1), transf=exp, digits=2)
}
\keyword{datasets}
