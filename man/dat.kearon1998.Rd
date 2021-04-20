\name{dat.kearon1998}
\docType{data}
\alias{dat.kearon1998}
\title{Studies on the Accuracy of Venous Ultrasonography for the Diagonosis of Deep Venous Thrombosis}
\description{Results from diagnostic accuracy studies examining the accuracy of venous ultrasonography for the diagonosis of deep venous thrombosis.}
\usage{
dat.kearon1998
}
\format{The data frame contains the following columns:
\tabular{lll}{
\bold{id}       \tab \code{numeric}   \tab study id \cr
\bold{author}   \tab \code{character} \tab study author(s) \cr
\bold{year}     \tab \code{numeric}   \tab publication year \cr
\bold{patients} \tab \code{character} \tab patient group (either symptomatic or asymptomatic patients) \cr
\bold{tp}       \tab \code{numeric}   \tab number of true positives \cr
\bold{np}       \tab \code{numeric}   \tab number of positive patients (cases) \cr
\bold{tn}       \tab \code{numeric}   \tab number of true negatives \cr
\bold{nn}       \tab \code{numeric}   \tab number of negative patients (non-cases)
}
}
\details{
   The studies included in the dataset examined the accuracy of venous ultrasonography for the diagnossis of a first deep venous thrombosis in symptomatic and asymptomatic patients. Cases and non-cases were determined based on contrast venography. Venous ultrasonography was then used to make a diagnosis, leading to a given number of true positives and negatives.

   A subset of this dataset (using only the studies with asymptomatic patients) was used by Deeks et al. (2005) to illustrate methods for detecting publication bias (or small-study effects) in meta-analyses of diagnostic accuracy studies.
}
\source{
   Kearon, C., Julian, J. A., Math, M., Newman, T. E., & Ginsberg, J. S. (1998). Noninvasive diagnosis of deep venous thrombosis. \emph{Annals of Internal Medicine}, \bold{128}(8), 663--677. \verb{https://doi.org/10.7326/0003-4819-128-8-199804150-00011}
}
\references{
   Deeks, J. J., Macaskill, P., & Irwig, L. (2005). The performance of tests of publication bias and other sample size effects in systematic reviews of diagnostic test accuracy was assessed. \emph{Journal of Clinical Epidemiology}, \bold{58}(9), 882--893. \verb{https://doi.org/10.1016/j.jclinepi.2005.01.016}
}
\examples{
### copy data into 'dat'
dat <- dat.kearon1998

### calculate diagnostic log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=tp, n1i=np, ci=nn-tn, n2i=nn, data=dat, add=1/2, to="all")
dat

### fit random-effects model for the symptomatic patients
res <- rma(yi, vi, data=dat, subset=patients=="symptomatic")
res

### fit random-effects model for the asymptomatic patients
res <- rma(yi, vi, data=dat, subset=patients=="asymptomatic")
res

### estimated average diagnostic odds ratio (with 95% CI)
predict(res, transf=exp, digits=2)

### regression test for funnel plot asymmetry using SE as predictor
reg <- regtest(res, model="lm")
reg

### corresponding funnel plot
funnel(res, atransf=exp, xlim=c(0,7), at=log(c(1,10,100,1000)), ylim=c(0,1.5), steps=4)
ys <- seq(0, 2, length=100)
lines(coef(reg$fit)[1] + coef(reg$fit)[2]*ys, ys, lwd=2, lty=3)

### regression test for funnel plot asymmetry using total sample size as predictor
reg <- regtest(res, model="lm", predictor="ni")
reg

### corresponding funnel plot
funnel(res, yaxis="ni", atransf=exp, xlim=c(0,7), at=log(c(1,10,100,1000)), ylim=c(0,300), steps=4)
ys <- seq(0, 300, length=100)
lines(coef(reg$fit)[1] + coef(reg$fit)[2]*ys, ys, lwd=2, lty=3)

### regression test for funnel plot asymmetry using 1/sqrt(ESS) as predictor (Deeks et al., 2005)
dat$invessi <- 1/(4*dat$np) + 1/(4*dat$nn)
tmp <- rma(yi, invessi, data=dat, subset=patients=="asymptomatic")
reg <- regtest(tmp, model="lm")
reg

### corresponding funnel plot
funnel(tmp, atransf=exp, xlim=c(0,7), at=log(c(1,10,100,1000)), ylim=c(0,.15), steps=4,
       refline=coef(res), level=0, ylab="1/root(ess)")
ys <- seq(0, .20, length=100)
lines(coef(reg$fit)[1] + coef(reg$fit)[2]*ys, ys, lwd=2, lty=3)

### convert data to long format
dat <- to.long(measure="OR", ai=tp, n1i=np, ci=tn, n2i=nn,
               data=dat.kearon1998, subset=patients=="asymptomatic")
dat <- dat[9:12]
levels(dat$group) <- c("sensitivity", "specificity")
dat

### calculate logit-transformed sensitivities
dat <- escalc(measure="PLO", xi=out1, mi=out2, data=dat, add=1/2, to="all",
              include=group=="sensitivity")
dat

### calculate logit-transformed specificities
dat <- escalc(measure="PLO", xi=out1, mi=out2, data=dat, add=1/2, to="all",
              include=group=="specificity")
dat

### bivariate random-effects model for logit sensitivity and specificity
res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | study, struct="UN", data=dat)
res

### estimated average sensitivity and specificity based on the model
predict(res, newmods = rbind(c(1,0),c(0,1)), transf=transf.ilogit, tau2.levels=c(1,2), digits=2)

### estimated average diagnostic odds ratio based on the model
predict(res, newmods = c(1,1), transf=exp, digits=2)
}
\keyword{datasets}
