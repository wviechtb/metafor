### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: .setlab() function")

yi <- c(-.3, -.1, 0, .2, .2)
vi <- rep(.02, length(yi))

test_that(".setlab() works correctly together with forest().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   opar <- par(no.readonly=TRUE)

   par(mfrow=c(5,3), mar=c(5,6,0,4))
   xlim <- c(-3,5)
   cex.lab <- .5

   dat <- escalc(measure="GEN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="RR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="OR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="RD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="AS", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="PHI", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="YUQ", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="YUY", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="IRR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="IRD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="IRSD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="MD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="SMD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="ROM", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="CVR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="VR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="RPB", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="COR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="ZCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ztor)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ztor)

   dat <- escalc(measure="PR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="PLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="PLO", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ilogit)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ilogit)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="PAS", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.iarcsin)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.iarcsin)

   dat <- escalc(measure="PFT", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ipft.hm, targs=list(ni=rep(10,length(yi))))
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ipft.hm, targs=list(ni=rep(10,length(yi))))

   dat <- escalc(measure="IR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="IRLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="IRS", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.isqrt)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.isqrt)

   dat <- escalc(measure="IRFT", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="MN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="MNLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="CVLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="SDLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="MC", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="SMCC", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="ROMC", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp)

   dat <- escalc(measure="ARAW", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="AHW", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.iahw)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.iahw)

   dat <- escalc(measure="ABT", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.iahw)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.iahw)

   dat <- escalc(measure="PCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   dat <- escalc(measure="ZPCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ztor)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ztor)

   dat <- escalc(measure="SPCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab)

   par(opar)

})
