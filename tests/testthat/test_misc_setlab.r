### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: .setlab() function")

source("settings.r")

yi <- c(-.3, -.1, 0, .2, .2)
vi <- rep(.02, length(yi))

test_that(".setlab() works correctly together with forest().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png(filename="images/test_misc_setlab_test.png", res=300, width=5000, height=8000, type="cairo")

   par(mfrow=c(14,6), mar=c(5,4,0,4))
   xlim <- c(-3,5)
   cex.lab <- 0.5

   dat <- escalc(measure="GEN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="RR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="OR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="RD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="AS", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="PHI", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="YUQ", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="YUY", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="IRR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="IRD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="IRSD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="MD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="SMD", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="ROM", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="CVR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="VR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="RPB", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="COR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="ZCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ztor, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ztor, header=TRUE)

   dat <- escalc(measure="PR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="PLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="PLO", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ilogit, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ilogit, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="PAS", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.iarcsin, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.iarcsin, header=TRUE)

   dat <- escalc(measure="PFT", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ipft.hm, targs=list(ni=rep(10,length(yi))), header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ipft.hm, targs=list(ni=rep(10,length(yi))), header=TRUE)

   dat <- escalc(measure="IR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="IRLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="IRS", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.isqrt, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.isqrt, header=TRUE)

   dat <- escalc(measure="IRFT", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="MN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="MNLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="CVLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="SDLN", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="MC", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="SMCC", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="ROMC", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=exp, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=exp, header=TRUE)

   dat <- escalc(measure="ARAW", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="AHW", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.iahw, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.iahw, header=TRUE)

   dat <- escalc(measure="ABT", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.iabt, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.iabt, header=TRUE)

   dat <- escalc(measure="PCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dat <- escalc(measure="ZPCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, transf=transf.ztor, header=TRUE)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, atransf=transf.ztor, header=TRUE)

   dat <- escalc(measure="SPCOR", yi=yi, vi=vi)
   forest(dat$yi, dat$vi, xlim=xlim, cex.lab=cex.lab, header=TRUE)

   dev.off()

   expect_true(.vistest("images/test_misc_setlab_test.png", "images/test_misc_setlab.png"))

})

rm(list=ls())
