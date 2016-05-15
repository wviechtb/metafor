### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: transformation functions")

test_that("transformations work correctly.", {

   expect_equivalent(round(transf.rtoz(.5), 6), 0.549306)
   expect_equivalent(transf.ztor(transf.rtoz(.5)), .5)

   expect_equivalent(round(transf.logit(.1), 6), -2.197225)
   expect_equivalent(transf.ilogit(transf.logit(.1)), .1)

   expect_equivalent(round(transf.arcsin(.1), 6), 0.321751)
   expect_equivalent(transf.iarcsin(transf.arcsin(.1)), .1)

   expect_equivalent(round(transf.pft(.1,10), 6), 0.373394)
   expect_equivalent(transf.ipft(transf.pft(.1,10), 10), .1)
   expect_equivalent(transf.ipft.hm(transf.pft(.1,10), targs=list(ni=c(10))), .1)

   expect_equivalent(transf.isqrt(.1), 0.01)

   expect_equivalent(round(transf.irft(.1,10), 6), 0.381721)
   expect_equivalent(transf.iirft(transf.irft(.1,10), 10), .1)

   expect_equivalent(round(transf.ahw(.9), 6), 0.535841)
   expect_equivalent(transf.iahw(transf.ahw(.9)), .9)

   expect_equivalent(round(transf.abt(.9), 6), 2.302585)
   expect_equivalent(transf.iabt(transf.abt(.9)), .9)

   expect_equivalent(transf.ztor.int(transf.rtoz(.5), targs=list(tau2=0)), .5)
   expect_equivalent(round(transf.ztor.int(transf.rtoz(.5), targs=list(tau2=0.1)), 6), 0.46663)

   expect_equivalent(transf.exp.int(log(.5), targs=list(tau2=0)), .5)
   expect_equivalent(round(transf.exp.int(log(.5), targs=list(tau2=0.1)), 6), 0.525635)
   expect_equivalent(round(transf.exp.int(log(.5), targs=list(tau2=0.1, lower=-10, upper=10)), 6), round(exp(log(.5) + .1/2), 6))

   expect_equivalent(transf.ilogit.int(transf.logit(.1), targs=list(tau2=0)), .1)
   expect_equivalent(round(transf.ilogit.int(transf.logit(.1), targs=list(tau2=0.1)), 6), 0.103591)

})
