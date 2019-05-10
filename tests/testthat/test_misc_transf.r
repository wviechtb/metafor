### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: transformation functions")

source("tolerances.r") # read in tolerances

test_that("transformations work correctly.", {

   expect_equivalent(transf.rtoz(.5), 0.549306, tolerance=.tol[["est"]])
   expect_equivalent(transf.ztor(transf.rtoz(.5)), .5)

   expect_equivalent(transf.logit(.1), -2.197225, tolerance=.tol[["est"]])
   expect_equivalent(transf.ilogit(transf.logit(.1)), .1)

   expect_equivalent(transf.arcsin(.1), 0.321751, tolerance=.tol[["est"]])
   expect_equivalent(transf.iarcsin(transf.arcsin(.1)), .1)

   expect_equivalent(transf.pft(.1,10), 0.373394, tolerance=.tol[["est"]])
   expect_equivalent(transf.ipft(transf.pft(.1,10), 10), .1)
   expect_equivalent(transf.ipft.hm(transf.pft(.1,10), targs=list(ni=c(10))), .1)

   expect_equivalent(transf.isqrt(.1), 0.01)

   expect_equivalent(transf.irft(.1,10), 0.381721, tolerance=.tol[["est"]])
   expect_equivalent(transf.iirft(transf.irft(.1,10), 10), .1)

   expect_equivalent(transf.ahw(.9), 0.535841, tolerance=.tol[["est"]])
   expect_equivalent(transf.iahw(transf.ahw(.9)), .9)

   expect_equivalent(transf.abt(.9), 2.302585, tolerance=.tol[["est"]])
   expect_equivalent(transf.iabt(transf.abt(.9)), .9)

   expect_equivalent(transf.ztor.int(transf.rtoz(.5), targs=list(tau2=0)), .5)
   expect_equivalent(transf.ztor.int(transf.rtoz(.5), targs=list(tau2=0.1)), 0.46663, tolerance=.tol[["est"]])

   expect_equivalent(transf.exp.int(log(.5), targs=list(tau2=0)), .5)
   expect_equivalent(transf.exp.int(log(.5), targs=list(tau2=0.1)), 0.525635, tolerance=.tol[["est"]])
   expect_equivalent(transf.exp.int(log(.5), targs=list(tau2=0.1, lower=-10, upper=10)), round(exp(log(.5) + 0.1/2), 6), tolerance=.tol[["est"]])

   expect_equivalent(transf.ilogit.int(transf.logit(.1), targs=list(tau2=0)), .1)
   expect_equivalent(transf.ilogit.int(transf.logit(.1), targs=list(tau2=0.1)), 0.103591, tolerance=.tol[["est"]])

})
