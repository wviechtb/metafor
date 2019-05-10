### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: update() function")

source("tolerances.r") # read in tolerances

test_that("update() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi, vi, data=dat, method="FE")
   res2 <- update(res1, method="DL")
   res3 <- rma(yi, vi, data=dat, method="DL")
   res4 <- update(res3, ~ ablat)
   res5 <- rma(yi, vi, mods = ~ ablat, data=dat, method="DL")
   expect_equivalent(res2, res3)
   expect_equivalent(res4, res5)

})

test_that("update() works for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma.mv(yi, vi, data=dat, method="FE")
   res2 <- update(res1, random = ~ 1 | trial, method="REML")
   res3 <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, method="REML")
   res4 <- update(res3, ~ ablat)
   res5 <- rma.mv(yi, vi, random = ~ 1 | trial, mods = ~ ablat, data=dat, method="REML")
   expect_equivalent(res2, res3)
   expect_equivalent(res4, res5)

})

test_that("update() works for rma.glmm().", {

   skip_on_cran()

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, method="FE")
   res2 <- update(res1, method="ML")
   res3 <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, method="ML")
   res4 <- update(res3, mods = ~ ablat)
   res5 <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, mods = ~ ablat, data=dat.bcg, method="ML")
   expect_equivalent(res2, res3)
   expect_equivalent(res4, res5)

})
