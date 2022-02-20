### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: update() function")

source("settings.r")

test_that("update() works for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi, vi, data=dat, method="EE")
   res2 <- update(res1, method="DL")
   res3 <- rma(yi, vi, data=dat, method="DL")
   res4 <- update(res3, ~ ablat)
   res5 <- rma(yi, vi, mods = ~ ablat, data=dat, method="DL")
   res2$time <- NULL
   res3$time <- NULL
   res4$time <- NULL
   res5$time <- NULL
   expect_equivalent(res2, res3)
   expect_equivalent(res4, res5)

})

test_that("update() works for rma.mv().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma.mv(yi, vi, data=dat, method="EE", sparse=sparse)
   res2 <- update(res1, random = ~ 1 | trial, method="REML")
   res3 <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, method="REML", sparse=sparse)
   res4 <- update(res3, ~ ablat)
   res5 <- rma.mv(yi, vi, random = ~ 1 | trial, mods = ~ ablat, data=dat, method="REML", sparse=sparse)
   res2$time <- NULL
   res3$time <- NULL
   res4$time <- NULL
   res5$time <- NULL
   expect_equivalent(res2, res3)
   expect_equivalent(res4, res5)

})

test_that("update() works for rma.glmm().", {

   skip_on_cran()

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, method="EE")
   res2 <- update(res1, method="ML")
   res3 <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, method="ML")
   res4 <- update(res3, mods = ~ ablat)
   res5 <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, mods = ~ ablat, data=dat.bcg, method="ML")
   res2$time <- NULL
   res3$time <- NULL
   res4$time <- NULL
   res5$time <- NULL
   expect_equivalent(res2, res3)
   expect_equivalent(res4, res5)

})

rm(list=ls())
