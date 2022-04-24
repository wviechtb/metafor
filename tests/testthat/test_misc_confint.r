### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: confint() function")

source("settings.r")

test_that("confint() works correctly for 'rma.uni' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat, method="DL")
   sav <- confint(res, fixed=TRUE, transf=exp)

   expect_equivalent(sav$fixed, c(0.4896, 0.3449, 0.6950), tolerance=.tol[["ci"]])
   expect_equivalent(sav$random[1,], c(0.3088, 0.1197, 1.1115), tolerance=.tol[["var"]])
   expect_equivalent(sav$random[3,], c(92.1173, 81.9177, 97.6781), tolerance=.tol[["het"]])
   expect_equivalent(sav$random[4,], c(12.6861, 5.5303, 43.0680), tolerance=.tol[["het"]])

   sav <- round(as.data.frame(sav), 4)
   expect_equivalent(sav[,1], c(0.4896, 0.3088, 0.5557, 92.1173, 12.6861))

})

test_that("confint() works correctly for 'rma.mh' objects.", {

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   sav <- confint(res, transf=exp)

   expect_equivalent(sav$fixed, c(0.6353, 0.5881, 0.6862), tolerance=.tol[["ci"]])

})

test_that("confint() works correctly for 'rma.peto' objects.", {

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   sav <- confint(res, transf=exp)

   expect_equivalent(sav$fixed, c(0.6222, 0.5746, 0.6738), tolerance=.tol[["ci"]])

})

rm(list=ls())
