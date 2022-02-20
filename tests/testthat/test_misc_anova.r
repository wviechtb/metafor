### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: anova() function")

source("settings.r")

test_that("anova() works correctly for comparing nested models.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res1 <- rma(yi, vi, data=dat, method="ML")
   res2 <- rma(yi ~ ablat, vi, data=dat, method="ML")
   sav <- anova(res1, res2)
   out <- capture.output(print(sav))

   expect_equivalent(sav$LRT, 9.9588, tolerance=.tol[["test"]])

   res1 <- rma(yi, vi, data=dat, method="REML")
   res2 <- rma(yi ~ ablat, vi, data=dat, method="REML")
   expect_warning(sav <- anova(res1, res2))

   expect_equivalent(sav$LRT, 8.2301, tolerance=.tol[["test"]])

})

test_that("anova() works correctly when using the 'btt' argument.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)
   sav <- anova(res, btt=3:4)
   out <- capture.output(print(sav))

   expect_equivalent(sav$QM, 1.2850, tolerance=.tol[["test"]])
   expect_equivalent(sav$QMp, 0.5260, tolerance=.tol[["pval"]])

   sav <- anova(res, btt="alloc")
   out <- capture.output(print(sav))

   expect_equivalent(sav$QM, 1.2850, tolerance=.tol[["test"]])
   expect_equivalent(sav$QMp, 0.5260, tolerance=.tol[["pval"]])

   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat, test="knha")
   sav <- anova(res, btt=3:4)
   out <- capture.output(print(sav))

   expect_equivalent(sav$QM, 0.6007, tolerance=.tol[["test"]])
   expect_equivalent(sav$QMp, 0.5690, tolerance=.tol[["pval"]])

})

test_that("anova() works correctly when using the 'X' argument.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)

   sav <- anova(res, X=rbind(c(1, 10, 0, 0), c(1, 30, 0, 0), c(1, 50, 0, 0)))
   out <- capture.output(print(sav))

   expect_equivalent(sav$zval, c(0.0588, -1.7964, -3.1210), tolerance=.tol[["test"]])

   sav <- anova(res, X=rbind(c(1, 10, 0, 0), c(1, 30, 0, 0), c(1, 50, 0, 0)), rhs=-.10)

   expect_equivalent(sav$zval, c(0.3463, -1.4543, -2.8295), tolerance=.tol[["test"]])

   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat, test="knha")
   sav <- anova(res, X=rbind(c(1, 10, 0, 0), c(1, 10, 1, 0), c(1, 10, 0, 1)))
   out <- capture.output(print(sav))

   expect_equivalent(sav$zval, c(0.0568, -0.8252, 0.2517), tolerance=.tol[["test"]])
   expect_equivalent(sav$QM, 0.4230, tolerance=.tol[["test"]])
   expect_equivalent(sav$QMp, 0.7412, tolerance=.tol[["pval"]])

})

rm(list=ls())
