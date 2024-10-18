### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

### see: https://www.metafor-project.org/doku.php/analyses:henmi2010

source("settings.r")

context("Checking analysis example: henmi2010")

### load dataset
dat <- dat.lee2004

### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)

test_that("results are correct for the random-effects model.", {

   ### fit random-effects model with DL estimator
   res <- rma(yi, vi, data=dat, method="DL")

   ### compare with results on page 2978
   expect_equivalent(res$tau2,   0.3325, tolerance=.tol[["var"]])
   expect_equivalent(coef(res), -0.6787, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -1.0664, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, -0.2911, tolerance=.tol[["ci"]])

})

test_that("results are correct for the Henmi & Copas method.", {

   ### fit random-effects model with DL estimator
   res <- rma(yi, vi, data=dat, method="DL")

   ### apply Henmi & Copas method
   sav <- hc(res)
   out <- capture.output(print(sav)) ### so that print.hc.rma.uni() is run (at least once)

   ### compare with results on page 2978
   expect_equivalent(sav$beta,  -0.5145, tolerance=.tol[["coef"]])
   expect_equivalent(sav$ci.lb, -0.9994, tolerance=.tol[["ci"]])
   expect_equivalent(sav$ci.ub, -0.0295, tolerance=.tol[["ci"]])

})

rm(list=ls())
