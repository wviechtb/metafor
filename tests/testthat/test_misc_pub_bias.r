### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: regtest() and ranktest() functions")

source("tolerances.r") # read in tolerances

test_that("regtest() works correctly for rma().", {

   dat <- dat.egger2001
   dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)
   res <- rma(yi, vi, data=dat)
   sav <- regtest(res)
   expect_equivalent(sav$zval, -4.6686, tolerance=.tol[["test"]])

   out <- capture.output(print(sav)) ### so that print.regtest.rma() is run (at least once)

   sav <- regtest(dat$yi, dat$vi)
   expect_equivalent(sav$zval, -4.6686, tolerance=.tol[["test"]])

   sav <- regtest(res, model="lm", predictor="sqrtninv")
   expect_equivalent(sav$zval, -5.6083, tolerance=.tol[["test"]])

   sav <- regtest(dat$yi, dat$vi, model="lm", predictor="sqrtninv")
   expect_equivalent(sav$zval, -5.6083, tolerance=.tol[["test"]])

})

test_that("ranktest() works correctly for rma().", {

   dat <- dat.egger2001
   dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)
   res <- rma(yi, vi, data=dat)
   sav <- ranktest(res)
   expect_equivalent(sav$tau, 0.15)
   expect_equivalent(sav$pval, 0.4503, tolerance=.tol[["pval"]])

   sav <- ranktest(dat$yi, dat$vi)
   expect_equivalent(sav$tau, 0.15)
   expect_equivalent(sav$pval, 0.4503, tolerance=.tol[["pval"]])

   out <- capture.output(print(sav)) ### so that print.ranktest.rma() is run (at least once)

})
