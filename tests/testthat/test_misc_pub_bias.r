### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: regtest() and ranktest() functions")

test_that("regtest() works correctly for rma().", {

   dat <- get(data(dat.egger2001, package="metafor"))
   res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)
   sav <- regtest(res)
   expect_equivalent(round(sav$zval, 4), -4.6686)

   out <- capture.output(print(sav)) ### so that print.regtest.rma() is run (at least once)

   sav <- regtest(res, model="lm", predictor="sqrtninv")
   expect_equivalent(round(sav$zval, 4), -5.6083)

})

test_that("ranktest() works correctly for rma().", {

   dat <- get(data(dat.egger2001, package="metafor"))
   res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)
   sav <- ranktest(res)
   expect_equivalent(sav$tau, 0.15)
   expect_equivalent(round(sav$pval, 4), 0.4503)

   out <- capture.output(print(sav)) ### so that print.ranktest.rma() is run (at least once)

})
