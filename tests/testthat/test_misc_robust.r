### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: robust() function")

test_that("robust() works correctly for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)

   sav <- robust(res, cluster=dat$trial)
   expect_equivalent(round(c(vcov(sav)), 6), 0.032106)

   sav <- robust(res, cluster=dat$trial, adjust=FALSE)
   expect_equivalent(round(c(vcov(sav)), 6), 0.029636)

   res <- rma(yi, vi, weights=1, data=dat)
   sav <- robust(res, cluster=dat$trial)
   expect_equivalent(round(c(vcov(sav)), 6), 0.037028)

})

test_that("robust() works correctly for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat)

   sav <- robust(res, cluster=dat$trial)
   expect_equivalent(round(c(vcov(sav)), 6), 0.032106)

   sav <- robust(res, cluster=dat$trial, adjust=FALSE)
   expect_equivalent(round(c(vcov(sav)), 6), 0.029636)

   res <- rma.mv(yi, vi, W=1, random = ~ 1 | trial, data=dat)
   sav <- robust(res, cluster=dat$trial)
   expect_equivalent(round(c(vcov(sav)), 6), 0.037028)

})
