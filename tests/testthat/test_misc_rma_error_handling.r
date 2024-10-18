### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking misc: proper handling of errors in rma()")

source("settings.r")

test_that("rma() handles NAs correctly.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   dat$yi[1] <- NA
   dat$yi[2] <- NA

   expect_warning(res <- rma(yi, vi, data=dat, digits=3))

   expect_equivalent(res$k, 11)
   expect_equivalent(res$k.f, 13)
   expect_equivalent(length(res$yi), 11)
   expect_equivalent(length(res$yi.f), 13)
   expect_equivalent(res$not.na, rep(c(FALSE,TRUE),times=c(2,11)))

   dat$ablat[3] <- NA

   ### TODO: complete this ...

})

rm(list=ls())
