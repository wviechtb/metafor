### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking misc: reporter() function")

source("settings.r")

test_that("reporter() works correctly for 'rma.uni' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   expect_error(res <- rma(yi, vi, data=dat), NA) # to avoid this being an empty test

   skip_on_cran()

   reporter(res, open=FALSE)

})

rm(list=ls())
