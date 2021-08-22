### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: reporter() function")

source("tolerances.r") # read in tolerances

test_that("reporter() works correctly for 'rma.uni' objects.", {

   skip_on_cran()

   data(dat.bcg)
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)

   reporter(res, open=FALSE)

})
