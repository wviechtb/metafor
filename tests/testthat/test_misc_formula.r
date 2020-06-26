### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: formula() function")

source("tolerances.r") # read in tolerances

data(dat.bcg, package="metafor")
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

test_that("formula() works correctly for 'rma.uni' objects.", {

   res <- rma(yi, vi, data=dat, method="DL")
   expect_null(formula(res, type="mods"))
   expect_null(formula(res, type="yi"))

   res <- rma(yi, vi, mods = ~ ablat, data=dat, method="DL")
   expect_equal(~ablat, formula(res, type="mods"))
   expect_null(formula(res, type="yi"))

   res <- rma(yi ~ ablat, vi, data=dat, method="DL")
   expect_null(formula(res, type="mods"))
   expect_equal(yi~ablat, formula(res, type="yi"))

   expect_error(formula(res, type="scale"))

})
