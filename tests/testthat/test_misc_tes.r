### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: tes() function")

source("tolerances.r") # read in tolerances

test_that("tes() works correctly for 'dat.dorn2007'.", {

   dat <- escalc(measure="RR", ai=x.a, n1i=n.a, ci=x.p, n2i=n.p, data=dat.dorn2007)

   sav <- tes(dat$yi, dat$vi, test="chi2")
   out <- capture.output(print(sav))

   expect_identical(sav$O, 10L)
   expect_equivalent(sav$E, 4.923333, tolerance=.tol[["misc"]])
   expect_equivalent(sav$X2, 7.065648, tolerance=.tol[["test"]])
   expect_equivalent(sav$pval, 0.003928794, tolerance=.tol[["pval"]])

   sav <- tes(dat$yi, dat$vi, test="binom")
   expect_equivalent(sav$pval, 0.01159554, tolerance=.tol[["pval"]])

   sav <- tes(dat$yi, dat$vi, test="exact", progbar=FALSE)
   expect_equivalent(sav$pval, 0.007778529, tolerance=.tol[["pval"]])

   res <- rma(yi, vi, data=dat, method="FE")
   sav <- tes(res, test="chi2")

   expect_identical(sav$O, 10L)
   expect_equivalent(sav$E, 4.923333, tolerance=.tol[["misc"]])
   expect_equivalent(sav$X2, 7.065648, tolerance=.tol[["test"]])
   expect_equivalent(sav$pval, 0.003928794, tolerance=.tol[["pval"]])

})
