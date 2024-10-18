### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking misc: vcov() function")

source("settings.r")

test_that("vcov() works correctly for 'rma.uni' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi ~ ablat, vi, data=dat)

   expect_equivalent(vcov(res), structure(c(0.0621, -0.0016, -0.0016, 1e-04), .Dim = c(2L, 2L), .Dimnames = list(c("intrcpt", "ablat"), c("intrcpt", "ablat"))), tolerance=.tol[["var"]])
   expect_equivalent(diag(vcov(res, type="obs")), dat$vi + res$tau2)
   expect_equivalent(vcov(res, type="fitted")[1,], c(0.0197, 0.0269, 0.0184, 0.025, -0.0007, 0.0197, 0.0033, -0.0007, 0.0085, 0.0184, 0.0026, 0.0125, 0.0125), tolerance=.tol[["var"]])
   expect_equivalent(vcov(res, type="resid")[1,], c(0.3822, -0.0269, -0.0184, -0.025, 7e-04, -0.0197, -0.0033, 0.0007, -0.0085, -0.0184, -0.0026, -0.0125, -0.0125), tolerance=.tol[["var"]])

})

test_that("vcov() works correctly for 'rma.mv' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi ~ ablat, vi, random = ~ 1 | trial, data=dat, sparse=.sparse)

   expect_equivalent(vcov(res), structure(c(0.062, -0.0016, -0.0016, 1e-04), .Dim = c(2L, 2L), .Dimnames = list(c("intrcpt", "ablat"), c("intrcpt", "ablat"))), tolerance=.tol[["var"]])
   expect_equivalent(diag(vcov(res, type="obs")), dat$vi + res$sigma2)
   expect_equivalent(vcov(res, type="fitted")[1,], c(0.0197, 0.0269, 0.0184, 0.025, -0.0007, 0.0197, 0.0033, -0.0007, 0.0085, 0.0184, 0.0026, 0.0125, 0.0125), tolerance=.tol[["var"]])
   expect_equivalent(vcov(res, type="resid")[1,], c(0.3822, -0.0269, -0.0184, -0.025, 7e-04, -0.0197, -0.0033, 0.0007, -0.0085, -0.0184, -0.0026, -0.0125, -0.0125), tolerance=.tol[["var"]])

})

rm(list=ls())
