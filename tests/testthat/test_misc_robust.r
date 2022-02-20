### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: robust() function")

source("settings.r")

test_that("robust() works correctly for 'rma' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, data=dat)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.032106, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.98776, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, adjust=FALSE)
   expect_equivalent(c(vcov(sav)), 0.029636, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -4.150592, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, clubSandwich=TRUE)
   expect_equivalent(c(vcov(sav)), 0.03229357, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 11.04125, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.97616, tolerance=.tol[["test"]])

   res <- rma(yi, vi, weights=1, data=dat)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.037028, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.848996, tolerance=.tol[["test"]])

})

test_that("robust() works correctly for 'rma.mv' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, sparse=sparse)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.032106, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.98776, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, adjust=FALSE)
   expect_equivalent(c(vcov(sav)), 0.029636, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -4.150592, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, clubSandwich=TRUE)
   expect_equivalent(c(vcov(sav)), 0.03229357, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 11.04125, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.97616, tolerance=.tol[["test"]])

   res <- rma.mv(yi, vi, W=1, random = ~ 1 | trial, data=dat, sparse=sparse)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.037028, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.848996, tolerance=.tol[["test"]])

})

rm(list=ls())
