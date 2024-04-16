### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: rma.peto() against metan with 'dat.bcg'")

source("settings.r")

test_that("results match (EE model, measure='OR').", {

   ### compare results with: metan tpos tneg cpos cneg, peto nograph or log

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(res$beta,  -0.4744, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.5541, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, -0.3948, tolerance=.tol[["ci"]])
   expect_equivalent(res$zval,  -11.6689, tolerance=.tol[["test"]]) ### 11.67 in Stata
   expect_equivalent(res$QE,    167.7302, tolerance=.tol[["test"]])

   ### compare results with: metan tpos tneg cpos cneg, peto nograph or

   sav <- predict(res, transf=exp)

   expect_equivalent(sav$pred,  0.6222, tolerance=.tol[["pred"]])
   expect_equivalent(sav$ci.lb, 0.5746, tolerance=.tol[["ci"]])
   expect_equivalent(sav$ci.ub, 0.6738, tolerance=.tol[["ci"]])

})

rm(list=ls())
