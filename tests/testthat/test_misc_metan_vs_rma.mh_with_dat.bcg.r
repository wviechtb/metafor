### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma.mh() against metan with 'dat.bcg'")

source("settings.r")

test_that("results match (EE model, measure='RR').", {

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph rr log

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(res$beta,  -0.4537, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.5308, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, -0.3766, tolerance=.tol[["ci"]])
   expect_equivalent(res$zval,  -11.5338, tolerance=.tol[["test"]]) ### 11.53 in Stata
   expect_equivalent(res$QE,    152.5676, tolerance=.tol[["test"]])

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph rr

   sav <- predict(res, transf=exp)

   expect_equivalent(sav$pred,  0.6353, tolerance=.tol[["est"]])
   expect_equivalent(sav$ci.lb, 0.5881, tolerance=.tol[["ci"]])
   expect_equivalent(sav$ci.ub, 0.6862, tolerance=.tol[["ci"]])

})

test_that("results match (EE model, measure='OR').", {

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph or log

   res <- rma.mh(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(res$beta,  -0.4734, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.5538, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, -0.3930, tolerance=.tol[["ci"]])
   expect_equivalent(res$zval,  -11.5444, tolerance=.tol[["test"]]) ### 11.54 in Stata
   expect_equivalent(res$QE,    163.9426, tolerance=.tol[["test"]])

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph or

   sav <- predict(res, transf=exp)

   expect_equivalent(sav$pred,  0.6229, tolerance=.tol[["pred"]])
   expect_equivalent(sav$ci.lb, 0.5748, tolerance=.tol[["ci"]])
   expect_equivalent(sav$ci.ub, 0.6750, tolerance=.tol[["ci"]])

})

test_that("results match (EE model, measure='RD').", {

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph rd

   res <- rma.mh(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(res$beta,  -0.0033, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.0039, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, -0.0027, tolerance=.tol[["ci"]])
   expect_equivalent(res$zval,  -11.4708, tolerance=.tol[["test"]]) ### 11.56 in Stata
   expect_equivalent(res$QE,    386.7759, tolerance=.tol[["test"]])

   # zval is slightly different, as metan apparently computes the SE as
   # described in Greenland & Robins (1985) while metafor uses the equation
   # given in Sato, Greenland, & Robins (1989) (only the latter is
   # asymptotically correct in both the sparse-data and large-strata case)

})

rm(list=ls())
