### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma.mh() against metan with 'dat.bcg'")

test_that("results match (FE model, measure='RR').", {

   data(dat.bcg, package="metafor")

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph rr log

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(round(res$beta,  digits=3), -0.454)
   expect_equivalent(round(res$ci.lb, digits=3), -0.531)
   expect_equivalent(round(res$ci.ub, digits=3), -0.377)
   expect_equivalent(round(res$zval,  digits=2), -11.53) ### 11.53 in Stata
   expect_equivalent(round(res$QE,    digits=2), 152.57)

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph rr

   sav <- predict(res, transf=exp)

   expect_equivalent(round(sav$pred,  digits=3), 0.635)
   expect_equivalent(round(sav$ci.lb, digits=3), 0.588)
   expect_equivalent(round(sav$ci.ub, digits=3), 0.686)

})

test_that("results match (FE model, measure='OR').", {

   data(dat.bcg, package="metafor")

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph or log

   res <- rma.mh(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(round(res$beta,  digits=3), -0.473)
   expect_equivalent(round(res$ci.lb, digits=3), -0.554)
   expect_equivalent(round(res$ci.ub, digits=3), -0.393)
   expect_equivalent(round(res$zval,  digits=2), -11.54) ### 11.54 in Stata
   expect_equivalent(round(res$QE,    digits=2), 163.94)

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph or

   sav <- predict(res, transf=exp)

   expect_equivalent(round(sav$pred,  digits=3), 0.623)
   expect_equivalent(round(sav$ci.lb, digits=3), 0.575)
   expect_equivalent(round(sav$ci.ub, digits=3), 0.675)

})

test_that("results match (FE model, measure='RD').", {

   data(dat.bcg, package="metafor")

   ### compare results with: metan tpos tneg cpos cneg, fixed nograph rd

   res <- rma.mh(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(round(res$beta,  digits=3), -0.003)
   expect_equivalent(round(res$ci.lb, digits=3), -0.004)
   expect_equivalent(round(res$ci.ub, digits=3), -0.003)
   expect_equivalent(round(res$zval,  digits=2), -11.47) ### 11.56 in Stata
   expect_equivalent(round(res$QE,    digits=2), 386.78)

   ### zval is slightly different, as metan apparently computes the SE as described in Greenland & Robins (1985)
   ### while metafor uses the equation given in Sato, Greenland, & Robins (1989) (only the latter is asymptotically
   ### correct in both the sparse-data and large-strata case)

})
