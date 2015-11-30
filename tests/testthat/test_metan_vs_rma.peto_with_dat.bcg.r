### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Comparing rma.peto() against metan with 'dat.bcg'")

test_that("results match (FE model, measure='OR').", {

   data(dat.bcg, package="metafor")

   ### compare results with: metan tpos tneg cpos cneg, peto nograph or log

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(round(res$b, digits=3), -0.474)
   expect_equivalent(round(res$ci.lb, digits=3), -0.554)
   expect_equivalent(round(res$ci.ub, digits=3), -0.395)
   expect_equivalent(round(res$zval,  digits=2), -11.67) ### 11.67 in Stata
   expect_equivalent(round(res$QE,    digits=2), 167.73)

   ### compare results with: metan tpos tneg cpos cneg, peto nograph or

   sav <- predict(res, transf=exp)

   expect_equivalent(round(sav$pred,  digits=3), 0.622)
   expect_equivalent(round(sav$ci.lb, digits=3), 0.575)
   expect_equivalent(round(sav$ci.ub, digits=3), 0.674)

})
