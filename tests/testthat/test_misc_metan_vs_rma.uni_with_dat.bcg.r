### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma.uni() against metan with 'dat.bcg'")

source("tolerances.r") # read in tolerances

test_that("results match (FE model, measure='RR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph rr log

   res <- rma(yi, vi, data=dat, method="FE")

   expect_equivalent(round(c(res$beta), digits=4), -0.4303, tolerance=.tol[["coef"]])
   expect_equivalent(round(res$ci.lb,   digits=4), -0.5097, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$ci.ub,   digits=4), -0.3509, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$zval,    digits=4), -10.6247, tolerance=.tol[["test"]]) ### -10.62 in Stata
   expect_equivalent(round(res$QE,      digits=4), 152.2330, tolerance=.tol[["test"]])

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph rr

   sav <- predict(res, transf=exp)

   expect_equivalent(round(sav$pred,  digits=4), 0.6503, tolerance=.tol[["pred"]])
   expect_equivalent(round(sav$ci.lb, digits=4), 0.6007, tolerance=.tol[["ci"]])
   expect_equivalent(round(sav$ci.ub, digits=4), 0.7040, tolerance=.tol[["ci"]])

})

test_that("results match (RE model w/ DL estimator, measure='RR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph rr log

   res <- rma(yi, vi, data=dat, method="DL")

   expect_equivalent(round(c(res$beta), digits=4), -0.7141, tolerance=.tol[["coef"]])
   expect_equivalent(round(res$ci.lb,   digits=4), -1.0644, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$ci.ub,   digits=4), -0.3638, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$zval,    digits=4), -3.9952, tolerance=.tol[["test"]]) ### 4.00 in Stata
   expect_equivalent(round(res$tau2,    digits=4), 0.3088, tolerance=.tol[["var"]])
   expect_equivalent(round(res$I2,      digits=4), 92.1173, tolerance=.tol[["het"]])

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph rr

   sav <- predict(res, transf=exp)

   expect_equivalent(round(sav$pred,  digits=4), 0.4896, tolerance=.tol[["pred"]])
   expect_equivalent(round(sav$ci.lb, digits=4), 0.3449, tolerance=.tol[["ci"]])
   expect_equivalent(round(sav$ci.ub, digits=4), 0.6950, tolerance=.tol[["ci"]])

})

test_that("results match (FE model, measure='OR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph or log

   res <- rma(yi, vi, data=dat, method="FE")

   expect_equivalent(round(c(res$beta), digits=4), -0.4361, tolerance=.tol[["coef"]])
   expect_equivalent(round(res$ci.lb,   digits=4), -0.5190, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$ci.ub,   digits=4), -0.3533, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$zval,    digits=4), -10.3190, tolerance=.tol[["test"]]) ### -10.32 in Stata
   expect_equivalent(round(res$QE,      digits=4), 163.1649, tolerance=.tol[["test"]])

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph or

   sav <- predict(res, transf=exp)

   expect_equivalent(round(sav$pred,  digits=4), 0.6465, tolerance=.tol[["pred"]])
   expect_equivalent(round(sav$ci.lb, digits=4), 0.5951, tolerance=.tol[["ci"]])
   expect_equivalent(round(sav$ci.ub, digits=4), 0.7024, tolerance=.tol[["ci"]])

})

test_that("results match (RE model w/ DL estimator, measure='OR').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph or log

   res <- rma(yi, vi, data=dat, method="DL")

   expect_equivalent(round(c(res$beta), digits=4), -0.7474, tolerance=.tol[["coef"]])
   expect_equivalent(round(res$ci.lb,   digits=4), -1.1242, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$ci.ub,   digits=4), -0.3706, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$zval,    digits=4), -3.8873, tolerance=.tol[["test"]]) ### -3.89 in Stata
   expect_equivalent(round(res$tau2,    digits=4), 0.3663, tolerance=.tol[["var"]])
   expect_equivalent(round(res$I2,      digits=4), 92.6455, tolerance=.tol[["het"]])

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph or

   sav <- predict(res, transf=exp)

   expect_equivalent(round(sav$pred,  digits=4), 0.4736, tolerance=.tol[["pred"]])
   expect_equivalent(round(sav$ci.lb, digits=4), 0.3249, tolerance=.tol[["ci"]])
   expect_equivalent(round(sav$ci.ub, digits=4), 0.6903, tolerance=.tol[["ci"]])

})

test_that("results match (FE model, measure='RD').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, fixedi nograph rd

   res <- rma(yi, vi, data=dat, method="FE")

   expect_equivalent(round(c(res$beta), digits=4), -0.0009, tolerance=.tol[["coef"]])
   expect_equivalent(round(res$ci.lb,   digits=4), -0.0014, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$ci.ub,   digits=4), -0.0005, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$zval,    digits=4), -4.0448, tolerance=.tol[["test"]]) ### -4.04 in Stata
   expect_equivalent(round(res$QE,      digits=4), 276.4737, tolerance=.tol[["test"]])

})

test_that("results match (RE model w/ DL estimator, measure='RD').", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### compare results with: metan tpos tneg cpos cneg, randomi nograph rd

   res <- rma(yi, vi, data=dat, method="DL")

   expect_equivalent(round(c(res$beta), digits=4), -0.0071, tolerance=.tol[["coef"]])
   expect_equivalent(round(res$ci.lb,   digits=4), -0.0101, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$ci.ub,   digits=4), -0.0040, tolerance=.tol[["ci"]])
   expect_equivalent(round(res$zval,    digits=4), -4.5128, tolerance=.tol[["test"]]) ### -4.51 in Stata
   expect_equivalent(round(res$tau2,    digits=4), 0.0000, tolerance=.tol[["var"]])
   expect_equivalent(round(res$I2,      digits=4), 95.6596, tolerance=.tol[["het"]])

})

#expect_that(rma(yi ~ ablat, vi, data=dat, subset=1:2), throws_error("Number of parameters to be estimated is larger than the number of observations."))
