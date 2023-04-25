### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma.glmm() function")

source("settings.r")

dat <- dat.nielweise2007

test_that("rma.glmm() works correctly for 'UM.FS' model.", {

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", method="EE"))
   out <- capture.output(print(res))

   expect_equivalent(coef(res), -1.2286, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0, tolerance=.tol[["var"]])

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", test="t"))
   out <- capture.output(print(res))

   expect_equivalent(coef(res), -1.2370, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0.3198, tolerance=.tol[["var"]])

   ### check some (current) stop()'s

   expect_error(confint(res))
   expect_error(plot(res))
   expect_error(qqnorm(res))
   expect_error(weights(res))

   skip_on_cran()

   ### check GLMMadaptive and glmmTMB results

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", test="t", control=list(package="GLMMadaptive")))

   expect_equivalent(coef(res), -1.2370, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0.3227, tolerance=.tol[["var"]])

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", test="t", control=list(package="glmmTMB")))

   expect_equivalent(coef(res), -1.2372, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0.3312, tolerance=.tol[["var"]])

})

test_that("rma.glmm() works correctly for 'UM.RS' model.", {

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.RS", method="EE"))
   out <- capture.output(print(res))

   expect_equivalent(coef(res), -1.2207, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0, tolerance=.tol[["var"]])
   expect_equivalent(res$sigma2, 0.6155, tolerance=.tol[["var"]])

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.RS", test="t"))
   out <- capture.output(print(res))

   expect_equivalent(coef(res), -1.2812, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0.7258, tolerance=.tol[["var"]])
   expect_equivalent(res$sigma2, 0.5212, tolerance=.tol[["var"]])

   ### check some (current) stop()'s

   expect_error(confint(res))
   expect_error(plot(res))
   expect_error(qqnorm(res))
   expect_error(weights(res))

   skip_on_cran()

   ### check GLMMadaptive and glmmTMB results

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.RS", test="t", control=list(package="GLMMadaptive")))

   expect_equivalent(coef(res), -1.2795, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0.7301, tolerance=.tol[["var"]])
   expect_equivalent(res$sigma2, 0.5364, tolerance=.tol[["var"]])

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.RS", test="t", control=list(package="glmmTMB")))

   expect_equivalent(coef(res), -1.2812, tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2,   0.7258, tolerance=.tol[["var"]])
   expect_equivalent(res$sigma2, 0.5212, tolerance=.tol[["var"]])

})

test_that("rma.glmm() works correctly when using 'clogit' or 'clogistic'.", {

   skip_on_cran()

   expect_warning(res1 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", method="EE"))
   expect_warning(res2 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", method="EE", control=list(optimizer="clogit")))
   expect_warning(res3 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", method="EE", control=list(optimizer="clogistic")))

   expect_equivalent(coef(res1), -1.2236, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res2), -1.2236, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res3), -1.2236, tolerance=.tol[["coef"]])

   expect_equivalent(c(vcov(res1)), 0.0502, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res2)), 0.0502, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res3)), 0.0502, tolerance=.tol[["var"]])

})

test_that("rma.glmm() works correctly for 'CM.EL' model.", {

   skip_on_cran()

   expect_warning(res1  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL"))
   expect_warning(res2  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="Nelder-Mead", hessianCtrl=list(d=0.00001))))
   expect_warning(res3  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="BFGS")))
   expect_warning(res4  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="bobyqa")))
   expect_warning(res5  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="nloptr")))
   expect_warning(res6  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="hjk")))
   expect_warning(res7  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="nmk", hessianCtrl=list(r=4))))
   expect_warning(res8  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="mads", hessianCtrl=list(r=4))))
   expect_warning(res9  <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="ucminf", optCtrl=list(xtol=1e-6))))
   expect_warning(res10 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="lbfgsb3c")))
   expect_warning(res11 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="subplex", hessianCtrl=list(r=4))))
   expect_warning(res12 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL", control=list(optimizer="BBoptim")))

   expect_equivalent(coef(res1),  -1.3532, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res2),  -1.3540, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res3),  -1.3532, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res4),  -1.3532, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res5),  -1.3536, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res6),  -1.3532, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res7),  -1.3593, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res8),  -1.3542, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res9),  -1.3532, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res10), -1.3532, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res11), -1.3542, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res12), -1.3532, tolerance=.tol[["coef"]])

   expect_equivalent(c(vcov(res1)),  0.1232, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res2)),  0.1228, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res3)),  0.1232, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res4)),  0.1232, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res5)),  0.1231, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res6)),  0.1232, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res7)),  0.0413, tolerance=.tol[["var"]]) # :(
   expect_equivalent(c(vcov(res8)),  0.0405, tolerance=.tol[["var"]]) # :(
   expect_equivalent(c(vcov(res9)),  0.1232, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res10)), 0.1232, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res11)), 0.0405, tolerance=.tol[["var"]]) # :(
   expect_equivalent(c(vcov(res12)), 0.1233, tolerance=.tol[["var"]])

   expect_equivalent(res1$tau2,  0.6935, tolerance=.tol[["var"]])
   expect_equivalent(res2$tau2,  0.6945, tolerance=.tol[["var"]])
   expect_equivalent(res3$tau2,  0.6935, tolerance=.tol[["var"]])
   expect_equivalent(res4$tau2,  0.6935, tolerance=.tol[["var"]])
   expect_equivalent(res5$tau2,  0.6937, tolerance=.tol[["var"]])
   expect_equivalent(res6$tau2,  0.6935, tolerance=.tol[["var"]])
   expect_equivalent(res7$tau2,  0.7043, tolerance=.tol[["var"]])
   expect_equivalent(res8$tau2,  0.6944, tolerance=.tol[["var"]])
   expect_equivalent(res9$tau2,  0.6935, tolerance=.tol[["var"]])
   expect_equivalent(res10$tau2, 0.6935, tolerance=.tol[["var"]])
   expect_equivalent(res11$tau2, 0.6944, tolerance=.tol[["var"]])
   expect_equivalent(res12$tau2, 0.6935, tolerance=.tol[["var"]])

})

rm(list=ls())
