### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma.glmm() function")

source("settings.r")

dat <- dat.nielweise2007

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

})

test_that("rma.glmm() works correctly when using 'clogit' or 'clogistic'.", {

   skip_on_cran()

   expect_warning(res1 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", method="EE"))
   expect_warning(res2 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", method="EE", control=list(optimizer="clogit")))
   expect_warning(res3 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", method="EE", control=list(optimizer="clogistic")))

   expect_equivalent(coef(res1), -1.2286, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res2), -1.2286, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res3), -1.2286, tolerance=.tol[["coef"]])

   expect_equivalent(c(vcov(res1)), 0.0504, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res2)), 0.0504, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res3)), 0.0504, tolerance=.tol[["var"]])

})

test_that("rma.glmm() works correctly when using 'nlminb' or 'minqa'.", {

   skip_on_cran()

   expect_warning(res1 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", method="ML"))
   expect_warning(res2 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", method="ML", control=list(optimizer="nlminb")))
   expect_warning(res3 <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.FS", method="ML", control=list(optimizer="bobyqa")))

   expect_equivalent(coef(res1), -1.2369, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res2), -1.2369, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res3), -1.2369, tolerance=.tol[["coef"]])

   expect_equivalent(c(vcov(res1)), 0.0786, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res2)), 0.0786, tolerance=.tol[["var"]])
   expect_equivalent(c(vcov(res3)), 0.0786, tolerance=.tol[["var"]])

})

rm(list=ls())
