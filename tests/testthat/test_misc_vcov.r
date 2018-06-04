### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: vcov() function")

test_that("vcov() works correctly for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi ~ ablat, vi, data=dat)

   expect_equivalent(round(vcov(res), 4), structure(c(0.0621, -0.0016, -0.0016, 1e-04), .Dim = c(2L, 2L), .Dimnames = list(c("intrcpt", "ablat"), c("intrcpt", "ablat"))))
   expect_equivalent(diag(vcov(res, type="obs")), dat$vi + res$tau2)
   expect_equivalent(round(vcov(res, type="fitted")[1,], 4), c(0.0197, 0.0269, 0.0184, 0.025, -0.0007, 0.0197, 0.0033, -0.0007, 0.0085, 0.0184, 0.0026, 0.0125, 0.0125))
   expect_equivalent(round(vcov(res, type="resid")[1,], 4), c(0.3822, -0.0269, -0.0184, -0.025, 7e-04, -0.0197, -0.0033, 0.0007, -0.0085, -0.0184, -0.0026, -0.0125, -0.0125))

})

test_that("vcov() works correctly for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi ~ ablat, vi, random = ~ 1 | trial, data=dat)

   expect_equivalent(round(vcov(res), 4), structure(c(0.062, -0.0016, -0.0016, 1e-04), .Dim = c(2L, 2L), .Dimnames = list(c("intrcpt", "ablat"), c("intrcpt", "ablat"))))
   expect_equivalent(diag(vcov(res, type="obs")), dat$vi + res$sigma2)
   expect_equivalent(round(vcov(res, type="fitted")[1,], 4), c(0.0197, 0.0269, 0.0184, 0.025, -0.0007, 0.0197, 0.0033, -0.0007, 0.0085, 0.0184, 0.0026, 0.0125, 0.0125))
   expect_equivalent(round(vcov(res, type="resid")[1,], 4), c(0.3822, -0.0269, -0.0184, -0.025, 7e-04, -0.0197, -0.0033, 0.0007, -0.0085, -0.0184, -0.0026, -0.0125, -0.0125))

})
