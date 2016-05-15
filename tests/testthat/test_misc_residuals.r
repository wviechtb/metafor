### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: residuals() function")

test_that("residuals are correct for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma(yi, vi, data=dat, method="DL")
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(round(rstandard(res)$z, 4), c(0.1379, -1.0289, -0.4854, -1.1195, 1.7221, 0.5048))
   expect_equivalent(round(rstudent(res)$z, 4),  c(0.1394, -1.0279, -0.4796, -1.4694,  1.8682, 0.3783))

})

test_that("residuals are correct for rma.mh().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(round(rstandard(res)$z, 4), c(0.1068, -1.4399, -0.6173, -3.4733, 3.2377, 1.9749))
   expect_equivalent(round(rstudent(res)$z, 4),  c(0.1076, -1.4668, -0.6219, -4.2413, 3.3947, 2.7908))

})

test_that("residuals are correct for rma.peto().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="PETO", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(round(rstandard(res)$z, 4), c(0.2684, -1.1482, -0.4142, -2.3440, 3.4961, 0.8037))
   expect_equivalent(round(rstudent(res)$z, 4),  c(0.2705, -1.1700, -0.4173, -2.8891, 3.6614, 1.1391))

})

test_that("residuals are correct for rma.glmm().", {

   skip_on_cran()

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))

})
