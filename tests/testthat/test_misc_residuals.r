### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: residuals() function")

test_that("residuals are correct for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma(yi, vi, data=dat)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(round(rstandard(res)$z, 4), c(0.1401, -0.9930, -0.4719, -1.0475, 1.6462, 0.4825))
   expect_equivalent(round(rstudent(res)$z, 4),  c(0.1426, -0.9957, -0.4591, -1.1949, 2.0949, 0.4330))

   res <- rma(yi, vi, data=dat, method="FE")
   expect_equivalent(round(sum(residuals(res, type="pearson")^2), 4), round(res$QE, 4))
   expect_equivalent(round(sum(residuals(res, type="cholesky")^2), 4), round(res$QE, 4))

})

test_that("residuals are correct for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(round(rstandard(res)$z, 4), c(0.1401, -0.9930, -0.4719, -1.0476, 1.6462, 0.4825))

})

test_that("residuals are correct for rma.mh().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(round(residuals(res, type="rstandard"), 4), c(0.1068, -1.4399, -0.6173, -3.4733, 3.2377, 1.9749))
   expect_equivalent(round(residuals(res, type="rstudent"), 4),  c(0.1076, -1.4668, -0.6219, -4.2413, 3.3947, 2.7908))

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
