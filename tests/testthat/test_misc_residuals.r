### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: residuals() function")

source("tolerances.r") # read in tolerances

test_that("residuals are correct for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma(yi, vi, data=dat)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(rstandard(res)$z, c(0.1401, -0.9930, -0.4719, -1.0475, 1.6462, 0.4825), tolerance=.tol[["pred"]])
   expect_equivalent(rstudent(res)$z,  c(0.1426, -0.9957, -0.4591, -1.1949, 2.0949, 0.4330), tolerance=.tol[["test"]])

   res <- rma(yi, vi, data=dat, method="EE")
   expect_equivalent(sum(residuals(res, type="pearson")^2), res$QE, tolerance=.tol[["test"]])
   expect_equivalent(sum(residuals(res, type="cholesky")^2), res$QE, tolerance=.tol[["test"]])

})

test_that("rstudent() yields the same results as a mean shift outlier model for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   dat$trial1 <- ifelse(dat$trial == 1, 1, 0)

   res <- rma(yi, vi, data=dat)
   sav <- rstudent(res)
   res <- rma(yi, vi, mods = ~ trial1, data=dat)

   expect_equivalent(coef(res)[2], sav$resid[1], tolerance=.tol[["coef"]])
   expect_equivalent(res$se[2], sav$se[1], tolerance=.tol[["se"]])

   res <- rma(yi, vi, data=dat, test="knha")
   sav <- rstudent(res)
   res <- rma(yi, vi, mods = ~ trial1, data=dat, test="knha")

   expect_equivalent(coef(res)[2], sav$resid[1], tolerance=.tol[["pred"]])
   expect_equivalent(res$se[2], sav$se[1], tolerance=.tol[["se"]])

})

test_that("residuals are correct for rma.mv().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(rstandard(res)$z, c(0.1401, -0.9930, -0.4719, -1.0476, 1.6462, 0.4825), tolerance=.tol[["test"]])
   expect_equivalent(rstandard(res, cluster=dat$alloc)$cluster$X2, c(3.7017, 3.6145), tolerance=.tol[["test"]])
   expect_equivalent(rstudent(res)$z, c(0.1426, -0.9957, -0.4591, -1.1949, 2.0949, 0.4330), tolerance=.tol[["test"]])
   expect_equivalent(rstudent(res, cluster=dat$alloc)$cluster$X2, c(27.4717, 5.2128), tolerance=.tol[["test"]])
   expect_equivalent(rstudent(res, cluster=dat$alloc, reestimate=FALSE)$cluster$X2, c(3.7017, 3.6145), tolerance=.tol[["test"]])

})

test_that("residuals are correct for rma.mh().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(residuals(res, type="rstandard"), c(0.1068, -1.4399, -0.6173, -3.4733, 3.2377, 1.9749), tolerance=.tol[["pred"]])
   expect_equivalent(residuals(res, type="rstudent"),  c(0.1076, -1.4668, -0.6219, -4.2413, 3.3947, 2.7908), tolerance=.tol[["pred"]])

})

test_that("residuals are correct for rma.peto().", {

   dat <- escalc(measure="PETO", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))
   expect_equivalent(rstandard(res)$z, c(0.2684, -1.1482, -0.4142, -2.3440, 3.4961, 0.8037), tolerance=.tol[["test"]])
   expect_equivalent(rstudent(res)$z,  c(0.2705, -1.1700, -0.4173, -2.8891, 3.6614, 1.1391), tolerance=.tol[["test"]])

})

test_that("residuals are correct for rma.glmm().", {

   skip_on_cran()

   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)

   res <- rma.glmm(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, subset=1:6)
   expect_equivalent(c(residuals(res)), c(dat$yi - coef(res)))

})
