### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: location-scale models")

source("tolerances.r") # read in tolerances

dat <- dat.bangertdrowns2004

test_that("location-scale model works correctly for an intercept-only model", {

   res1 <- rma(yi, vi, data=dat)
   res2 <- rma.mv(yi, vi, random = ~ 1 | id, data=dat)
   res3 <- rma(yi, vi, data=dat, scale = ~ 1)
   res4 <- rma(yi, vi, data=dat, scale = res3$Z)

   expect_equivalent(res1$tau2, res2$sigma2, tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, exp(res3$alpha[1]), tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, exp(res4$alpha[1]), tolerance=.tol[["var"]])

})

test_that("location-scale model works correctly for two subgroups with different tau^2 values", {

   res1 <- rma.mv(yi, vi, data=dat, random = ~ factor(meta) | id, struct="DIAG", subset=!is.na(meta), control=list(hessian=TRUE, vctransf=TRUE))
   expect_warning(res2 <- rma(yi, vi, data=dat, scale = ~ meta))
   expect_warning(res3 <- rma(yi, vi, data=dat, scale = res2$Z.f))

   expect_equivalent(res1$tau2, c(exp(res2$alpha[1]), exp(res2$alpha[1] + res2$alpha[2])), tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, c(exp(res3$alpha[1]), exp(res3$alpha[1] + res3$alpha[2])), tolerance=.tol[["var"]])

   expect_warning(res4 <- rma(yi, vi, data=dat, scale = ~  0 + factor(meta)))

   expect_equivalent(unname(sqrt(diag(solve(res1$hessian[1:2, 1:2])))), res4$se.alpha, tolerance=.tol[["se"]])

})

test_that("profiling works correctly for location-scale models", {

   res1 <- rma(yi, vi, data=dat)
   sav1 <- profile(res1, plot=FALSE, progbar=FALSE, steps=2)

   res2 <- rma(yi, vi, data=dat, scale = ~ 1)
   sav2 <- profile(res2, xlim=log(range(sav1$tau2)), plot=FALSE, progbar=FALSE, steps=2)

   expect_equivalent(sav1$ll, sav2$ll, tolerance=.tol[["fit"]])

})

test_that("location-scale model works correctly for a continuous predictor", {

   res1 <- rma(yi, vi, data=dat, scale = ~ grade)
   expect_equivalent(res1$beta, 0.2220791, tolerance=.tol[["coef"]])
   expect_equivalent(res1$alpha, c(-3.10513013522415, 0.041361925354706), tolerance=.tol[["coef"]])

   res2 <- rma(yi, vi, data=dat, scale = ~ grade, link="identity")
   expect_equivalent(res1$tau2, res2$tau2, tolerance=.tol[["var"]])

   res3 <- rma.mv(yi, vi, data=dat, random = ~ sqrt(grade) | id, rho=0, struct="GEN", control=list(hessian=TRUE, vctransf=FALSE))

   expect_equivalent(c(res2$alpha), diag(res3$G), tolerance=.tol[["coef"]])
   expect_equivalent(diag(res2$M),  diag(res3$M), tolerance=.tol[["var"]])

   expect_equivalent(unname(sqrt(diag(solve(res3$hessian[1:2, 1:2])))), res2$se.alpha, tolerance=.tol[["se"]])

})

test_that("location-scale model works correctly for multiple predictors", {

   expect_warning(res1 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni)))
   expect_equivalent(res1$beta, 0.1110317, tolerance=.tol[["coef"]])
   expect_equivalent(res1$alpha, c(-1.08826059, -0.03429344, 2.09197456, -0.28439165), tolerance=.tol[["coef"]])

   expect_warning(res2 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(scaleZ=FALSE)))
   expect_equivalent(res2$beta, 0.1110317, tolerance=.tol[["coef"]])
   expect_equivalent(res2$alpha, c(-1.08826210, -0.03429332, 2.09197501, -0.28439156), tolerance=.tol[["coef"]])

})
