### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma() function with location-scale models")

source("settings.r")

test_that("location-scale model results are correct for in intercept-only model", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi, vi, data=dat, test="t")
   res2 <- rma(yi, vi, scale = ~ 1, data=dat, test="t", control=list(optimizer="optim"))
   res3 <- suppressWarnings(rma(yi, vi, scale = ~ 1, link="identity", data=dat, test="t", control=list(optimizer="optim", optmethod="Nelder-Mead")))
   expect_equivalent(res1$tau2, as.vector(exp(coef(res2)$alpha)), tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, as.vector(coef(res3)$alpha), tolerance=.tol[["var"]])

})

test_that("location-scale model results are correct for a categorical predictor", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi ~ alloc, vi, scale = ~ alloc - 1, data=dat)
   res2 <- rma(yi ~ alloc, vi, scale = ~ alloc - 1, link = "identity", data=dat)
   res3 <- rma.mv(yi ~ alloc, vi, random = ~ alloc | trial, struct="DIAG", data=dat, sparse=sparse)
   expect_equivalent(as.vector(exp(coef(res1)$alpha)), as.vector(coef(res2)$alpha), tolerance=.tol[["var"]])
   expect_equivalent(as.vector(exp(coef(res1)$alpha)), res3$tau2, tolerance=.tol[["var"]])

})

test_that("location-scale model results are correct for a continuous predictor", {

   dat <- escalc(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.laopaiboon2015)

   dat$ni <- dat$n1i + dat$n2i
   dat$ni[dat$study == "Whitlock"] <- dat$ni[dat$study == "Whitlock"] + 2

   res <- suppressWarnings(rma(yi, vi, scale = ~ I(1/ni) - 1, link="identity", data=dat, method="ML"))
   expect_equivalent(as.vector(coef(res)$alpha), 79.1084, tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(0.8539, 0.5482, 1.3302), tolerance=.tol[["coef"]])

   res <- rma(yi, vi, scale = ~ I(1/ni), link="identity", data=dat, method="ML")
   expect_equivalent(as.vector(coef(res)$alpha), c(0.2750, 31.5127), tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(1.0163, 0.6215, 1.6618), tolerance=.tol[["coef"]])

   res <- rma(yi, vi, scale = ~ I(1/ni) - 1, data=dat)
   expect_equivalent(as.vector(coef(res)$alpha), -34.5187, tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(1.1251, 0.6381, 1.9839), tolerance=.tol[["coef"]])

   res <- rma(yi, vi, scale = ~ I(1/ni), data=dat)
   expect_equivalent(as.vector(coef(res)$alpha), c(-0.8868, 42.4065), tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(1.0474, 0.6242, 1.7577), tolerance=.tol[["coef"]])

   sav <- coef(summary(res))

   expected <- list(beta = structure(list(estimate = 0.0463, se = 0.2641, zval = 0.1755, pval = 0.8607, ci.lb = -0.4713, ci.ub = 0.564), row.names = "intrcpt", class = "data.frame"), alpha = structure(list(estimate = c(-0.8868, 42.4065), se = c(1.2392, 118.6932), zval = c(-0.7156, 0.3573), pval = c(0.4742, 0.7209 ), ci.lb = c(-3.3156, -190.228), ci.ub = c(1.542, 275.041 )), row.names = c("intrcpt", "I(1/ni)"), class = "data.frame"))

   expect_equivalent(sav, expected, tolerance=.tol[["misc"]])

   sav <- model.matrix(res)$scale
   expect_equivalent(sav, cbind(1, 1/dat$ni))

   sav <- fitted(res)$scale
   expect_equivalent(sav, c(-0.479, -0.588, -0.831, -0.711, -0.494, -0.254, -0.661, -0.458, -0.542, -0.039, -0.039, -0.13, -0.405, -0.764, -0.357), tolerance=.tol[["var"]])

})

rm(list=ls())
