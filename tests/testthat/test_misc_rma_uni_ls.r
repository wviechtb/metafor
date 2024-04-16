### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: rma() function with location-scale models")

source("settings.r")

test_that("location-scale model results are correct for in intercept-only model", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi, vi, data=dat, test="t")
   res2 <- rma(yi, vi, scale = ~ 1, data=dat, test="t", control=list(optimizer="optim"))
   res3 <- suppressWarnings(rma(yi, vi, scale = ~ 1, link="identity", data=dat, test="t", control=list(optimizer="Nelder-Mead")))
   expect_equivalent(res1$tau2, as.vector(exp(coef(res2)$alpha)), tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, as.vector(coef(res3)$alpha), tolerance=.tol[["var"]])

})

test_that("location-scale model results are correct for a categorical predictor", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi ~ alloc, vi, scale = ~ alloc - 1, data=dat)
   res2 <- rma(yi ~ alloc, vi, scale = ~ alloc - 1, link = "identity", data=dat, control=list(optimizer="solnp"))
   res3 <- rma.mv(yi ~ alloc, vi, random = ~ alloc | trial, struct="DIAG", data=dat, sparse=.sparse)
   expect_equivalent(as.vector(exp(coef(res1)$alpha)), as.vector(coef(res2)$alpha), tolerance=.tol[["var"]])
   expect_equivalent(as.vector(exp(coef(res1)$alpha)), res3$tau2, tolerance=.tol[["var"]])

})

test_that("location-scale model results are correct for a continuous predictor", {

   dat <- escalc(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.laopaiboon2015)

   dat$ni <- dat$n1i + dat$n2i
   dat$ni[dat$study == "Whitlock"] <- dat$ni[dat$study == "Whitlock"] + 2

   res <- suppressWarnings(rma(yi, vi, scale = ~ I(1/ni) - 1, link="identity", data=dat, method="ML"))
   expect_equivalent(as.vector(coef(res)$alpha), 79.07531, tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(0.8539, 0.5482, 1.3302), tolerance=.tol[["coef"]])

   res <- rma(yi, vi, scale = ~ I(1/ni), link="identity", data=dat, method="ML")
   expect_equivalent(as.vector(coef(res)$alpha), c(0.274623, 31.523043), tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(1.0161589, 0.6214663, 1.6615205), tolerance=.tol[["coef"]])

   res <- rma(yi, vi, scale = ~ I(1/ni) - 1, data=dat)
   expect_equivalent(as.vector(coef(res)$alpha), -34.5187, tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(1.1251, 0.6381, 1.9839), tolerance=.tol[["coef"]])

   res <- rma(yi, vi, scale = ~ I(1/ni), data=dat)
   expect_equivalent(as.vector(coef(res)$alpha), c(-0.8868, 42.4065), tolerance=.tol[["var"]])
   expect_equivalent(exp(c(res$beta, res$ci.lb, res$ci.ub)), c(1.0474, 0.6242, 1.7577), tolerance=.tol[["coef"]])

   sav <- coef(summary(res))

   expected <- list(beta = structure(list(estimate = 0.0463401794422422, se = 0.264116077624852, zval = 0.175453837793485, pval = 0.86072304016451, ci.lb = -0.471317820440453, ci.ub = 0.563998179324937), class = "data.frame", row.names = "intrcpt"), alpha = structure(list(estimate = c(-0.886827277584096, 42.4065282951426 ), se = c(1.23920300372018, 118.69324661881), zval = c(-0.715643260161388, 0.357278358315816), pval = c(0.474211654391012, 0.720883429839682 ), ci.lb = c(-3.31562053440951, -190.227960285855), ci.ub = c(1.54196597924132, 275.04101687614)), class = "data.frame", row.names = c("intrcpt", "I(1/ni)")))

   expect_equivalent(sav, expected, tolerance=.tol[["misc"]])

   sav <- model.matrix(res)$scale
   expect_equivalent(sav, cbind(1, 1/dat$ni))

   sav <- fitted(res)$scale
   expect_equivalent(sav, c(-0.4790722, -0.58818975, -0.8305852, -0.71086658, -0.49417424, -0.25389402, -0.66126064, -0.45847851, -0.54205875, -0.03869671, -0.03869671, -0.12956784, -0.40493491, -0.76426506, -0.35674567), tolerance=.tol[["var"]])

})

rm(list=ls())
