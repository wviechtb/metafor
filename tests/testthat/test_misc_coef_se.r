### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking misc: coef() and se() functions")

source("settings.r")

test_that("coef() and se() works correctly.", {

   dat <- dat.baskerville2012
   res <- rma(smd, se^2, data=dat, method="ML", digits=3)
   sel <- selmodel(res, type="beta")

   tmp <- list(beta = c(intrcpt = 0.114740253923052), delta = c(delta.1 = 0.473113053609697, delta.2 = 4.46131624677985))
   expect_equivalent(coef(sel)$beta,  tmp$beta,  tolerance=.tol[["coef"]])
   expect_equivalent(coef(sel)$delta, tmp$delta, tolerance=.tol[["coef"]])
   expect_equivalent(coef(sel, type="beta"),  tmp$beta,  tolerance=.tol[["coef"]])
   expect_equivalent(coef(sel, type="delta"), tmp$delta, tolerance=.tol[["coef"]])

   tmp <- list(beta = c(intrcpt = 0.166413798184622), delta = c(delta.1 = 0.235248084207613, delta.2 = 2.18419833595518))
   expect_equivalent(se(sel)$beta,  tmp$beta,  tolerance=.tol[["se"]])
   expect_equivalent(se(sel)$delta, tmp$delta, tolerance=.tol[["se"]])
   expect_equivalent(se(sel, type="beta"),  tmp$beta,  tolerance=.tol[["se"]])
   expect_equivalent(se(sel, type="delta"), tmp$delta, tolerance=.tol[["se"]])

   dat <- dat.bangertdrowns2004
   dat$ni100 <- dat$ni/100
   res <- rma(yi, vi, mods = ~ ni100, scale = ~ ni100, data=dat)

   tmp <- list(beta = c(intrcpt = 0.301681362709591, ni100 = -0.0552663301809239), alpha = c(intrcpt = -1.92087854601148, ni100 = -0.917428772771085))
   expect_equivalent(coef(res)$beta,  tmp$beta,  tolerance=.tol[["coef"]])
   expect_equivalent(coef(res)$alpha, tmp$alpha, tolerance=.tol[["coef"]])
   expect_equivalent(coef(res, type="beta"),  tmp$beta,  tolerance=.tol[["coef"]])
   expect_equivalent(coef(res, type="alpha"), tmp$alpha, tolerance=.tol[["coef"]])

   tmp <- list(beta = c(intrcpt = 0.0661161560867381, ni100 = 0.0197546220146866), alpha = c(intrcpt = 0.668982417863205, ni100 = 0.514064772257437))
   expect_equivalent(se(res)$beta,  tmp$beta,  tolerance=.tol[["se"]])
   expect_equivalent(se(res)$alpha, tmp$alpha, tolerance=.tol[["se"]])
   expect_equivalent(se(res, type="beta"),  tmp$beta,  tolerance=.tol[["se"]])
   expect_equivalent(se(res, type="alpha"), tmp$alpha, tolerance=.tol[["se"]])

})

rm(list=ls())
