### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: rma.uni() against direct computations")

source("settings.r")

test_that("results match (FE model).", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="FE")

   X <- cbind(1, dat$ablat, dat$year)
   W <- diag(1/dat$vi)
   y <- cbind(dat$yi)

   beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
   vb <- solve(t(X) %*% W %*% X)

   expect_equivalent(res$beta, beta)
   expect_equivalent(res$vb, vb)

   yhat <- c(X %*% beta)

   expect_equivalent(fitted(res), yhat)

   H <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

   expect_equivalent(hatvalues(res, type="matrix"), H)

   ei <- (diag(res$k) - H) %*% y

   expect_equivalent(resid(res), c(ei))

})

rm(list=ls())
