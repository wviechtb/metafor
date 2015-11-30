### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Comparing rma.uni() against direct computations")

test_that("results match (FE model).", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, mods = ~ ablat + year, data=dat, method="FE")

   X <- cbind(1, dat$ablat, dat$year)
   W <- diag(1/dat$vi)
   y <- cbind(dat$yi)

   b  <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
   vb <- solve(t(X) %*% W %*% X)

   expect_equivalent(res$b, b)
   expect_equivalent(res$vb, vb)

   yhat <- c(X %*% b)

   expect_equivalent(fitted(res), yhat)

   H <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

   expect_equivalent(hatvalues(res, type="matrix"), H)

   ei <- (diag(res$k) - H) %*% y

   expect_equivalent(resid(res), c(ei))

})
