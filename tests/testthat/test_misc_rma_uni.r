### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma() function")

data(dat.bcg, package="metafor")
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

test_that("rma() correctly handles a formula for the 'yi' argument", {

   res1 <- rma(yi ~ ablat, vi, data=dat)
   res2 <- rma(yi, vi, mods = ~ ablat, data=dat)
   expect_equivalent(round(coef(res1), 4), round(coef(res2), 4))

})

test_that("rma() correctly handles an 'escalc' object", {

   res1 <- rma(yi, vi, data=dat)
   res2 <- rma(dat)
   expect_equivalent(round(coef(res1), 4), round(coef(res2), 4))

})

test_that("rma() works with method='DLIT' and method='SJIT'", {

   res <- rma(yi, vi, data=dat, method="DLIT", control=list(maxiter=500))
   expect_equivalent(round(res$tau2, 4), 0.1576)
   res <- rma(yi, vi, data=dat, method="SJIT")
   expect_equivalent(round(res$tau2, 4), 0.3181)

})
