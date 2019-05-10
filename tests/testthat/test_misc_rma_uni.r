### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma() function")

source("tolerances.r") # read in tolerances

data(dat.bcg, package="metafor")
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

test_that("rma() correctly handles a formula for the 'yi' argument", {

   res1 <- rma(yi ~ ablat, vi, data=dat)
   res2 <- rma(yi, vi, mods = ~ ablat, data=dat)
   expect_equivalent(coef(res1), coef(res2))

})

test_that("rma() correctly handles an 'escalc' object", {

   res1 <- rma(yi, vi, data=dat)
   res2 <- rma(dat)
   expect_equivalent(coef(res1), coef(res2))

})

test_that("rma() works with method='DLIT' and method='SJIT'", {

   res <- rma(yi, vi, data=dat, method="DLIT", control=list(maxiter=500))
   expect_equivalent(res$tau2, 0.1576, tolerance=.tol[["var"]])
   res <- rma(yi, vi, data=dat, method="SJIT")
   expect_equivalent(res$tau2, 0.3181, tolerance=.tol[["var"]])

})

test_that("rma() works directly with input for measure='SMD'", {

   dat <- dat.normand1999
   dat <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   res1 <- rma(yi, vi, data=dat)
   res2 <- rma(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)

   expect_equivalent(res1$tau2, 1.0090, tolerance=.tol[["var"]])
   expect_equivalent(res2$tau2, 1.0090, tolerance=.tol[["var"]])

})

test_that("rma() works directly with input for measure='PCOR'", {

   ### data from Aloe and Thompson (2013)
   dat <- data.frame(ti = c(4.61, 6.19, 4.07, -0.77, 1.16),
                     ni = c(218, 232, 156, 382, 259),
                     mi = c(4, 7, 6, 19, 15),
                     r2i = c(.240, .455, .500, .327, .117))
   dat <- escalc(measure="PCOR", ti=ti, ni=ni, mi=mi, data=dat)
   res1 <- rma(yi, vi, data=dat)
   res2 <- rma(measure="PCOR", ti=ti, ni=ni, mi=mi, data=dat)

   expect_equivalent(res1$tau2, 0.0297, tolerance=.tol[["var"]])
   expect_equivalent(res2$tau2, 0.0297, tolerance=.tol[["var"]])

})

test_that("rma() works directly with input for measure='MN'", {

   dat <- dat.normand1999
   dat <- escalc(measure="MN", mi=m1i, sdi=sd1i, ni=n1i, data=dat)
   res1 <- rma(yi, vi, data=dat)
   res2 <- rma(measure="MN", mi=m1i, sdi=sd1i, ni=n1i, data=dat)

   expect_equivalent(res1$tau2, 408.9277, tolerance=.tol[["var"]])
   expect_equivalent(res2$tau2, 408.9277, tolerance=.tol[["var"]])

})

test_that("rma() works directly with input for measure='SMCR'", {

   datT <- data.frame(
   m_pre   = c(30.6, 23.5, 0.5, 53.4, 35.6),
   m_post  = c(38.5, 26.8, 0.7, 75.9, 36.0),
   sd_pre  = c(15.0, 3.1, 0.1, 14.5, 4.7),
   sd_post = c(11.6, 4.1, 0.1, 4.4, 4.6),
   ni      = c(20, 50, 9, 10, 14),
   ri      = c(.47, .64, .77, .89, .44))

   dat <- escalc(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ni=ni, ri=ri, data=datT)
   res1 <- rma(yi, vi, data=dat)
   res2 <- rma(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ni=ni, ri=ri, data=datT)

   expect_equivalent(res1$tau2, 0.3164, tolerance=.tol[["var"]])
   expect_equivalent(res2$tau2, 0.3164, tolerance=.tol[["var"]])

})

test_that("rma() works directly with input for measure='AHW'", {

   dat <- dat.bonett2010
   dat <- escalc(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat)
   res1 <- rma(yi, vi, data=dat)
   res2 <- rma(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat)

   expect_equivalent(res1$tau2, 0.0011, tolerance=.tol[["var"]])
   expect_equivalent(res2$tau2, 0.0011, tolerance=.tol[["var"]])

})
