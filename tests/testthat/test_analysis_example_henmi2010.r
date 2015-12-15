### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:henmi2010

context("Checking analysis example henmi2010")

### load dataset
dat <- get(data(dat.lee2004, package="metafor"))

### calculate log odds ratios and corresponding sampling variances
dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)

test_that("results are correct for the random-effects model.", {

   ### fit random-effects model with DL estimator
   res <- rma(yi, vi, data=dat, method="DL")

   ### compare with results on page 2978
   expect_equivalent(round(res$tau2,3), 0.333)
   expect_equivalent(round(coef(res),3), -0.679)
   expect_equivalent(round(res$ci.lb,3), -1.066)
   expect_equivalent(round(res$ci.ub,3), -0.291)

})

test_that("results are correct for the Henmi & Copas method.", {

   ### fit random-effects model with DL estimator
   res <- rma(yi, vi, data=dat, method="DL")

   ### apply Henmi & Copas method
   sav <- hc(res)

   ### compare with results on page 2978
   expect_equivalent(round(sav$b,3), -0.514)
   expect_equivalent(round(sav$ci.lb,3), -0.999)
   expect_equivalent(round(sav$ci.ub,3), -0.030)

})
