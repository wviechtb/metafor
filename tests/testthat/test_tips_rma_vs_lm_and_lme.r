### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking tip: rma() results match up with those from lm() and lme()")

source("settings.r")

### this is essentially checking the equivalence of the results as explained here:
### https://www.metafor-project.org/doku.php/tips:rma_vs_lm_and_lme

test_that("results for rma() and lm() match for method='FE'.", {

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)

   res.ee <- rma(yi, vi, data=dat, method="EE")
   res.lm <- lm(yi ~ 1, weights = 1/vi, data=dat)

   ### coefficients should be the same
   expect_equivalent(coef(res.ee), coef(res.lm), tolerance=.tol[["coef"]])

   ### standard errors should be the same after adjusting the 'lm' one for sigma
   expect_equivalent(se(res.ee), se(res.lm) / sigma(res.lm), tolerance=.tol[["se"]])

   ### fit the same model as is fitted by lm() with rma() function
   res.ee <- rma(yi, vi*sigma(res.lm)^2, data=dat, method="EE")

   ### coefficients should still be the same
   expect_equivalent(coef(res.ee), coef(res.lm), tolerance=.tol[["coef"]])

   ### standard errors should be the same
   expect_equivalent(se(res.ee), se(res.lm), tolerance=.tol[["se"]])

})

test_that("results for rma() and lme() match for method='ML'.", {

   library("nlme")

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
   dat$study <- 1:nrow(dat)
   res.lme <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), data=dat, method="ML")

   res.re <- rma(yi, vi*sigma(res.lme)^2, data=dat, method="ML")

   ### coefficients should be the same
   expect_equivalent(coef(res.re), fixef(res.lme), tolerance=.tol[["coef"]])

   ### standard errors should be the same after adjusting the 'rma' one by the factor sqrt(k/(k-p))
   expect_equivalent(se(res.re) * sqrt(res.re$k / (res.re$k - res.re$p)), summary(res.lme)$tTable[1,2], tolerance=.tol[["se"]])

   ### check that BLUPs are the same
   expect_equivalent(blup(res.re)$pred, coef(res.lme)$"(Intercept)", tolerance=.tol[["pred"]])

})

test_that("results for rma() and lme() match for method='REML'.", {

   library("nlme")

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
   dat$study <- 1:nrow(dat)
   res.lme <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), data=dat, method="REML")

   res.re <- rma(yi, vi*sigma(res.lme)^2, data=dat, method="REML")

   ### coefficients should be the same
   expect_equivalent(coef(res.re), fixef(res.lme), tolerance=.tol[["coef"]])

   ### standard errors should be the same
   expect_equivalent(se(res.re), summary(res.lme)$tTable[1,2], tolerance=.tol[["se"]])

   ### check that BLUPs are the same
   expect_equivalent(blup(res.re)$pred, coef(res.lme)$"(Intercept)", tolerance=.tol[["pred"]])

})

rm(list=ls())
