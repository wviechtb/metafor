### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking tip: rma() results match up with those from lm() and lme()")

### this is essentially checking the equivalence of the results as explained here:
### http://www.metafor-project.org/doku.php/tips:rma_vs_lm_and_lme

test_that("results for rma() and lm() match for method='FE'.", {

   data(dat.molloy2014, package="metafor")

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)

   res.fe <- rma(yi, vi, data=dat, method="FE")
   res.lm <- lm(yi ~ 1, weights = 1/vi, data=dat)

   ### coefficients should be the same
   expect_equivalent(coef(res.fe), coef(res.lm))

   ### standard errors should be the same after adjusting the 'lm' one for sigma
   expect_equivalent(res.fe$se, coef(summary(res.lm))[1,2] / summary(res.lm)$sigma)

   ### fit the same model as is fitted by lm() with rma() function
   res.fe <- rma(yi, vi*summary(res.lm)$sigma^2, data=dat, method="FE")

   ### coefficients should still be the same
   expect_equivalent(coef(res.fe), coef(res.lm))

   ### standard errors should be the same
   expect_equivalent(res.fe$se, coef(summary(res.lm))[1,2])

})

test_that("results for rma() and lme() match for method='ML'.", {

   library("nlme")

   data(dat.molloy2014, package="metafor")

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
   dat$study <- 1:nrow(dat)
   res.lme <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), data=dat, method="ML")

   res.re <- rma(yi, vi*res.lme$sigma^2, data=dat, method="ML")

   ### coefficients should be the same
   expect_equivalent(round(coef(res.re),4), round(fixef(res.lme),4))

   ### standard errors should be the same after adjusting the 'rma' one by the factor sqrt(k/(k-p))
   expect_equivalent(round(res.re$se * sqrt(res.re$k / (res.re$k - res.re$p)),4), round(summary(res.lme)$tTable[1,2],4))

   ### check that BLUPs are the same
   expect_equivalent(round(blup(res.re)$pred,4), round(coef(res.lme)$"(Intercept)",4))

})

test_that("results for rma() and lme() match for method='REML'.", {

   library("nlme")

   data(dat.molloy2014, package="metafor")

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
   dat$study <- 1:nrow(dat)
   res.lme <- lme(yi ~ 1, random = ~ 1 | study, weights = varFixed(~ vi), data=dat, method="REML")

   res.re <- rma(yi, vi*res.lme$sigma^2, data=dat, method="REML")

   ### coefficients should be the same
   expect_equivalent(round(coef(res.re),4), round(fixef(res.lme),4))

   ### standard errors should be the same
   expect_equivalent(round(res.re$se,4), round(summary(res.lme)$tTable[1,2],4))

   ### check that BLUPs are the same
   expect_equivalent(round(blup(res.re)$pred,4), round(coef(res.lme)$"(Intercept)",4))

})
