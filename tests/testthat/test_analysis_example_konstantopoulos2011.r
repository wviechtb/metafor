### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011

context("Checking analysis example konstantopoulos2011")

dat <- get(data(dat.konstantopoulos2011, package="metafor"))

test_that("results are correct for the two-level random-effects model fitted with rma().", {

   res <- rma(yi, vi, data=dat)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(round(coef(res),3), 0.128)
   expect_equivalent(round(res$se,3), 0.044)
   expect_equivalent(round(res$tau2,3), 0.088)
   expect_equivalent(round(res$se.tau2,3), 0.020)

   ### CI for tau^2 based on the Q-profile method (CI in paper is based on a Satterthwaite approximation)
   tmp <- confint(res, digits=3)
   expect_equivalent(round(tmp$random[1,2],3), 0.056)
   expect_equivalent(round(tmp$random[1,3],3), 0.139)

})

test_that("results are correct for the two-level mixed-effects model fitted with rma().", {

   res <- rma(yi, vi, mods = ~ I(year-mean(year)), data=dat)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(round(coef(res),3), c(0.126, 0.005))
   expect_equivalent(round(res$se,3), c(0.044, 0.004)) ### 0.043 in paper
   expect_equivalent(round(res$tau2,3), 0.089) ### 0.088 in paper
   expect_equivalent(round(res$se.tau2,3), 0.020)

   ### CI for tau^2 based on the Q-profile method (CI in paper is based on a Satterthwaite approximation)
   tmp <- confint(res, digits=3)
   expect_equivalent(round(tmp$random[1,2],3), 0.056)
   expect_equivalent(round(tmp$random[1,3],3), 0.138)

})

test_that("results are correct for the two-level random-effects model fitted with rma.mv().", {

   res <- rma.mv(yi, vi, random = ~ 1 | study, data=dat)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(round(coef(res),3), 0.128)
   expect_equivalent(round(res$se,3), 0.044)
   expect_equivalent(round(res$sigma2,3), 0.088)

})

test_that("results are correct for the three-level random-effects model fitted with rma.mv() using ML estimation.", {

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, method="ML")

   ### compare with results on page 71 (Table 5)
   expect_equivalent(round(coef(res.ml),3), 0.184)
   expect_equivalent(round(res.ml$se,3), 0.080)
   expect_equivalent(round(res.ml$sigma2,3), c(0.058, 0.033))

})

test_that("results are correct for the three-level mixed-effects model fitted with rma.mv() using ML estimation.", {

   ### three-level model (multilevel parameterization)
   res.ml <- rma.mv(yi, vi, mods = ~ I(year-mean(year)), random = ~ 1 | district/study, data=dat, method="ML")

   ### compare with results on page 71 (Table 5)
   expect_equivalent(round(coef(res.ml),3), c(0.178, 0.005)) ### intercept is given as 0.183 in paper, but this seems to be a misprint
   expect_equivalent(round(res.ml$se,3), c(0.080, 0.009))
   expect_equivalent(round(res.ml$sigma2,3), c(0.056, 0.033))

})

test_that("results are correct for the three-level random-effects model fitted with rma.mv() using REML estimation.", {

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat)

   ### (results for this not given in paper)
   expect_equivalent(round(coef(res.ml),3), 0.185)
   expect_equivalent(round(res.ml$se,3), 0.085)
   expect_equivalent(round(res.ml$sigma2,3), c(0.065, 0.033))

   ### ICC
   expect_equivalent(round(res.ml$sigma2[1] / sum(res.ml$sigma2), 3), 0.665)

   ### total amount of heterogeneity
   expect_equivalent(round(sum(res.ml$sigma2), 3), 0.098)

   ### log likelihood
   expect_equivalent(round(c(logLik(res.ml)), 4), -7.9587)

})

test_that("profiling works for the three-level random-effects model (multilevel parameterization).", {

   skip_on_cran()

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat)

   ### profile variance components
   opar <- par()
   par(mfrow=c(2,1))
   profile(res.ml, sigma2=1, progbar=FALSE)
   profile(res.ml, sigma2=2, progbar=FALSE)
   par(opar)

})

test_that("results are correct for the three-level random-effects model when using the multivariate parameterization.", {

   ### three-level model (mv = multivariate parameterization)
   res.mv <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat)

   ### (results for this not given in paper)
   expect_equivalent(round(coef(res.mv),3), 0.185)
   expect_equivalent(round(res.mv$se,3), 0.085)
   expect_equivalent(round(res.mv$tau2,3), 0.098)
   expect_equivalent(round(res.mv$rho,3), 0.665)

   ### log likelihood
   expect_equivalent(round(c(logLik(res.mv)), 4), -7.9587)

})

test_that("profiling works for the three-level random-effects model (multivariate parameterization).", {

   skip_on_cran()

   ### three-level model (mv = multivariate parameterization)
   res.mv <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat)

   ### profile variance components
   opar <- par()
   par(mfrow=c(2,1))
   profile(res.mv, tau2=1, progbar=FALSE)
   profile(res.mv, rho=1, progbar=FALSE)
   par(opar)

})
