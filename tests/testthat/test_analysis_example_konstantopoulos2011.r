### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011

context("Checking analysis example: konstantopoulos2011")

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
   out <- capture.output(print(tmp)) ### so that print.confint.rma() is run (at least once)
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
   out <- capture.output(print(res.ml))
   out <- capture.output(print(summary(res.ml)))

   ### compare with results on page 71 (Table 5)
   expect_equivalent(round(coef(res.ml),3), 0.184)
   expect_equivalent(round(res.ml$se,3), 0.080)
   expect_equivalent(round(res.ml$sigma2,3), c(0.058, 0.033))

   sav <- predict(res.ml)
   expect_equivalent(round(c(sav$cr.lb, sav$cr.ub),3), c(-0.426, 0.795))

})

test_that("results are correct for the three-level mixed-effects model fitted with rma.mv() using ML estimation.", {

   ### three-level model (multilevel parameterization)
   res.ml <- rma.mv(yi, vi, mods = ~ I(year-mean(year)), random = ~ 1 | district/study, data=dat, method="ML")
   out <- capture.output(print(res.ml))

   ### compare with results on page 71 (Table 5)
   expect_equivalent(round(coef(res.ml),3), c(0.178, 0.005)) ### intercept is given as 0.183 in paper, but this seems to be a misprint
   expect_equivalent(round(res.ml$se,3), c(0.080, 0.009))
   expect_equivalent(round(res.ml$sigma2,3), c(0.056, 0.033))

})

test_that("results are correct for the three-level random-effects model fitted with rma.mv() using REML estimation.", {

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat)
   out <- capture.output(print(res.ml))

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
   opar <- par(no.readonly=TRUE)
   par(mfrow=c(2,1))
   profile(res.ml, progbar=FALSE)
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
   opar <- par(no.readonly=TRUE)
   par(mfrow=c(2,1))
   #profile(res.mv, progbar=FALSE)
   profile(res.mv, progbar=FALSE, parallel="snow")
   par(opar)

})

test_that("BLUPs are calculated correctly for the three-level random-effects model (multilevel parameterization).", {

   skip_on_cran()

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat)

   sav <- lapply(ranef(res.ml), round, 3)

   expect_equivalent(sav[[1]]$intrcpt, c(-0.190, -0.085,  0.141,  0.241, -0.107, -0.237,  0.534, -0.200,  0.057, -0.142, -0.012))
   expect_equivalent(sav[[1]]$se,      c( 0.167,  0.124,  0.137,  0.119,  0.119,  0.101,  0.130,  0.101,  0.111,  0.125,  0.150))
   expect_equivalent(sav[[2]]$intrcpt, c(-0.038, -0.047, 0.044, -0.055, 0.021, -0.252, 0.062, 0.127, 0.073, 0.024, -0.026, -0.165, 0.2, -0.058, 0.144, 0.002, -0.031, 0.098, -0.122, -0.08, 0.033, 0.033, -0.136, 0.007, -0.151, 0.103, 0.043, 0.084, -0.023, -0.031, -0.287, 0.195, 0.361, -0.053, -0.033, 0.006, 0.035, -0.014, 0.015, 0.025, -0.082, 0.198, 0.313, -0.032, -0.19, -0.137, -0.123, -0.289, 0.337, -0.038, 0.118, -0.2, -0.014, 0.125, -0.044, -0.073))
   expect_equivalent(sav[[2]]$se,      c(0.164, 0.164, 0.166, 0.166, 0.122, 0.122, 0.123, 0.132, 0.137, 0.146, 0.129, 0.126, 0.103, 0.103, 0.109, 0.125, 0.109, 0.105, 0.103, 0.118, 0.115, 0.117, 0.121, 0.118, 0.119, 0.089, 0.092, 0.092, 0.092, 0.092, 0.127, 0.123, 0.122, 0.069, 0.069, 0.069, 0.069, 0.069, 0.069, 0.069, 0.069, 0.107, 0.109, 0.107, 0.106, 0.113, 0.113, 0.137, 0.137, 0.136, 0.136, 0.136, 0.159, 0.158, 0.155, 0.155))

})
