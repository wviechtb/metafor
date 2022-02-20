### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:konstantopoulos2011

context("Checking analysis example: konstantopoulos2011")

source("settings.r")

dat <- dat.konstantopoulos2011

test_that("results are correct for the two-level random-effects model fitted with rma().", {

   res <- rma(yi, vi, data=dat)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(coef(res), 0.1279, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.0439, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.0884, tolerance=.tol[["var"]])
   expect_equivalent(res$se.tau2, 0.0202, tolerance=.tol[["sevar"]])

   ### CI for tau^2 based on the Q-profile method (CI in paper is based on a Satterthwaite approximation)
   tmp <- confint(res, digits=3)
   out <- capture.output(print(tmp)) ### so that print.confint.rma() is run (at least once)
   expect_equivalent(tmp$random[1,2], 0.0564, tolerance=.tol[["var"]])
   expect_equivalent(tmp$random[1,3], 0.1388, tolerance=.tol[["var"]])

})

test_that("results are correct for the two-level mixed-effects model fitted with rma().", {

   res <- rma(yi, vi, mods = ~ I(year-mean(year)), data=dat)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(coef(res), c(0.1258, 0.0052), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.0440, 0.0044), tolerance=.tol[["se"]]) ### 0.043 in paper
   expect_equivalent(res$tau2, 0.0889, tolerance=.tol[["var"]]) ### 0.088 in paper
   expect_equivalent(res$se.tau2, 0.0205, tolerance=.tol[["sevar"]])

   ### CI for tau^2 based on the Q-profile method (CI in paper is based on a Satterthwaite approximation)
   tmp <- confint(res, digits=3)
   expect_equivalent(tmp$random[1,2], 0.0560, tolerance=.tol[["var"]])
   expect_equivalent(tmp$random[1,3], 0.1376, tolerance=.tol[["var"]])

})

test_that("results are correct for the two-level random-effects model fitted with rma.mv().", {

   res <- rma.mv(yi, vi, random = ~ 1 | study, data=dat, sparse=sparse)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(coef(res), 0.1279, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.0439, tolerance=.tol[["se"]])
   expect_equivalent(res$sigma2, 0.0884, tolerance=.tol[["var"]])

})

test_that("results are correct for the three-level random-effects model fitted with rma.mv() using ML estimation.", {

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, method="ML", sparse=sparse)
   out <- capture.output(print(res.ml))
   out <- capture.output(print(summary(res.ml)))

   ### compare with results on page 71 (Table 5)
   expect_equivalent(coef(res.ml), 0.1845, tolerance=.tol[["coef"]])
   expect_equivalent(res.ml$se, 0.0805, tolerance=.tol[["se"]])
   expect_equivalent(res.ml$sigma2, c(0.0577, 0.0329), tolerance=.tol[["var"]])

   sav <- predict(res.ml)
   expect_equivalent(c(sav$pi.lb, sav$pi.ub), c(-0.4262, 0.7951), tolerance=.tol[["pred"]])

})

test_that("results are correct for the three-level mixed-effects model fitted with rma.mv() using ML estimation.", {

   ### three-level model (multilevel parameterization)
   res.ml <- rma.mv(yi, vi, mods = ~ I(year-mean(year)), random = ~ 1 | district/study, data=dat, method="ML", sparse=sparse)
   out <- capture.output(print(res.ml))

   ### compare with results on page 71 (Table 5)
   expect_equivalent(coef(res.ml), c(0.1780, 0.0051), tolerance=.tol[["coef"]]) ### intercept is given as 0.183 in paper, but this seems to be a misprint
   expect_equivalent(res.ml$se, c(0.0805, 0.0085), tolerance=.tol[["se"]])
   expect_equivalent(res.ml$sigma2, c(0.0565, 0.0329), tolerance=.tol[["var"]])

})

test_that("results are correct for the three-level random-effects model fitted with rma.mv() using REML estimation.", {

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=sparse)
   out <- capture.output(print(res.ml))

   ### (results for this not given in paper)
   expect_equivalent(coef(res.ml), 0.1847, tolerance=.tol[["coef"]])
   expect_equivalent(res.ml$se, 0.0846, tolerance=.tol[["se"]])
   expect_equivalent(res.ml$sigma2, c(0.0651, 0.0327), tolerance=.tol[["var"]])

   ### ICC
   expect_equivalent(res.ml$sigma2[1] / sum(res.ml$sigma2), 0.6653, tolerance=.tol[["cor"]])

   ### total amount of heterogeneity
   expect_equivalent(sum(res.ml$sigma2), 0.0978, tolerance=.tol[["var"]])

   ### log likelihood
   expect_equivalent(c(logLik(res.ml)), -7.9587, tolerance=.tol[["fit"]])

})

test_that("profiling works for the three-level random-effects model (multilevel parameterization).", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=sparse)

   ### profile variance components
   opar <- par(no.readonly=TRUE)
   par(mfrow=c(2,1))
   sav <- profile(res.ml, progbar=FALSE)
   out <- capture.output(print(sav))
   par(opar)

})

test_that("results are correct for the three-level random-effects model when using the multivariate parameterization.", {

   ### three-level model (mv = multivariate parameterization)
   res.mv <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat, sparse=sparse)

   ### (results for this not given in paper)
   expect_equivalent(coef(res.mv), 0.1847, tolerance=.tol[["coef"]])
   expect_equivalent(res.mv$se, 0.0846, tolerance=.tol[["se"]])
   expect_equivalent(res.mv$tau2, 0.0978, tolerance=.tol[["var"]])
   expect_equivalent(res.mv$rho, 0.6653, tolerance=.tol[["cor"]])

   ### log likelihood
   expect_equivalent(c(logLik(res.mv)), -7.9587, tolerance=.tol[["fit"]])

})

test_that("profiling works for the three-level random-effects model (multivariate parameterization).", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### three-level model (mv = multivariate parameterization)
   res.mv <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat, sparse=sparse)

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
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=sparse)

   sav <- ranef(res.ml)

   expect_equivalent(sav[[1]]$intrcpt, c(-0.190, -0.085,  0.141,  0.241, -0.107, -0.237,  0.534, -0.200,  0.057, -0.142, -0.012), tolerance=.tol[["pred"]])
   expect_equivalent(sav[[1]]$se,      c( 0.167,  0.124,  0.137,  0.119,  0.119,  0.101,  0.130,  0.101,  0.111,  0.125,  0.150), tolerance=.tol[["se"]])
   expect_equivalent(sav[[2]]$intrcpt, c(-0.038, -0.047,  0.044, -0.055,  0.021, -0.252,  0.062,  0.127,  0.073,  0.024, -0.026, -0.165, 0.200, -0.058, 0.144, 0.002, -0.031, 0.098, -0.122, -0.080, 0.033, 0.033, -0.136, 0.007, -0.151, 0.103, 0.043, 0.084, -0.023, -0.031, -0.287, 0.195, 0.361, -0.053, -0.033, 0.006, 0.035, -0.014, 0.015, 0.025, -0.082, 0.198, 0.313, -0.032, -0.19, -0.137, -0.123, -0.289, 0.337, -0.038, 0.118, -0.2, -0.014, 0.125, -0.044, -0.073), tolerance=.tol[["pred"]])
   expect_equivalent(sav[[2]]$se,      c( 0.164,  0.164,  0.166,  0.166,  0.122,  0.122,  0.123,  0.132,  0.137,  0.146,  0.129,  0.126, 0.103,  0.103, 0.109, 0.125,  0.109, 0.105,  0.103,  0.118, 0.115, 0.117, 0.121, 0.118, 0.119, 0.089, 0.092, 0.092, 0.092, 0.092, 0.127, 0.123, 0.122, 0.069, 0.069, 0.069, 0.069, 0.069, 0.069, 0.069, 0.069, 0.107, 0.109, 0.107, 0.106, 0.113, 0.113, 0.137, 0.137, 0.136, 0.136, 0.136, 0.159, 0.158, 0.155, 0.155), tolerance=.tol[["se"]])

})

test_that("results are correct when allowing for different tau^2 per district.", {

   skip_on_cran()

   ### shuffle up dat to make sure that this does not affect things
   set.seed(1234)
   dat <- dat[sample(nrow(dat)),]

   res <- rma.mv(yi, vi, random = list(~ 1 | district, ~ factor(district) | study), struct="DIAG", data=dat, control=list(optimizer="optim"), sparse=sparse)
   out <- capture.output(print(res, digits=4))
   out <- capture.output(print(summary(res, digits=4)))

   expect_equivalent(coef(res), 0.1270, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.0588, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, c(0.0000, 0.0402, 0.0000, 0.0582, 0.0082, 0.0000, 0.5380, 0.0008, 0.0606, 0.1803, 0.0000), tolerance=.tol[["var"]])

   ### check that output is also correct
   tau2 <- as.numeric(substr(out[grep("tau", out)], 13, 18))
   expect_equivalent(res$tau2, c(0.0000, 0.0402, 0.0000, 0.0582, 0.0082, 0.0000, 0.5380, 0.0008, 0.0606, 0.1803, 0.0000), tolerance=.tol[["var"]])
   k.lvl <- as.numeric(substr(out[grep("tau", out)], 32, 33))
   expect_equivalent(k.lvl, c(4, 4, 3, 4, 4, 11, 3, 8, 6, 5, 4))
   level <- as.numeric(substr(out[grep("tau", out)], 45, 47))
   expect_equivalent(level, c(11, 12, 18, 27, 56, 58, 71, 86, 91, 108, 644))

})

rm(list=ls())
