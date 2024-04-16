### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see also: https://www.metafor-project.org/doku.php/analyses:konstantopoulos2011

context("Checking analysis example: konstantopoulos2011")

source("settings.r")

dat <- dat.konstantopoulos2011

test_that("results are correct for the two-level random-effects model fitted with rma().", {

   res <- rma(yi, vi, data=dat)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(coef(res), 0.1279, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.0439, tolerance=.tol[["se"]])
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
   expect_equivalent(se(res), c(0.0440, 0.0044), tolerance=.tol[["se"]]) ### 0.043 in paper
   expect_equivalent(res$tau2, 0.0889, tolerance=.tol[["var"]]) ### 0.088 in paper
   expect_equivalent(res$se.tau2, 0.0205, tolerance=.tol[["sevar"]])

   ### CI for tau^2 based on the Q-profile method (CI in paper is based on a Satterthwaite approximation)
   tmp <- confint(res, digits=3)
   expect_equivalent(tmp$random[1,2], 0.0560, tolerance=.tol[["var"]])
   expect_equivalent(tmp$random[1,3], 0.1376, tolerance=.tol[["var"]])

})

test_that("results are correct for the two-level random-effects model fitted with rma.mv().", {

   res <- rma.mv(yi, vi, random = ~ 1 | study, data=dat, sparse=.sparse)

   ### compare with results on page 70 (Table 4)
   expect_equivalent(coef(res), 0.1279, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.0439, tolerance=.tol[["se"]])
   expect_equivalent(res$sigma2, 0.0884, tolerance=.tol[["var"]])

})

test_that("results are correct for the three-level random-effects model fitted with rma.mv() using ML estimation.", {

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, method="ML", sparse=.sparse)
   out <- capture.output(print(res.ml))
   out <- capture.output(print(summary(res.ml)))

   ### compare with results on page 71 (Table 5)
   expect_equivalent(coef(res.ml), 0.1845, tolerance=.tol[["coef"]])
   expect_equivalent(se(res.ml), 0.0805, tolerance=.tol[["se"]])
   expect_equivalent(res.ml$sigma2, c(0.0577, 0.0329), tolerance=.tol[["var"]])

   sav <- predict(res.ml)
   expect_equivalent(c(sav$pi.lb, sav$pi.ub), c(-0.4262, 0.7951), tolerance=.tol[["pred"]])

})

test_that("results are correct for the three-level mixed-effects model fitted with rma.mv() using ML estimation.", {

   ### three-level model (multilevel parameterization)
   res.ml <- rma.mv(yi, vi, mods = ~ I(year-mean(year)), random = ~ 1 | district/study, data=dat, method="ML", sparse=.sparse)
   out <- capture.output(print(res.ml))

   ### compare with results on page 71 (Table 5)
   expect_equivalent(coef(res.ml), c(0.1780, 0.0051), tolerance=.tol[["coef"]]) ### intercept is given as 0.183 in paper, but this seems to be a misprint
   expect_equivalent(se(res.ml), c(0.0805, 0.0085), tolerance=.tol[["se"]])
   expect_equivalent(res.ml$sigma2, c(0.0565, 0.0329), tolerance=.tol[["var"]])

})

test_that("results are correct for the three-level random-effects model fitted with rma.mv() using REML estimation.", {

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=.sparse)
   out <- capture.output(print(res.ml))

   ### (results for this not given in paper)
   expect_equivalent(coef(res.ml), 0.1847, tolerance=.tol[["coef"]])
   expect_equivalent(se(res.ml), 0.0846, tolerance=.tol[["se"]])
   expect_equivalent(res.ml$sigma2, c(0.0651, 0.0327), tolerance=.tol[["var"]])

   ### ICC
   expect_equivalent(res.ml$sigma2[1] / sum(res.ml$sigma2), 0.6653, tolerance=.tol[["cor"]])

   ### total amount of heterogeneity
   expect_equivalent(sum(res.ml$sigma2), 0.0978, tolerance=.tol[["var"]])

   ### log likelihood
   expect_equivalent(c(logLik(res.ml)), -7.9587, tolerance=.tol[["fit"]])

   ### CIs for variance components
   sav <- confint(res.ml)
   sav <- round(as.data.frame(sav), 4)
   expected <- structure(c(0.0651, 0.2551, 0.0327, 0.1809, 0.0222, 0.1491, 0.0163, 0.1276, 0.2072, 0.4552, 0.0628, 0.2507), .Dim = 4:3, .Dimnames = list(c("sigma^2.1", "sigma.1", "sigma^2.2", "sigma.2"), c("estimate", "ci.lb", "ci.ub")))
   expect_equivalent(sav, expected, tolerance=.tol[["var"]])

})

test_that("profiling works for the three-level random-effects model (multilevel parameterization).", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=.sparse)

   ### profile variance components
   png("images/test_analysis_example_konstantopoulos2011_profile_1_light_test.png", res=200, width=1800, height=2000, type="cairo")
   par(mfrow=c(2,1))
   sav <- profile(res.ml, progbar=FALSE)
   dev.off()

   expect_true(.vistest("images/test_analysis_example_konstantopoulos2011_profile_1_light_test.png", "images/test_analysis_example_konstantopoulos2011_profile_1_light.png"))

   ### profile variance components (dark theme)
   png("images/test_analysis_example_konstantopoulos2011_profile_1_dark_test.png", res=200, width=1800, height=2000, type="cairo")
   setmfopt(theme="dark")
   par(mfrow=c(2,1))
   sav <- profile(res.ml, progbar=FALSE)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_konstantopoulos2011_profile_1_dark_test.png", "images/test_analysis_example_konstantopoulos2011_profile_1_dark.png"))

   out <- capture.output(print(sav))

})

test_that("results are correct for the three-level random-effects model when using the multivariate parameterization.", {

   ### three-level model (mv = multivariate parameterization)
   res.mv <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat, sparse=.sparse)

   ### (results for this not given in paper)
   expect_equivalent(coef(res.mv), 0.1847, tolerance=.tol[["coef"]])
   expect_equivalent(se(res.mv), 0.0846, tolerance=.tol[["se"]])
   expect_equivalent(res.mv$tau2, 0.0978, tolerance=.tol[["var"]])
   expect_equivalent(res.mv$rho, 0.6653, tolerance=.tol[["cor"]])

   ### log likelihood
   expect_equivalent(c(logLik(res.mv)), -7.9587, tolerance=.tol[["fit"]])

})

test_that("profiling works for the three-level random-effects model (multivariate parameterization).", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### three-level model (mv = multivariate parameterization)
   res.mv <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat, sparse=.sparse)

   ### profile variance components
   png("images/test_analysis_example_konstantopoulos2011_profile_2_light_test.png", res=200, width=1800, height=2000, type="cairo")
   par(mfrow=c(2,1))
   #profile(res.mv, progbar=FALSE)
   profile(res.mv, progbar=FALSE, parallel="snow")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_konstantopoulos2011_profile_2_light_test.png", "images/test_analysis_example_konstantopoulos2011_profile_2_light.png"))

   ### profile variance components (dark theme)
   png("images/test_analysis_example_konstantopoulos2011_profile_2_dark_test.png", res=200, width=1800, height=2000, type="cairo")
   setmfopt(theme="dark")
   par(mfrow=c(2,1))
   #profile(res.mv, progbar=FALSE)
   profile(res.mv, progbar=FALSE, parallel="snow")
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_konstantopoulos2011_profile_2_dark_test.png", "images/test_analysis_example_konstantopoulos2011_profile_2_dark.png"))

})

test_that("BLUPs are calculated correctly for the three-level random-effects model (multilevel parameterization).", {

   skip_on_cran()

   ### three-level model (ml = multilevel parameterization)
   res.ml <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=.sparse)

   sav <- ranef(res.ml)

   expect_equivalent(sav[[1]]$intrcpt, c(-0.18998596, -0.08467077, 0.1407273, 0.24064814, -0.1072942, -0.23650899, 0.5342778, -0.2004695, 0.05711692, -0.14168396, -0.01215679), tolerance=.tol[["pred"]])
   expect_equivalent(sav[[1]]$se,      c(0.16653966, 0.12407891, 0.13724053, 0.11885896, 0.11895233, 0.10112845, 0.1297891, 0.101322, 0.11104458, 0.12485549, 0.15042221), tolerance=.tol[["se"]])
   expect_equivalent(sav[[2]]$intrcpt, c(-0.03794675, -0.04663383, 0.04357906, -0.05459167, 0.02098376, -0.25219111, 0.06169069, 0.12691378, 0.07315932, 0.02358293, -0.02593401, -0.16472466, 0.20017925, -0.05824454, 0.14387428, 0.00163316, -0.03082723, 0.09766431, -0.12245631, -0.07958353, 0.03342001, 0.03277405, -0.13648311, 0.00732233, -0.15120705, 0.10293055, 0.04267145, 0.08386343, -0.02323572, -0.03147411, -0.28733359, 0.19536367, 0.36079672, -0.0526358, -0.03322863, 0.00558571, 0.03469647, -0.01382146, 0.0152893, 0.02499288, -0.08174655, 0.19776024, 0.31299764, -0.03204218, -0.18968221, -0.13730492, -0.12298966, -0.28918454, 0.33743506, -0.03810734, 0.11843554, -0.19986832, -0.01436916, 0.12481101, -0.04350898, -0.07304968), tolerance=.tol[["pred"]])
   expect_equivalent(sav[[2]]$se,      c(0.16388194, 0.16388194, 0.16603559, 0.16603559, 0.12233812, 0.12233812, 0.12342216, 0.13171712, 0.13653182, 0.14617064, 0.12941105, 0.12588568, 0.10313659, 0.10313659, 0.10868276, 0.12489868, 0.10877088, 0.10517399, 0.10324522, 0.11803445, 0.11512181, 0.11661284, 0.12068892, 0.11803445, 0.11939164, 0.08878259, 0.09186319, 0.09186319, 0.09186319, 0.09186319, 0.12687757, 0.12311091, 0.12210943, 0.06873404, 0.06873404, 0.06873404, 0.06873404, 0.06873404, 0.06873404, 0.06873404, 0.06873404, 0.10744609, 0.10928134, 0.10744609, 0.10550931, 0.11267925, 0.11267925, 0.13697347, 0.13697347, 0.13632667, 0.13632667, 0.13632667, 0.1589217, 0.1581043, 0.15527374, 0.15527374), tolerance=.tol[["se"]])

})

test_that("restarting with 'restart=TRUE' works.", {

   skip_on_cran()

   expect_error(res <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, control=list(maxiter=4)))
   expect_error(res <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, control=list(maxiter=4), restart=TRUE))
   res <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, control=list(maxiter=4), restart=TRUE)

   expect_equivalent(coef(res), 0.1847132, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.08455592, tolerance=.tol[["se"]])
   expect_equivalent(res$sigma2, c(0.06506194, 0.03273652), tolerance=.tol[["var"]])

})

test_that("results are correct when allowing for different tau^2 per district.", {

   skip_on_cran()

   ### shuffle up dat to make sure that this does not affect things
   set.seed(1234)
   dat <- dat[sample(nrow(dat)),]

   res <- rma.mv(yi, vi, random = list(~ 1 | district, ~ factor(district) | study), struct="DIAG", data=dat, control=list(optimizer="optim"), sparse=.sparse)
   out <- capture.output(print(res, digits=4))
   out <- capture.output(print(summary(res, digits=4)))

   expect_equivalent(coef(res), 0.1270, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.0588, tolerance=.tol[["se"]])
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
