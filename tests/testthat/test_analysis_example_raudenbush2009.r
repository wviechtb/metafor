### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:raudenbush2009

context("Checking analysis example: raudenbush2009")

source("tolerances.r") # read in tolerances

### load data
dat <- dat.raudenbush1985

test_that("results are correct for the fixed-effects model.", {

   ### fixed-effects model
   res.FE <- rma(yi, vi, data=dat, digits=3, method="FE")

   ### compare with results on page 301 (Table 16.2) and page 302
   expect_equivalent(coef(res.FE), 0.0604, tolerance=.tol[["coef"]])
   expect_equivalent(res.FE$QE,   35.8295, tolerance=.tol[["test"]])
   expect_equivalent(res.FE$zval,  1.6553, tolerance=.tol[["test"]]) ### 1.65 in chapter

})

test_that("results are correct for the random-effects model.", {

   ### random-effects model
   res.RE <- rma(yi, vi, data=dat, digits=3)

   ### compare with results on page 301 (Table 16.2) and page 302
   expect_equivalent(coef(res.RE), 0.0837, tolerance=.tol[["coef"]]) ### 0.083 in chapter
   expect_equivalent(res.RE$zval,  1.6208, tolerance=.tol[["test"]])
   expect_equivalent(res.RE$tau2,  0.0188, tolerance=.tol[["var"]])

   ### prediction interval
   tmp <- predict(res.RE)

   ### compare with results on page 301 (Table 16.2) and page 302
   expect_equivalent(tmp$pi.lb, -0.2036, tolerance=.tol[["ci"]]) ### -0.19 in chapter but computed in a slightly different way
   expect_equivalent(tmp$pi.ub,  0.3711, tolerance=.tol[["ci"]]) ###  0.35 in chapter but computed in a slightly different way

   ### range of BLUPs
   tmp <- range(blup(res.RE)$pred)

   ### compare with results on page 301 (Table 16.2)
   expect_equivalent(tmp, c(-0.0293, 0.2485), tolerance=.tol[["pred"]])

})

test_that("results are correct for the mixed-effects model.", {

   ### recode weeks variable
   dat$weeks.c <- ifelse(dat$weeks > 3, 3, dat$weeks)

   ### mixed-effects model
   res.ME <- rma(yi, vi, mods = ~ weeks.c, data=dat, digits=3)

   ### compare with results on page 301 (Table 16.2)
   expect_equivalent(res.ME$tau2, 0.0000, tolerance=.tol[["var"]])
   expect_equivalent(coef(res.ME), c(0.4072, -0.1572), tolerance=.tol[["coef"]])
   expect_equivalent(res.ME$QE, 16.5708, tolerance=.tol[["test"]])
   expect_equivalent(res.ME$zval, c(4.6782, -4.3884), tolerance=.tol[["test"]])

   ### range of BLUPs
   tmp <- range(blup(res.ME)$pred)

   ### compare with results on page 301 (Table 16.2)
   expect_equivalent(tmp, c(-0.0646, 0.4072), tolerance=.tol[["pred"]]) ### -0.07 in chapter

})

test_that("results are correct for the random-effects model (conventional approach).", {

   res.std <- list()

   res.std$FE   <- rma(yi, vi, data=dat, digits=3, method="FE")
   res.std$ML   <- rma(yi, vi, data=dat, digits=3, method="ML")
   res.std$REML <- rma(yi, vi, data=dat, digits=3, method="REML")
   res.std$DL   <- rma(yi, vi, data=dat, digits=3, method="DL")
   res.std$HE   <- rma(yi, vi, data=dat, digits=3, method="HE")

   tmp <- t(sapply(res.std, function(x) c(tau2=x$tau2, mu=x$beta, se=x$se, z=x$zval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)))

   expected <- structure(c(0, 0.0126, 0.0188, 0.0259, 0.0804, 0.0604, 0.0777, 0.0837, 0.0893, 0.1143, 0.0365, 0.0475, 0.0516, 0.0558, 0.0792, 1.6553, 1.6368, 1.6208, 1.6009, 1.4432, -0.0111, -0.0153, -0.0175, -0.02, -0.0409, 0.1318, 0.1708, 0.1849, 0.1987, 0.2696),
                         .Dim = 5:6, .Dimnames = list(c("FE", "ML", "REML", "DL", "HE"), c("tau2", "mu", "se", "z", "ci.lb", "ci.ub")))

   ### compare with results on page 309 (Table 16.3)
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

})

test_that("results are correct for the random-effects model (Knapp & Hartung method).", {

   res.knha <- list()

   expect_warning(res.knha$FE <- rma(yi, vi, data=dat, digits=3, method="FE", test="knha"))
   res.knha$ML   <- rma(yi, vi, data=dat, digits=3, method="ML", test="knha")
   res.knha$REML <- rma(yi, vi, data=dat, digits=3, method="REML", test="knha")
   res.knha$DL   <- rma(yi, vi, data=dat, digits=3, method="DL", test="knha")
   res.knha$HE   <- rma(yi, vi, data=dat, digits=3, method="HE", test="knha")

   tmp <- t(sapply(res.knha, function(x) c(tau2=x$tau2, mu=x$beta, se=x$se, z=x$zval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)))

   expected <- structure(c(0, 0.0126, 0.0188, 0.0259, 0.0804, 0.0604, 0.0777, 0.0837, 0.0893, 0.1143, 0.0515, 0.0593, 0.0616, 0.0636, 0.0711, 1.1733, 1.311, 1.3593, 1.405, 1.6078, -0.0477, -0.0468, -0.0457, -0.0442, -0.0351, 0.1685, 0.2023, 0.2131, 0.2229, 0.2637),
                         .Dim = 5:6, .Dimnames = list(c("FE", "ML", "REML", "DL", "HE"), c("tau2", "mu", "se", "z", "ci.lb", "ci.ub")))

   ### compare with results on page 309 (Table 16.3)
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

})

test_that("results are correct for the random-effects model (Huber-White method).", {

   res.std <- list()

   res.std$FE   <- rma(yi, vi, data=dat, digits=3, method="FE")
   res.std$ML   <- rma(yi, vi, data=dat, digits=3, method="ML")
   res.std$REML <- rma(yi, vi, data=dat, digits=3, method="REML")
   res.std$DL   <- rma(yi, vi, data=dat, digits=3, method="DL")
   res.std$HE   <- rma(yi, vi, data=dat, digits=3, method="HE")

   res.hw <- list()

   res.hw$FE   <- robust(res.std$FE,   cluster=dat$study, adjust=FALSE)
   res.hw$ML   <- robust(res.std$ML,   cluster=dat$study, adjust=FALSE)
   res.hw$REML <- robust(res.std$REML, cluster=dat$study, adjust=FALSE)
   res.hw$DL   <- robust(res.std$DL,   cluster=dat$study, adjust=FALSE)
   res.hw$HE   <- robust(res.std$HE,   cluster=dat$study, adjust=FALSE)

   out <- capture.output(print(res.hw$REML)) ### so that print.robust.rma() is run (at least once)

   tmp <- t(sapply(res.hw, function(x) c(tau2=x$tau2, mu=x$beta, se=x$se, t=x$zval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)))

   expected <- structure(c(0, 0.0126, 0.0188, 0.0259, 0.0804, 0.0604, 0.0777, 0.0837, 0.0893, 0.1143, 0.0398, 0.0475, 0.05, 0.0522, 0.0618, 1.5148, 1.6369, 1.6756, 1.7105, 1.8503, -0.0234, -0.022, -0.0213, -0.0204, -0.0155, 0.1441, 0.1775, 0.1887, 0.199, 0.2441),
                         .Dim = 5:6, .Dimnames = list(c("FE", "ML", "REML", "DL", "HE"), c("tau2", "mu", "se", "t", "ci.lb", "ci.ub")))

   ### compare with results on page 309 (Table 16.3)
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

})
