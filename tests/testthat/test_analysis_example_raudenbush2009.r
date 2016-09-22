### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:raudenbush2009

context("Checking analysis example: raudenbush2009")

### load data
dat <- get(data(dat.raudenbush1985, package="metafor"))

test_that("results are correct for the fixed-effects model.", {

   ### fixed-effects model
   res.FE <- rma(yi, vi, data=dat, digits=3, method="FE")

   ### compare with results on page 301 (Table 16.2) and page 302
   expect_equivalent(round(coef(res.FE),3), 0.060)
   expect_equivalent(round(res.FE$QE,2), 35.83)
   expect_equivalent(round(res.FE$zval,2), 1.66) ### 1.65 in chapter

})

test_that("results are correct for the random-effects model.", {

   ### random-effects model
   res.RE <- rma(yi, vi, data=dat, digits=3)

   ### compare with results on page 301 (Table 16.2) and page 302
   expect_equivalent(round(coef(res.RE),3), 0.084) ### 0.083 in chapter
   expect_equivalent(round(res.RE$zval,2), 1.62)
   expect_equivalent(round(res.RE$tau2,3), 0.019)

   ### credibility/prediction interval
   tmp <- predict(res.RE)

   ### compare with results on page 301 (Table 16.2) and page 302
   expect_equivalent(round(tmp$cr.lb,3), -0.204) ### -0.19 in chapter but computed in a slightly different way
   expect_equivalent(round(tmp$cr.ub,3),  0.371) ###  0.35 in chapter but computed in a slightly different way

   ### range of BLUPs
   tmp <- round(range(blup(res.RE)$pred), 2)

   ### compare with results on page 301 (Table 16.2)
   expect_equivalent(tmp, c(-0.03, 0.25))

})

test_that("results are correct for the mixed-effects model.", {

   ### recode weeks variable
   dat$weeks.c <- ifelse(dat$weeks > 3, 3, dat$weeks)

   ### mixed-effects model
   res.ME <- rma(yi, vi, mods = ~ weeks.c, data=dat, digits=3)

   ### compare with results on page 301 (Table 16.2)
   expect_equivalent(round(res.ME$tau2,3), 0)
   expect_equivalent(round(coef(res.ME),3), c(0.407, -0.157))
   expect_equivalent(round(res.ME$QE,2), 16.57)
   expect_equivalent(round(res.ME$zval,2), c(4.68, -4.39))

   ### range of BLUPs
   tmp <- round(range(blup(res.ME)$pred), 2)

   ### compare with results on page 301 (Table 16.2)
   expect_equivalent(tmp, c(-0.06, 0.41)) ### -0.07 in chapter

})

test_that("results are correct for the random-effects model (conventional approach).", {

   res.std <- list()

   res.std$FE   <- rma(yi, vi, data=dat, digits=3, method="FE")
   res.std$ML   <- rma(yi, vi, data=dat, digits=3, method="ML")
   res.std$REML <- rma(yi, vi, data=dat, digits=3, method="REML")
   res.std$DL   <- rma(yi, vi, data=dat, digits=3, method="DL")
   res.std$HE   <- rma(yi, vi, data=dat, digits=3, method="HE")

   tmp <- round(t(sapply(res.std, function(x) c(tau2=x$tau2, mu=x$b, se=x$se, z=x$zval, ci.lb=x$ci.lb, ci.ub=x$ci.ub))), 3)

   expected <- structure(c(0, 0.013, 0.019, 0.026, 0.08, 0.06, 0.078, 0.084, 0.089, 0.114, 0.036, 0.047, 0.052, 0.056, 0.079, 1.655, 1.637, 1.621, 1.601, 1.443, -0.011, -0.015, -0.018, -0.02, -0.041, 0.132, 0.171, 0.185, 0.199, 0.27),
                         .Dim = 5:6, .Dimnames = list(c("FE", "ML", "REML", "DL", "HE"), c("tau2", "mu", "se", "z", "ci.lb", "ci.ub")))

   ### compare with results on page 309 (Table 16.3)
   expect_equivalent(tmp, expected)

})

test_that("results are correct for the random-effects model (Knapp & Hartung method).", {

   res.knha <- list()

   expect_warning(res.knha$FE   <- rma(yi, vi, data=dat, digits=3, method="FE", test="knha"))
   res.knha$ML   <- rma(yi, vi, data=dat, digits=3, method="ML", test="knha")
   res.knha$REML <- rma(yi, vi, data=dat, digits=3, method="REML", test="knha")
   res.knha$DL   <- rma(yi, vi, data=dat, digits=3, method="DL", test="knha")
   res.knha$HE   <- rma(yi, vi, data=dat, digits=3, method="HE", test="knha")

   tmp <- round(t(sapply(res.knha, function(x) c(tau2=x$tau2, mu=x$b, se=x$se, z=x$zval, ci.lb=x$ci.lb, ci.ub=x$ci.ub))), 3)

   expected <- structure(c(0, 0.013, 0.019, 0.026, 0.08, 0.06, 0.078, 0.084, 0.089, 0.114, 0.051, 0.059, 0.062, 0.064, 0.071, 1.173, 1.311, 1.359, 1.405, 1.608, -0.048, -0.047, -0.046, -0.044, -0.035, 0.168, 0.202, 0.213, 0.223, 0.264),
                         .Dim = 5:6, .Dimnames = list(c("FE", "ML", "REML", "DL", "HE"), c("tau2", "mu", "se", "z", "ci.lb", "ci.ub")))

   ### compare with results on page 309 (Table 16.3)
   expect_equivalent(tmp, expected)

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

   tmp <- round(t(sapply(res.hw, function(x) c(tau2=x$tau2, mu=x$b, se=x$se, t=x$tval, ci.lb=x$ci.lb, ci.ub=x$ci.ub))), 3)

   expected <- structure(c(0, 0.013, 0.019, 0.026, 0.08, 0.06, 0.078, 0.084, 0.089, 0.114, 0.04, 0.047, 0.05, 0.052, 0.062, 1.515, 1.637, 1.676, 1.71, 1.85, -0.023, -0.022, -0.021, -0.02, -0.015, 0.144, 0.178, 0.189, 0.199, 0.244),
                         .Dim = 5:6, .Dimnames = list(c("FE", "ML", "REML", "DL", "HE"), c("tau2", "mu", "se", "t", "ci.lb", "ci.ub")))

   ### compare with results on page 309 (Table 16.3)
   expect_equivalent(tmp, expected)

})
