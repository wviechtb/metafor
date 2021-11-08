### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:normand1999

context("Checking analysis example: normand1999")

source("tolerances.r") # read in tolerances

test_that("results are correct for the first example (using dat.hine1989).", {

   ### calculate risk differences and corresponding sampling variances
   dat <- escalc(measure="RD", n1i=n1i, n2i=n2i, ai=ai, ci=ci, data=dat.hine1989)

   ### transform into percentage points
   dat$yi <- dat$yi * 100
   dat$vi <- dat$vi * 100^2

   out <- capture.output(print(dat)) ### so that print.escalc() is run (at least once)

   ### compare with results on page 330 (Table III)
   expect_equivalent(dat$yi, c(2.8026, 0.0000, 1.9711, 1.7961, 3.5334, 4.4031), tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, c(17.7575, 37.5657, 8.1323, 10.8998, 8.0114, 6.1320), tolerance=.tol[["var"]])

   ### CIs for individual studies
   tmp <- summary(dat)

   ### compare with results on page 330 (Table III)
   expect_equivalent(tmp$ci.lb, c(-5.4566, -12.0128, -3.6182, -4.6747, -2.0141, -0.4503), tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, c(11.0618, 12.0128, 7.5604, 8.2669, 9.0810, 9.2566), tolerance=.tol[["ci"]])

   ### fit equal-effects model
   res <- rma(yi, vi, data=dat, method="EE", digits=2)

   ### compare with results on page 349 (Table VII)
   expect_equivalent(coef(res), 2.9444, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.3831, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 5.5058, tolerance=.tol[["ci"]])

   ### fit random-effects model (REML estimator)
   res <- rma(yi, vi, data=dat, digits=2)

   ### compare with results on page 349 (Table VII)
   expect_equivalent(coef(res), 2.9444, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.3831, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 5.5058, tolerance=.tol[["ci"]])
   expect_equivalent(res$tau2,  0.0000, tolerance=.tol[["var"]])

   ### fit random-effects model (DL estimator)
   res <- rma(yi, vi, data=dat, method="DL", digits=2)

   ### compare with results on page 349 (Table VII)
   expect_equivalent(coef(res), 2.9444, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.3831, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 5.5058, tolerance=.tol[["ci"]])
   expect_equivalent(res$tau2,  0.0000, tolerance=.tol[["var"]])

})

test_that("results are correct for the second example (using dat.normand1999).", {

   ### compute pooled SD
   dat.normand1999$sdpi <- with(dat.normand1999, sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/(n1i+n2i-2)))

   ### calculate mean differences and corresponding sampling variances
   dat <- escalc(m1i=m1i, sd1i=sdpi, n1i=n1i, m2i=m2i, sd2i=sdpi, n2i=n2i, measure="MD", data=dat.normand1999, digits=2)

   ### compare with results on page 351 (Table VIII)
   expect_equivalent(dat$yi, c(-20, -2, -55, -71, -4, 1, 11, -10, 7))
   expect_equivalent(dat$vi, c(40.5863, 2.0468, 15.2809, 150.2222, 20.1923, 1.2235, 95.3756, 8.0321, 20.6936), tolerance=.tol[["var"]])

   ### CIs for individual studies
   tmp <- summary(dat)

   ### (results for this not given in paper)
   expect_equivalent(tmp$ci.lb, c(-32.4864, -4.8041, -62.6616, -95.0223, -12.8073, -1.168, -8.1411, -15.5547, -1.9159), tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, c(-7.5136, 0.8041, -47.3384, -46.9777, 4.8073, 3.168, 30.1411, -4.4453, 15.9159), tolerance=.tol[["ci"]])

   ### fit equal-effects model
   res <- rma(yi, vi, data=dat, method="EE", digits=2)

   ### compare with results on page 352 (Table IX)
   expect_equivalent(coef(res), -3.4939, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -5.0265, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, -1.9613, tolerance=.tol[["ci"]])

   ### fit random-effects model (DL estimator)
   res <- rma(yi, vi, data=dat, method="DL", digits=2)

   ### compare with results on page 352 (Table IX)
   expect_equivalent(coef(res), -14.0972, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -24.4454, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub,  -3.7490, tolerance=.tol[["ci"]])
   expect_equivalent(res$tau2,  218.7216, tolerance=.tol[["var"]])

   ### fit random-effects model (REML estimator)
   res <- rma(yi, vi, data=dat, digits=2)

   ### compare with results on page 352 (Table IX)
   expect_equivalent(coef(res), -15.1217, tolerance=.tol[["est"]])
   expect_equivalent(res$ci.lb, -32.6716, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub,   2.4282, tolerance=.tol[["ci"]])
   expect_equivalent(res$tau2,  685.1965, tolerance=.tol[["var"]])

})
