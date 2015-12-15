### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:normand1999

context("Checking analysis example normand1999")

test_that("results are correct for the first example (using dat.hine1989).", {

   ### load data
   data(dat.hine1989, package="metafor")

   ### calculate risk differences and corresponding sampling variances
   dat <- escalc(measure="RD", n1i=n1i, n2i=n2i, ai=ai, ci=ci, data=dat.hine1989)

   ### transform into percentage points
   dat$yi <- dat$yi * 100
   dat$vi <- dat$vi * 100^2

   print(dat) ### so that print.escalc() is run (at least once)

   ### compare with results on page 330 (Table III)
   expect_equivalent(round(dat$yi,4), c(2.8026, 0.0000, 1.9711, 1.7961, 3.5334, 4.4031))
   expect_equivalent(round(dat$vi,4), c(17.7575, 37.5657, 8.1323, 10.8998, 8.0114, 6.1320))

   ### CIs for individual studies
   tmp <- summary(dat)

   ### compare with results on page 330 (Table III)
   expect_equivalent(round(tmp$ci.lb,1), c(-5.5, -12.0, -3.6, -4.7, -2.0, -0.5))
   expect_equivalent(round(tmp$ci.ub,1), c(11.1, 12.0, 7.6, 8.3, 9.1, 9.3))

   ### fit fixed-effects model
   res <- rma(yi, vi, data=dat, method="FE", digits=2)

   ### compare with results on page 349 (Table VII)
   expect_equivalent(round(coef(res),2), 2.94)
   expect_equivalent(round(res$ci.lb,1), 0.4)
   expect_equivalent(round(res$ci.ub,1), 5.5)

   ### fit random-effects model (REML estimator)
   res <- rma(yi, vi, data=dat, digits=2)

   ### compare with results on page 349 (Table VII)
   expect_equivalent(round(coef(res),2), 2.94)
   expect_equivalent(round(res$ci.lb,1), 0.4)
   expect_equivalent(round(res$ci.ub,1), 5.5)
   expect_equivalent(round(res$tau2,3), 0.000)

   ### fit random-effects model (DL estimator)
   res <- rma(yi, vi, data=dat, method="DL", digits=2)

   ### compare with results on page 349 (Table VII)
   expect_equivalent(round(coef(res),2), 2.94)
   expect_equivalent(round(res$ci.lb,1), 0.4)
   expect_equivalent(round(res$ci.ub,1), 5.5)
   expect_equivalent(round(res$tau2,3), 0.000)

})

test_that("results are correct for the second example (using dat.normand1999).", {

   ### load data
   data(dat.normand1999, package="metafor")

   ### compute pooled SD
   dat.normand1999$sdpi <- with(dat.normand1999, sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/(n1i+n2i-2)))

   ### calculate mean differences and corresponding sampling variances
   dat <- escalc(m1i=m1i, sd1i=sdpi, n1i=n1i, m2i=m2i, sd2i=sdpi, n2i=n2i, measure="MD", data=dat.normand1999, digits=2)

   ### compare with results on page 351 (Table VIII)
   expect_equivalent(dat$yi, c(-20, -2, -55, -71, -4, 1, 11, -10, 7))
   expect_equivalent(round(dat$vi,2), c(40.59, 2.05, 15.28, 150.22, 20.19, 1.22, 95.38, 8.03, 20.69))

   ### CIs for individual studies
   tmp <- summary(dat)

   ### (results for this not given in paper)
   expect_equivalent(round(tmp$ci.lb,1), c(-32.5, -4.8, -62.7, -95.0, -12.8, -1.2, -8.1, -15.6, -1.9))
   expect_equivalent(round(tmp$ci.ub,1), c(-7.5, 0.8, -47.3, -47.0, 4.8, 3.2, 30.1, -4.4, 15.9))

   ### fit fixed-effects model
   res <- rma(yi, vi, data=dat, method="FE", digits=2)

   ### compare with results on page 352 (Table IX)
   expect_equivalent(round(coef(res),2), -3.49)
   expect_equivalent(round(res$ci.lb,2), -5.03)
   expect_equivalent(round(res$ci.ub,2), -1.96)

   ### fit random-effects model (DL estimator)
   res <- rma(yi, vi, data=dat, method="DL", digits=2)

   ### compare with results on page 352 (Table IX)
   expect_equivalent(round(coef(res),2), -14.10)
   expect_equivalent(round(res$ci.lb,2), -24.45)
   expect_equivalent(round(res$ci.ub,2),  -3.75)
   expect_equivalent(round(res$tau2,2), 218.72)

   ### fit random-effects model (REML estimator)
   res <- rma(yi, vi, data=dat, digits=2)

   ### compare with results on page 352 (Table IX)
   expect_equivalent(round(coef(res),2), -15.12)
   expect_equivalent(round(res$ci.lb,2), -32.67)
   expect_equivalent(round(res$ci.ub,2),   2.43)
   expect_equivalent(round(res$tau2,2), 685.20)

})
