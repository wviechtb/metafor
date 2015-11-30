### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:viechtbauer2005

context("Checking analysis example viechtbauer2005")

### create dataset for example 1
dat <- data.frame(
id=1:10,
yi = c(-0.581, 0.530, 0.771, 1.031, 0.553, 0.295, 0.078, 0.573, -0.176, -0.232),
vi = c(0.023, 0.052, 0.060, 0.115, 0.095, 0.203, 0.200, 0.211, 0.051, 0.040))

test_that("results are correct for example 1.", {

   res.HS   <- rma(yi, vi, data=dat, method="HS")
   res.HE   <- rma(yi, vi, data=dat, method="HE")
   res.DL   <- rma(yi, vi, data=dat, method="DL")
   res.ML   <- rma(yi, vi, data=dat, method="ML")
   res.REML <- rma(yi, vi, data=dat, method="REML")
   res.EB   <- rma(yi, vi, data=dat, method="EB")
   res.SJ   <- rma(yi, vi, data=dat, method="SJ")

   res <- list(res.HS, res.HE, res.DL, res.ML, res.REML, res.EB, res.SJ)
   res <- data.frame(method=sapply(res, function(x) x$method),
                     tau2=sapply(res, function(x) round(x$tau2,3)),
                     I2=sapply(res, function(x) round(x$I2,2)),
                     H2=sapply(res, function(x) round(x$H2,2)))

   ### compare with results on page 271
   expect_equivalent(res$tau2, c(0.228, 0.148, 0.277, 0.197, 0.223, 0.192, 0.199))
   expect_equivalent(res$I2,   c(77.23, 68.80, 80.44, 74.51, 76.84, 74.05, 74.75))
   expect_equivalent(res$H2,   c(4.39, 3.21, 5.11, 3.92, 4.32, 3.85, 3.96))

})

### create dataset for example 2
dat <- data.frame(
id=1:18,
yi = c(0.100, -0.162, -0.090, -0.049, -0.046, -0.010, -0.431, -0.261, 0.134, 0.019, 0.175, 0.056, 0.045, 0.103, 0.121, -0.482, 0.290, 0.342),
vi = c(0.016, 0.015, 0.050, 0.050, 0.032, 0.052, 0.036, 0.024, 0.034, 0.033, 0.031, 0.034, 0.039, 0.167, 0.134, 0.096, 0.016, 0.035))

test_that("results are correct for example 2.", {

   res.HS   <- rma(yi, vi, data=dat, method="HS")
   res.HE   <- rma(yi, vi, data=dat, method="HE")
   res.DL   <- rma(yi, vi, data=dat, method="DL")
   res.ML   <- rma(yi, vi, data=dat, method="ML")
   res.REML <- rma(yi, vi, data=dat, method="REML")
   res.EB   <- rma(yi, vi, data=dat, method="EB")
   res.SJ   <- rma(yi, vi, data=dat, method="SJ")

   res <- list(res.HS, res.HE, res.DL, res.ML, res.REML, res.EB, res.SJ)

   res <- data.frame(method=sapply(res, function(x) x$method),
                     tau2=sapply(res, function(x) round(x$tau2,3)),
                     I2=sapply(res, function(x) round(x$I2,2)),
                     H2=sapply(res, function(x) round(x$H2,2)))

   ### compare with results on page 272
   expect_equivalent(res$tau2, c(0.010, 0.000, 0.013, 0.013, 0.016, 0.010, 0.025))
   expect_equivalent(res$I2,   c(22.93, 0.00, 27.53, 28.45, 32.02, 23.72, 42.67))
   expect_equivalent(res$H2,   c(1.30, 1.00, 1.38, 1.40, 1.47, 1.31, 1.74))

})
