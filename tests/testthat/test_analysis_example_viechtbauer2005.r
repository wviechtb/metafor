### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:viechtbauer2005

context("Checking analysis example: viechtbauer2005")

source("settings.r")

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
                     tau2=sapply(res, function(x) x$tau2),
                     I2=sapply(res, function(x) x$I2),
                     H2=sapply(res, function(x) x$H2),
                     se.tau2=sapply(res, function(x) x$se.tau2))

   ### compare with results on page 271
   expect_equivalent(res$tau2,    c(0.2282, 0.1484, 0.2768, 0.1967, 0.2232, 0.192, 0.1992), tolerance=.tol[["var"]])
   expect_equivalent(res$I2,      c(77.2284, 68.7988, 80.4447, 74.5098, 76.8399, 74.0511, 74.7545), tolerance=.tol[["het"]])
   expect_equivalent(res$H2,      c(4.3914, 3.205, 5.1137, 3.9231, 4.3178, 3.8537, 3.9611), tolerance=.tol[["het"]])
   expect_equivalent(res$se.tau2, c(0.1328, 0.1234, 0.1841, 0.1255, 0.1464, 0.133, 0.0979), tolerance=.tol[["sevar"]])

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
                     tau2=sapply(res, function(x) x$tau2),
                     I2=sapply(res, function(x) x$I2),
                     H2=sapply(res, function(x) x$H2),
                     se.tau2=sapply(res, function(x) x$se.tau2))

   ### compare with results on page 272
   expect_equivalent(res$tau2,    c(0.0099, 0, 0.0126, 0.0132, 0.0157, 0.0104, 0.0248), tolerance=.tol[["var"]])
   expect_equivalent(res$I2,      c(22.9266, 0, 27.5275, 28.4505, 32.0203, 23.7198, 42.6734), tolerance=.tol[["het"]])
   expect_equivalent(res$H2,      c(1.2975, 1, 1.3798, 1.3976, 1.471, 1.311, 1.7444), tolerance=.tol[["het"]])
   expect_equivalent(res$se.tau2, c(0.0138, 0.0217, 0.0159, 0.0151, 0.0167, 0.0156, 0.0118), tolerance=.tol[["sevar"]])

})

rm(list=ls())
