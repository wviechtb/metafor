### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: dfround() function")

source("settings.r")

test_that("dfround() works correctly.", {

   dat <- as.data.frame(dat.raudenbush1985)
   dat$yi <- c(dat$yi)
   dat <- dfround(dat, c(rep(NA,8), 2, 3))
   expect_identical(dat$yi, c(0.03, 0.12, -0.14, 1.18, 0.26, -0.06, -0.02, -0.32, 0.27, 0.8, 0.54, 0.18, -0.02, 0.23, -0.18, -0.06, 0.3, 0.07, -0.07))
   expect_identical(dat$vi, c(0.016, 0.022, 0.028, 0.139, 0.136, 0.011, 0.011, 0.048, 0.027, 0.063, 0.091, 0.05, 0.084, 0.084, 0.025, 0.028, 0.019, 0.009, 0.03))

})

rm(list=ls())
