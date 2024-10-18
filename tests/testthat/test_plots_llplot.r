### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

source("settings.r")

context("Checking plots example: likelihood plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   png("images/test_plots_llplot_light_test.png", res=200, width=1800, height=1600, type="cairo")

   par(mar=c(5,4,2,2))
   llplot(measure="GEN", yi=yi, vi=vi, data=dat, lwd=1, refline=NA, xlim=c(-3,2))

   dev.off()

   expect_true(.vistest("images/test_plots_llplot_light_test.png", "images/test_plots_llplot_light.png"))

   png("images/test_plots_llplot_dark_test.png", res=200, width=1800, height=1600, type="cairo")

   setmfopt(theme="dark")

   par(mar=c(5,4,2,2))
   llplot(measure="GEN", yi=yi, vi=vi, data=dat, lwd=1, refline=NA, xlim=c(-3,2))

   setmfopt(theme="default")

   dev.off()

   expect_true(.vistest("images/test_plots_llplot_dark_test.png", "images/test_plots_llplot_dark.png"))

})

rm(list=ls())
