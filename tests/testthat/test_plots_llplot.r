### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

source("settings.r")

context("Checking plots example: likelihood plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("test_plots_llplot.png", res=200, width=1800, height=1600, type="cairo")

   ### adjust margins so the space is better used
   par(mar=c(5,4,2,2))

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### create likelihood plot
   llplot(measure="GEN", yi=yi, vi=vi, data=dat, lwd=1, refline=NA, xlim=c(-3,2))

   dev.off()

   expect_true(.vistest("test_plots_llplot.png", "images/test_plots_llplot.png"))

})

rm(list=ls())
