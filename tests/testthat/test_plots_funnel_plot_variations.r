### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:funnel_plot_variations

source("settings.r")

context("Checking plots example: funnel plot variations")

test_that("plot can be drawn.", {

   skip_on_cran()

   ### fit equal-effects model
   res <- rma(yi, vi, data=dat.hackshaw1998, measure="OR", method="EE")

   png("images/test_plots_funnel_plot_variations_light_test.png", res=200, width=1800, height=1800, type="cairo")

   par(mfrow=c(2,2))

   funnel(res, main="Standard Error")
   funnel(res, yaxis="vi", main="Sampling Variance")
   funnel(res, yaxis="seinv", main="Inverse Standard Error")
   funnel(res, yaxis="vinv", main="Inverse Sampling Variance")

   dev.off()

   expect_true(.vistest("images/test_plots_funnel_plot_variations_light_test.png", "images/test_plots_funnel_plot_variations_light.png"))

   png("images/test_plots_funnel_plot_variations_dark_test.png", res=200, width=1800, height=1800, type="cairo")

   setmfopt(theme="dark")

   par(mfrow=c(2,2))

   funnel(res, main="Standard Error")
   funnel(res, yaxis="vi", main="Sampling Variance")
   funnel(res, yaxis="seinv", main="Inverse Standard Error")
   funnel(res, yaxis="vinv", main="Inverse Sampling Variance")

   setmfopt(theme="default")

   dev.off()

   expect_true(.vistest("images/test_plots_funnel_plot_variations_dark_test.png", "images/test_plots_funnel_plot_variations_dark.png"))

})

rm(list=ls())
