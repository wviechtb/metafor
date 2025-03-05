### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:funnel_plot_with_trim_and_fill

source("settings.r")

context("Checking plots example: funnel plot with trim and fill")

test_that("plot can be drawn.", {

   skip_on_cran()

   res <- rma(yi, vi, data=dat.hackshaw1998, measure="OR")
   taf <- trimfill(res)
   out <- capture.output(print(taf))

   png("images/test_plots_funnel_plot_with_trim_and_fill_light_test.png", res=200, width=1800, height=1500, type="cairo")

   par(mar=c(5,4,1,2))
   funnel(taf, legend=list(show="cis"))

   dev.off()

   expect_true(.vistest("images/test_plots_funnel_plot_with_trim_and_fill_light_test.png", "images/test_plots_funnel_plot_with_trim_and_fill_light.png"))

   png("images/test_plots_funnel_plot_with_trim_and_fill_dark_test.png", res=200, width=1800, height=1500, type="cairo")

   setmfopt(theme="dark")

   par(mar=c(5,4,1,2))
   funnel(taf, legend=list(show="cis"))

   setmfopt(theme="default")

   dev.off()

   expect_true(.vistest("images/test_plots_funnel_plot_with_trim_and_fill_dark_test.png", "images/test_plots_funnel_plot_with_trim_and_fill_dark.png"))

})

rm(list=ls())
