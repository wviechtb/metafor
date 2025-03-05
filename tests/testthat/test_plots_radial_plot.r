### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:radial_plot

source("settings.r")

context("Checking plots example: radial (Galbraith) plot")

test_that("plot can be drawn.", {

   skip_on_cran()

   res <- rma(yi, vi, data=dat.hackshaw1998, method="EE")

   png("images/test_plots_radial_plot_light_test.png", res=200, width=1800, height=1800, type="cairo")

   par(mar=c(5,4,0,3))
   radial(res)

   dev.off()

   expect_true(.vistest("images/test_plots_radial_plot_light_test.png", "images/test_plots_radial_plot_light.png"))

   png("images/test_plots_radial_plot_dark_test.png", res=200, width=1800, height=1800, type="cairo")

   setmfopt(theme="dark")

   par(mar=c(5,4,0,3))
   radial(res)

   setmfopt(theme="default")

   dev.off()

   expect_true(.vistest("images/test_plots_radial_plot_dark_test.png", "images/test_plots_radial_plot_dark.png"))

})

rm(list=ls())
