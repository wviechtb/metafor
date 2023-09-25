### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:radial_plot

source("settings.r")

context("Checking plots example: radial (Galbraith) plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_radial_plot_test.png", res=200, width=1800, height=1800, type="cairo")

   ### adjust margins so the space is better used
   par(mar=c(5,4,0,3))

   ### fit equal-effects model
   res <- rma(yi, vi, data=dat.hackshaw1998, method="EE")

   ### draw radial plot
   radial(res)

   dev.off()

   expect_true(.vistest("images/test_plots_radial_plot_test.png", "images/test_plots_radial_plot.png"))

})

rm(list=ls())
