### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see also: https://www.metafor-project.org/doku.php/plots:plot_of_cumulative_results

source("settings.r")

context("Checking plots example: plot of cumulative results")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   tmp <- cumul(res, order=year)

   png("images/test_plots_plot_of_cumulative_results_light_test.png", res=200, width=1800, height=1600, type="cairo")

   par(mar=c(5,5,2,2))
   plot(tmp, transf=exp, xlim=c(0.25,0.5), lwd=3, cex=1.3)

   dev.off()

   expect_true(.vistest("images/test_plots_plot_of_cumulative_results_light_test.png", "images/test_plots_plot_of_cumulative_results_light.png"))

   png("images/test_plots_plot_of_cumulative_results_dark_test.png", res=200, width=1800, height=1600, type="cairo")

   setmfopt(theme="dark")

   par(mar=c(5,5,2,2))
   plot(tmp, transf=exp, xlim=c(0.25,0.5), lwd=3, cex=1.3)

   setmfopt(theme="default")

   dev.off()

   expect_true(.vistest("images/test_plots_plot_of_cumulative_results_dark_test.png", "images/test_plots_plot_of_cumulative_results_dark.png"))

})

rm(list=ls())
