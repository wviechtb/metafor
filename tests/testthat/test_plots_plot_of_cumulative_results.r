### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:plot_of_cumulative_results

source("settings.r")

context("Checking plots example: plot of cumulative results")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_plot_of_cumulative_results_test.png", res=200, width=1800, height=1600, type="cairo")

   ### decrease margins so the more space is used
   par(mar=c(5,5,2,2))

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit random-effects models
   res <- rma(yi, vi, data=dat)

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=year)

   ### plot of cumulative results
   plot(tmp, transf=exp, xlim=c(0.25,0.5), lwd=3, cex=1.3)

   dev.off()

   expect_true(.vistest("images/test_plots_plot_of_cumulative_results_test.png", "images/test_plots_plot_of_cumulative_results.png"))

})

rm(list=ls())
