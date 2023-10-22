### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:labbe_plot

source("settings.r")

context("Checking plots example: L'Abbe plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_labbe_plot_test.png", res=200, width=1800, height=1600, type="cairo")

   ### decrease margins so the full space is used
   par(mar=c(5,4,1,2))

   ### fit random-effects model
   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR")

   ### draw L'Abbé plot
   labbe(res, las=1, bty="l")

   dev.off()

   expect_true(.vistest("images/test_plots_labbe_plot_test.png", "images/test_plots_labbe_plot.png"))

})

rm(list=ls())
