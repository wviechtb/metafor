### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see also: https://www.metafor-project.org/doku.php/plots:labbe_plot

source("settings.r")

context("Checking plots example: L'Abbe plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR")

   png("images/test_plots_labbe_plot_1_light_test.png", res=200, width=1800, height=1600, type="cairo")
   par(mar=c(5,4,1,2))
   labbe(res, las=1, bty="l")
   dev.off()

   expect_true(.vistest("images/test_plots_labbe_plot_1_light_test.png", "images/test_plots_labbe_plot_1_light.png"))

   png("images/test_plots_labbe_plot_1_dark_test.png", res=200, width=1800, height=1600, type="cairo")
   setmfopt(theme="dark")
   par(mar=c(5,4,1,2))
   labbe(res, las=1, bty="l")
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_labbe_plot_1_dark_test.png", "images/test_plots_labbe_plot_1_dark.png"))

   png("images/test_plots_labbe_plot_2_light_test.png", res=200, width=1800, height=1600, type="cairo")
   par(mar=c(5,4,1,2))
   labbe(res, las=1, bty="l", ci=TRUE, grid=TRUE)
   dev.off()

   expect_true(.vistest("images/test_plots_labbe_plot_2_light_test.png", "images/test_plots_labbe_plot_2_light.png"))

   png("images/test_plots_labbe_plot_2_dark_test.png", res=200, width=1800, height=1600, type="cairo")
   setmfopt(theme="dark")
   par(mar=c(5,4,1,2))
   labbe(res, las=1, bty="l", ci=TRUE, grid=TRUE)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_labbe_plot_2_dark_test.png", "images/test_plots_labbe_plot_2_dark.png"))

})

rm(list=ls())
