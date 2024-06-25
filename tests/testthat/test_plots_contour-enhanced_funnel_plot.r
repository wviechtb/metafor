### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:contour_enhanced_funnel_plot

source("settings.r")

context("Checking plots example: contour-enhanced funnel plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
              slab=paste(author, year, sep=", "), method="REML")

   png("images/test_plots_contour_enhanced_funnel_plot_light_test.png", res=200, width=1800, height=1500, type="cairo")

   par(mar=c(5,4,1,2))

   funnel(res, level=c(90, 95, 99), refline=0, legend=TRUE)

   dev.off()

   expect_true(.vistest("images/test_plots_contour_enhanced_funnel_plot_light_test.png", "images/test_plots_contour_enhanced_funnel_plot_light.png"))

   png("images/test_plots_contour_enhanced_funnel_plot_dark_test.png", res=200, width=1800, height=1500, type="cairo")

   setmfopt(theme="dark")

   par(mar=c(5,4,1,2))

   funnel(res, level=c(90, 95, 99), refline=0, legend=TRUE)

   setmfopt(theme="default")

   dev.off()

   expect_true(.vistest("images/test_plots_contour_enhanced_funnel_plot_dark_test.png", "images/test_plots_contour_enhanced_funnel_plot_dark.png"))

})

rm(list=ls())
