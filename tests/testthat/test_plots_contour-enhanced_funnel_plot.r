### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:contour_enhanced_funnel_plot

source("settings.r")

context("Checking plots example: contour-enhanced funnel plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_contour_enhanced_funnel_plot_test.png", res=200, width=1800, height=1500, type="cairo")

   ### decrease margins so the full space is used
   par(mar=c(5,4,1,2))

   ### fit random-effects model
   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
              slab=paste(author, year, sep=", "), method="REML")

   ### create contour enhanced funnel plot (with funnel centered at 0)
   funnel(res, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=TRUE)

   dev.off()

   expect_true(.vistest("images/test_plots_contour_enhanced_funnel_plot_test.png", "images/test_plots_contour_enhanced_funnel_plot.png"))

})

rm(list=ls())
