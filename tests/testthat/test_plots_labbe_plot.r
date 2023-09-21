### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:labbe_plot

source("settings.r")

context("Checking plots example: L'Abbe plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("test_plots_labbe_plot.png", res=200, width=1800, height=1600, type="cairo")

   ### decrease margins so the full space is used
   par(mar=c(5,4,1,2))

   ### fit random-effects model
   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR")

   ### draw L'Abbé plot
   labbe(res)

   dev.off()

   expect_true(.vistest("test_plots_labbe_plot.png", "images/test_plots_labbe_plot.png"))

})

rm(list=ls())
