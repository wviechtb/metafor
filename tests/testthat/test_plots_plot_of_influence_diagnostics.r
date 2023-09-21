### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:plot_of_influence_diagnostics

source("settings.r")

context("Checking plots example: plot of influence diagnostics")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### fit random-effects model with r-to-z transformed correlations
   res <- rma(ri=ri, ni=ni, measure="ZCOR", data=dat.mcdaniel1994)

   ### calculate influence diagnostics
   inf <- influence(res)

   ### plot the influence diagnostics
   png("test_plots_plot_of_influence_diagnostics_1.png", res=200, width=1800, height=3600, type="cairo")
   plot(inf, layout=c(8,1))
   dev.off()

   expect_true(.vistest("test_plots_plot_of_influence_diagnostics_1.png", "images/test_plots_plot_of_influence_diagnostics_1.png"))

   png("test_plots_plot_of_influence_diagnostics_2.png", res=200, width=1800, height=1800, type="cairo")
   plot(inf, plotinf=FALSE, plotdfbs=TRUE)
   dev.off()

   expect_true(.vistest("test_plots_plot_of_influence_diagnostics_2.png", "images/test_plots_plot_of_influence_diagnostics_2.png"))

   out <- capture.output(print(inf)) # so that print.infl.rma.uni() is run (at least once)

})

rm(list=ls())
