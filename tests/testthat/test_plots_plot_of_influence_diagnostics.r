### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:plot_of_influence_diagnostics

source("settings.r")

context("Checking plots example: plot of influence diagnostics")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma(ri=ri, ni=ni, measure="ZCOR", data=dat.mcdaniel1994)
   inf <- influence(res)
   out <- capture.output(print(inf)) # so that print.infl.rma.uni() is run (at least once)

   png("images/test_plots_plot_of_influence_diagnostics_1_light_test.png", res=200, width=1800, height=3600, type="cairo")
   plot(inf, layout=c(8,1))
   dev.off()

   expect_true(.vistest("images/test_plots_plot_of_influence_diagnostics_1_light_test.png", "images/test_plots_plot_of_influence_diagnostics_1_light.png"))

   png("images/test_plots_plot_of_influence_diagnostics_1_dark_test.png", res=200, width=1800, height=3600, type="cairo")
   setmfopt(theme="dark")
   plot(inf, layout=c(8,1))
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_plot_of_influence_diagnostics_1_dark_test.png", "images/test_plots_plot_of_influence_diagnostics_1_dark.png"))

   png("images/test_plots_plot_of_influence_diagnostics_2_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(inf, plotinf=FALSE, plotdfbs=TRUE)
   dev.off()

   expect_true(.vistest("images/test_plots_plot_of_influence_diagnostics_2_light_test.png", "images/test_plots_plot_of_influence_diagnostics_2_light.png"))

   png("images/test_plots_plot_of_influence_diagnostics_2_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(inf, plotinf=FALSE, plotdfbs=TRUE)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_plot_of_influence_diagnostics_2_dark_test.png", "images/test_plots_plot_of_influence_diagnostics_2_dark.png"))

})

rm(list=ls())
