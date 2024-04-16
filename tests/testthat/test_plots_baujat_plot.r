### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see also: https://www.metafor-project.org/doku.php/plots:baujat_plot

source("settings.r")

context("Checking plots example: Baujat plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   doplot <- function() {

      par(mar=c(5,4,2,2))

      dat <- dat.pignon2000
      dat$yi <- with(dat, OmE/V)
      dat$vi <- with(dat, 1/V)

      res <- rma(yi, vi, data=dat, method="EE", slab=id)

      baujat(res, xlim=c(0,20), ylim=c(0,0.2), bty="l", las=1)

   }

   png("images/test_plots_baujat_plot_light_test.png", res=200, width=1800, height=1800, type="cairo")
   doplot()
   dev.off()

   expect_true(.vistest("images/test_plots_baujat_plot_light_test.png", "images/test_plots_baujat_plot_light.png"))

   png("images/test_plots_baujat_plot_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   doplot()
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_baujat_plot_dark_test.png", "images/test_plots_baujat_plot_dark.png"))

})

rm(list=ls())
