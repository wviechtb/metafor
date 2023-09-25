### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:baujat_plot

source("settings.r")

context("Checking plots example: Baujat plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### create Baujat plot
   png("images/test_plots_baujat_plot_test.png", res=200, width=1800, height=1800, type="cairo")

   ### adjust margins so the space is better used
   par(mar=c(5,4,2,2))

   ### load data from Pignon et al. (2000)
   dat <- dat.pignon2000

   ### compute estimated log hazard ratios and sampling variances
   dat$yi <- with(dat, OmE/V)
   dat$vi <- with(dat, 1/V)

   ### meta-analysis based on all 65 trials
   res <- rma(yi, vi, data=dat, method="EE", slab=id)

   baujat(res, xlim=c(0,20), ylim=c(0,0.2))

   dev.off()

   expect_true(.vistest("images/test_plots_baujat_plot_test.png", "images/test_plots_baujat_plot.png"))

})

rm(list=ls())
