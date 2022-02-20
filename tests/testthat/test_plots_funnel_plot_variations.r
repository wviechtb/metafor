### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:funnel_plot_variations

source("settings.r")

context("Checking plots example: funnel plot variations")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### fit equal-effects model
   res <- rma(yi, vi, data=dat.hackshaw1998, measure="OR", method="EE")

   ### set up 2x2 array for plotting
   par(mfrow=c(2,2))

   ### draw funnel plots
   funnel(res, main="Standard Error")
   funnel(res, yaxis="vi", main="Sampling Variance")
   funnel(res, yaxis="seinv", main="Inverse Standard Error")
   funnel(res, yaxis="vinv", main="Inverse Sampling Variance")

   par(opar)

})

rm(list=ls())
