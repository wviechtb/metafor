### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:caterpillar_plot

source("settings.r")

context("Checking plots example: caterpillar plot")

test_that("plot can be drawn.", {

   skip_on_cran()

   ### simulate some data
   set.seed(5132)
   k <- 250
   vi <- rchisq(k, df=1) * .03
   yi <- rnorm(k, rnorm(k, 0.5, 0.4), sqrt(vi))

   ### fit RE model
   res <- rma(yi, vi)

   doplot <- function() {

      par(mar=c(5,1,2,1))

      forest(yi, vi, header=FALSE,
             xlim=c(-2.5,3.5), ylim=c(-8, 254),
             order=yi,
             slab=NA, annotate=FALSE,
             efac=0,
             pch=19,
             col="gray40",
             psize=2,
             cex.lab=1, cex.axis=1,
             lty=c("solid","blank"))

      points(sort(yi), k:1, pch=19, cex=0.5)

      addpoly(res, mlab="", cex=1)
      text(-2, -2, "RE Model", pos=4, offset=0, cex=1)

   }

   png("images/test_plots_caterpillar_plot_light_test.png", res=200, width=1800, height=1500, type="cairo")
   doplot()
   dev.off()

   expect_true(.vistest("images/test_plots_caterpillar_plot_light_test.png", "images/test_plots_caterpillar_plot_light.png"))

   png("images/test_plots_caterpillar_plot_dark_test.png", res=200, width=1800, height=1500, type="cairo")
   setmfopt(theme="dark")
   doplot()
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_caterpillar_plot_dark_test.png", "images/test_plots_caterpillar_plot_dark.png"))

})

rm(list=ls())
