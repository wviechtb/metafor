### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:gosh_plot

source("settings.r")

context("Checking plots example: GOSH plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### meta-analysis of all trials including ISIS-4 using an equal-effects model
   res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, method="EE")

   ### fit EE model to all possible subsets
   sav <- gosh(res, progbar=FALSE)
   out <- capture.output(print(sav)) # so that print.gosh.rma() is run (at least once)

   ### create GOSH plot
   ### red points for subsets that include and blue points
   ### for subsets that exclude study 16 (the ISIS-4 trial)
   png("images/test_plots_gosh_1_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(sav, out=16, breaks=100)
   dev.off()

   expect_true(.vistest("images/test_plots_gosh_1_light_test.png", "images/test_plots_gosh_1_light.png"))

   png("images/test_plots_gosh_1_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(sav, out=16, breaks=100)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_gosh_1_dark_test.png", "images/test_plots_gosh_1_dark.png"))

   ### fit EE model to random subsets (with parallel processing)
   sav <- gosh(res, progbar=FALSE, parallel="snow", subsets=1000)

   ### meta-analysis using MH method (using subset to speed things up)
   res <- rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, subset=c(1:7,16))
   sav <- gosh(res, progbar=FALSE)

   ### create GOSH plot
   png("images/test_plots_gosh_2_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(sav, out=8, breaks=40)
   dev.off()

   expect_true(.vistest("images/test_plots_gosh_2_light_test.png", "images/test_plots_gosh_2_light.png"))

   png("images/test_plots_gosh_2_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(sav, out=8, breaks=40)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_gosh_2_dark_test.png", "images/test_plots_gosh_2_dark.png"))

   ### fit EE model to all possible subsets (with parallel processing)
   sav <- gosh(res, progbar=FALSE, parallel="snow", subsets=1000)

   ### meta-analysis using Peto's method (using subset to speed things up)
   res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, subset=c(1:7,16))
   sav <- gosh(res, progbar=FALSE)

   ### create GOSH plot
   png("images/test_plots_gosh_3_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(sav, out=8, breaks=40)
   dev.off()

   expect_true(.vistest("images/test_plots_gosh_3_light_test.png", "images/test_plots_gosh_3_light.png"))

   png("images/test_plots_gosh_3_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(sav, out=8, breaks=40)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_plots_gosh_3_dark_test.png", "images/test_plots_gosh_3_dark.png"))

   ### fit EE model to all possible subsets (with parallel processing)
   sav <- gosh(res, progbar=FALSE, parallel="snow", subsets=1000)

})

rm(list=ls())
