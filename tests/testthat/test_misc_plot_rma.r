### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: plot() function")

source("settings.r")

test_that("plot can be drawn for rma().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)

   png(filename="images/test_misc_plot_rma_1_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_1_test.png", "images/test_misc_plot_rma_1.png"))

   res <- rma(yi ~ ablat, vi, data=dat)

   png(filename="images/test_misc_plot_rma_2_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_2_test.png", "images/test_misc_plot_rma_2.png"))

})

test_that("plot can be drawn for rma.mh().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   png(filename="images/test_misc_plot_rma_3_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_3_test.png", "images/test_misc_plot_rma_3.png"))

})

test_that("plot can be drawn for rma.peto().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   png(filename="images/test_misc_plot_rma_4_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_4_test.png", "images/test_misc_plot_rma_4.png"))

})

rm(list=ls())
