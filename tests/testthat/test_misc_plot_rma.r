### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: plot() function")

source("settings.r")

test_that("plot can be drawn for rma().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)

   png(filename="test_misc_plot_rma_1.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("test_misc_plot_rma_1.png", "images/test_misc_plot_rma_1.png"))

   res <- rma(yi ~ ablat, vi, data=dat)

   png(filename="test_misc_plot_rma_2.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("test_misc_plot_rma_2.png", "images/test_misc_plot_rma_2.png"))

})

test_that("plot can be drawn for rma.mh().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   png(filename="test_misc_plot_rma_3.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("test_misc_plot_rma_3.png", "images/test_misc_plot_rma_3.png"))

})

test_that("plot can be drawn for rma.peto().", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   png(filename="test_misc_plot_rma_4.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("test_misc_plot_rma_4.png", "images/test_misc_plot_rma_4.png"))

})

rm(list=ls())
