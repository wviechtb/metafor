### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: plot() function")

source("settings.r")

test_that("plot can be drawn for rma().", {

   skip_on_cran()

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)

   png(filename="images/test_misc_plot_rma_1_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_1_light_test.png", "images/test_misc_plot_rma_1_light.png"))

   png(filename="images/test_misc_plot_rma_1_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(res)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_1_dark_test.png", "images/test_misc_plot_rma_1_dark.png"))

   res <- rma(yi ~ ablat, vi, data=dat)

   png(filename="images/test_misc_plot_rma_2_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_2_light_test.png", "images/test_misc_plot_rma_2_light.png"))

   png(filename="images/test_misc_plot_rma_2_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(res)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_2_dark_test.png", "images/test_misc_plot_rma_2_dark.png"))

})

test_that("plot can be drawn for rma.mh().", {

   skip_on_cran()

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   png(filename="images/test_misc_plot_rma_3_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_3_light_test.png", "images/test_misc_plot_rma_3_light.png"))

   png(filename="images/test_misc_plot_rma_3_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(res)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_3_dark_test.png", "images/test_misc_plot_rma_3_dark.png"))

})

test_that("plot can be drawn for rma.peto().", {

   skip_on_cran()

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   png(filename="images/test_misc_plot_rma_4_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(res)
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_4_light_test.png", "images/test_misc_plot_rma_4_light.png"))

   png(filename="images/test_misc_plot_rma_4_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(res)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_misc_plot_rma_4_dark_test.png", "images/test_misc_plot_rma_4_dark.png"))

})

rm(list=ls())
