### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:normal_qq_plots

source("settings.r")

context("Checking plots example: normal QQ plots")

test_that("plot can be drawn for 'rma.uni' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_normal_qq_plots_1_test.png", res=200, width=1800, height=1800, type="cairo")

   ### set up 2x2 array for plotting
   par(mfrow=c(2,2))

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit equal- and random-effects models
   res1 <- rma(yi, vi, data=dat, method="EE")
   res2 <- rma(yi, vi, data=dat)

   ### fit fixed- and random-effects models with absolute latitude moderator
   res3 <- rma(yi, vi, mods=~ablat, data=dat, method="FE")
   res4 <- rma(yi, vi, mods=~ablat, data=dat)

   ### normal QQ plots for the various models
   qqnorm(res1, seed=1234, grid=TRUE, main="Equal-Effects Model")
   qqnorm(res2, seed=1234, grid=TRUE, main="Random-Effects Model")
   qqnorm(res3, seed=1234, grid=TRUE, main="Fixed-Effects with Moderators Model")
   qqnorm(res4, seed=1234, grid=TRUE, main="Mixed-Effects Model")

   dev.off()

   expect_true(.vistest("images/test_plots_normal_qq_plots_1_test.png", "images/test_plots_normal_qq_plots_1.png"))

   ### draw plot with studentized residuals and labels
   png("images/test_plots_normal_qq_plots_2_test.png", res=200, width=1800, height=1800, type="cairo")
   qqnorm(res2, type="rstudent", grid=TRUE, label=TRUE, seed=1234)
   dev.off()

   expect_true(.vistest("images/test_plots_normal_qq_plots_2_test.png", "images/test_plots_normal_qq_plots_2.png"))

})

test_that("plot can be drawn for 'rma.mh' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_normal_qq_plots_3_test.png", res=200, width=1800, height=1800, type="cairo")

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   qqnorm(res)
   qqnorm(res, type="rstudent", label=TRUE)

   dev.off()

   expect_true(.vistest("images/test_plots_normal_qq_plots_3_test.png", "images/test_plots_normal_qq_plots_3.png"))

})

test_that("plot can be drawn for 'rma.peto' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_normal_qq_plots_4_test.png", res=200, width=1800, height=1800, type="cairo")

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   qqnorm(res)
   qqnorm(res, type="rstudent", label=TRUE)

   dev.off()

   expect_true(.vistest("images/test_plots_normal_qq_plots_4_test.png", "images/test_plots_normal_qq_plots_4.png"))

})

test_that("plot cannot be drawn for 'rma.mv' object.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, sparse=.sparse)
   expect_error(qqnorm(res))

})

rm(list=ls())
