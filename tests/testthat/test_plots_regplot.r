### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:baujat_plot

source("settings.r")

context("Checking plots example: scatter/bubble plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   png("images/test_plots_regplot_light_test.png", res=200, width=1800, height=1500, type="cairo")

   par(mar=c(5,5,1,2))

   sav <- regplot(res, xlim=c(10,60), predlim=c(10,60), xlab="Absolute Latitude", refline=0,
                  atransf=exp, at=log(seq(0.2,1.6,by=0.2)), digits=1, las=1, bty="l",
                  label=c(4,7,12,13), offset=c(1.6,0.8), labsize=0.9,
                  pi=TRUE, legend=TRUE, grid=TRUE)
   points(sav)

   dev.off()

   expect_true(.vistest("images/test_plots_regplot_light_test.png", "images/test_plots_regplot_light.png"))

   png("images/test_plots_regplot_dark_test.png", res=200, width=1800, height=1500, type="cairo")

   setmfopt(theme="dark")

   par(mar=c(5,5,1,2))

   sav <- regplot(res, xlim=c(10,60), predlim=c(10,60), xlab="Absolute Latitude", refline=0,
                  atransf=exp, at=log(seq(0.2,1.6,by=0.2)), digits=1, las=1, bty="l",
                  label=c(4,7,12,13), offset=c(1.6,0.8), labsize=0.9,
                  pi=TRUE, legend=TRUE, grid=TRUE)
   points(sav)

   setmfopt(theme="default")

   dev.off()

   expect_true(.vistest("images/test_plots_regplot_dark_test.png", "images/test_plots_regplot_dark.png"))

})

rm(list=ls())
