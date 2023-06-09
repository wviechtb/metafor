### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:baujat_plot

source("settings.r")

context("Checking plots example: scatter/bubble plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### calculate (log) risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit mixed-effects model with absolute latitude as predictor
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   ### draw plot
   sav <- regplot(res, xlim=c(10,60), predlim=c(10,60), xlab="Absolute Latitude", refline=0,
                  atransf=exp, at=log(seq(0.2,1.6,by=0.2)), digits=1, las=1, bty="l",
                  label=c(4,7,12,13), offset=c(1.6,0.8), labsize=0.9,
                  pi=TRUE, legend=TRUE, grid=TRUE)
   points(sav)

})

rm(list=ls())
