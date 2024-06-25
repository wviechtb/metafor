### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/plots:cumulative_forest_plot

source("settings.r")

context("Checking plots example: cumulative forest plot")

test_that("plot can be drawn for 'rma.uni' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_cumulative_forest_plot_1_test.png", res=240, width=1800, height=1400, type="cairo")

   ### decrease margins so the full space is used
   par(mar=c(4,4,1,2))

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit random-effects models
   res <- rma(yi, vi, data=dat, slab=paste(author, year, sep=", "))

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=year)

   ### cumulative forest plot
   forest(tmp, xlim=c(-4,2), at=log(c(0.125, 0.25, 0.5, 1, 2)),
          atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")

   dev.off()

   expect_true(.vistest("images/test_plots_cumulative_forest_plot_1_test.png", "images/test_plots_cumulative_forest_plot_1.png"))


})

test_that("plot can be drawn for 'rma.mh' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_cumulative_forest_plot_2_test.png", res=240, width=1800, height=1400, type="cairo")

   ### decrease margins so the full space is used
   par(mar=c(4,4,1,2))

   ### fit equal-effects models using the Mantel-Haenszel method
   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, slab=paste(author, year, sep=", "))

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=dat.bcg$year)

   ### cumulative forest plot
   forest(tmp, xlim=c(-4,2), at=log(c(0.125, 0.25, 0.5, 1, 2)),
          atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")

   dev.off()

   expect_true(.vistest("images/test_plots_cumulative_forest_plot_2_test.png", "images/test_plots_cumulative_forest_plot_2.png"))

})

test_that("plot can be drawn for 'rma.peto' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   png("images/test_plots_cumulative_forest_plot_3_test.png", res=240, width=1800, height=1400, type="cairo")

   ### decrease margins so the full space is used
   par(mar=c(4,4,1,2))

   ### fit equal-effects models using Peto's method
   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, slab=paste(author, year, sep=", "))

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=dat.bcg$year)

   ### cumulative forest plot
   forest(tmp, xlim=c(-4,2), at=log(c(0.125, 0.25, 0.5, 1, 2)),
          atransf=exp, digits=c(2L,3L), cex=0.85, header="Author(s) and Year")

   dev.off()

   expect_true(.vistest("images/test_plots_cumulative_forest_plot_3_test.png", "images/test_plots_cumulative_forest_plot_3.png"))

})

rm(list=ls())
