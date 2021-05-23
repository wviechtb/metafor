### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:baujat_plot

context("Checking plots example: Likelihood plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### adjust margins so the space is better used
   par(mar=c(5,4,2,2))

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### create likelihood plot
   llplot(measure="GEN", yi=yi, vi=vi, data=dat, lwd=1, refline=NA, xlim=c(-3,2))

   par(opar)

})
