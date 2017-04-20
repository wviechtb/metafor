### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:plot_of_cumulative_results

context("Checking plots example: plot of cumulative results")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the more space is used
   par(mar=c(5,5,2,2))

   ### load BCG vaccine data
   data(dat.bcg, package="metafor")

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit random-effects models
   res <- rma(yi, vi, data=dat)

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=order(dat$year))

   ### plot of cumulative results
   plot(tmp, transf=exp, xlim=c(.25,.5), lwd=3, cex=1.3)

   par(opar)

})
