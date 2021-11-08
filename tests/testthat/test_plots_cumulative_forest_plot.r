### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:cumulative_forest_plot

context("Checking plots example: cumulative forest plot")

test_that("plot can be drawn for 'rma.uni' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(4,4,1,2))

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit random-effects models
   res <- rma(yi, vi, data=dat, slab=paste(author, year, sep=", "))

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=dat$year)

   ### cumulative forest plot
   forest(tmp, xlim=c(-4,2), at=log(c(.125, .25, .5, 1, 2)),
          atransf=exp, digits=c(2,3), cex=.75, header="Author(s) and Year")

   par(opar)

})

test_that("plot can be drawn for 'rma.mh' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(4,4,1,2))

   ### fit equal-effects models using the Mantel-Haenszel method
   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, slab=paste(author, year, sep=", "))

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=dat.bcg$year)

   ### cumulative forest plot
   forest(tmp, xlim=c(-4,2), at=log(c(.125, .25, .5, 1, 2)),
          atransf=exp, digits=c(2,3), cex=.75, header="Author(s) and Year")

   par(opar)

})

test_that("plot can be drawn for 'rma.peto' object.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(4,4,1,2))

   ### fit equal-effects models using Peto's method
   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, slab=paste(author, year, sep=", "))

   ### cumulative meta-analysis (in the order of publication year)
   tmp <- cumul(res, order=dat.bcg$year)

   ### cumulative forest plot
   forest(tmp, xlim=c(-4,2), at=log(c(.125, .25, .5, 1, 2)),
          atransf=exp, digits=c(2,3), cex=.75, header="Author(s) and Year")

   par(opar)

})
