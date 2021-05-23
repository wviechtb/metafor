### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:gosh_plot

context("Checking plots example: GOSH plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### load data
   data(dat.egger2001, package="metafor")

   ### meta-analysis of all trials including ISIS-4 using a FE model
   res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, method="FE")

   ### fit FE model to all possible subsets
   sav <- gosh(res, progbar=FALSE)
   out <- capture.output(print(sav)) ### so that print.gosh.rma() is run (at least once)

   ### create GOSH plot
   ### red points for subsets that include and blue points
   ### for subsets that exclude study 16 (the ISIS-4 trial)
   plot(sav, out=16, breaks=100)

   ### fit FE model to random subsets (with parallel processing)
   sav <- gosh(res, progbar=FALSE, parallel="snow", subsets=1000)

   ### meta-analysis using MH method (using subset to speed things up)
   res <- rma.mh(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, subset=c(1:7,16))
   sav <- gosh(res, progbar=FALSE)
   plot(sav, out=8, breaks=40)

   ### fit FE model to all possible subsets (with parallel processing)
   sav <- gosh(res, progbar=FALSE, parallel="snow", subsets=1000)

   ### meta-analysis using Peto's method (using subset to speed things up)
   res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, subset=c(1:7,16))
   sav <- gosh(res, progbar=FALSE)
   plot(sav, out=8, breaks=40)

   ### fit FE model to all possible subsets (with parallel processing)
   sav <- gosh(res, progbar=FALSE, parallel="snow", subsets=1000)

   par(opar)

})
