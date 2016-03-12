### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:gosh_plot

context("Checking plots example: GOSH plot")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### load data
   data(dat.egger2001, package="metafor")

   ### meta-analysis of all trials including ISIS-4 using a FE model
   res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.egger2001, method="FE")

   ### fit FE model to all possible subsets
   sav <- gosh(res, progbar=FALSE)

   ### create GOSH plot
   ### red points for subsets that include and blue points
   ### for subsets that exclude study 16 (the ISIS-4 trial)
   plot(sav, out=16, breaks=100, adjust=.5)

   par(opar)

})
