### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:labbe_plot

context("Checking plots example: L'Abbé plot")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(5,4,1,2))

   ### load BCG vaccine data
   data(dat.bcg)

   ### fit random-effects model
   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR")

   ### draw L'Abbé plot
   labbe(res)

   par(opar)

})
