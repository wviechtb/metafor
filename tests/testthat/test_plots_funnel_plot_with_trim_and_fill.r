### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:funnel_plot_with_trim_and_fill

context("Checking plots example: funnel plot with trim and fill")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(5,4,1,2))

   ### load ETS data
   data(dat.hackshaw1998, package="metafor")

   ### fit random-effects model
   res <- rma(yi, vi, data=dat.hackshaw1998, measure="OR")

   ### carry out trim-and-fill analysis
   taf <- trimfill(res)

   ### draw funnel plot with missing studies filled in
   funnel(taf)

   par(opar)

})
