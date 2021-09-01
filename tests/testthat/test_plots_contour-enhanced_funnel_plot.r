### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:contour_enhanced_funnel_plot

context("Checking plots example: contour-enhanced funnel plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(5,4,1,2))

   ### fit random-effects model
   res <- rma(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, measure="RR",
              slab=paste(author, year, sep=", "), method="REML")

   ### create contour enhanced funnel plot (with funnel centered at 0)
   funnel(res, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=TRUE)

   par(opar)

})
