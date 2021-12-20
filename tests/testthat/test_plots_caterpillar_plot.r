### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:caterpillar_plot

context("Checking plots example: Caterpillar plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### simulate some data
   set.seed(5132)
   k <- 250
   vi <- rchisq(k, df=1) * .03
   yi <- rnorm(k, rnorm(k, 0.5, 0.4), sqrt(vi))

   ### fit RE model
   res <- rma(yi, vi)

   opar <- par(no.readonly=TRUE)

   ### decrease margins so the full space is used
   par(mar=c(5,1,1,1))

   ### create plot
   forest(yi, vi,
          xlim=c(-2.5,3.5),        ### adjust horizontal plot region limits
          order=yi,                ### order by size of yi
          slab=NA, annotate=FALSE, ### remove study labels and annotations
          efac=0,                  ### remove vertical bars at end of CIs
          pch=19,                  ### changing point symbol to filled circle
          col="gray40",            ### change color of points/CIs
          psize=2,                 ### increase point size
          cex.lab=1, cex.axis=1,   ### increase size of x-axis title/labels
          lty=c("solid","blank"))  ### remove horizontal line at top of plot

   ### draw points one more time to make them easier to see
   points(sort(yi), k:1, pch=19, cex=0.5)

   ### add summary polygon at bottom and text
   addpoly(res, mlab="", cex=1)
   text(-2, -2, "RE Model", pos=4, offset=0, cex=1)

   par(opar)

})
