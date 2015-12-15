### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:plot_of_influence_diagnostics

context("Checking plots example: plot of influence diagnostics")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### load validity of employment interviews data
   data(dat.mcdaniel1994)

   ### fit random-effects model with r-to-z transformed correlations
   res <- rma(ri=ri, ni=ni, measure="ZCOR", data=dat.mcdaniel1994)

   ### calculate influence diagnostics
   inf <- influence(res)

   print(inf) ### so that print.infl.rma.uni() is run (at least once)

   ### plot the influence diagnostics
   plot(inf, layout=c(8,1))

   par(opar)

})
