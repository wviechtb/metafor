### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:normal_qq_plots

context("Checking plots example: normal QQ plots")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### set up 2x2 array for plotting
   par(mfrow=c(2,2))

   ### load BCG vaccine data
   data(dat.bcg)

   ### calculate (log) relative risks and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit fixed- and random-effects models
   res1 <- rma(yi, vi, data=dat, method="FE")
   res2 <- rma(yi, vi, data=dat)

   ### fit fixed- and random-effects models with absolute latitude moderator
   res3 <- rma(yi, vi, mods=~ablat, data=dat, method="FE")
   res4 <- rma(yi, vi, mods=~ablat, data=dat)

   ### normal QQ plots for the various models
   qqnorm(res1, main="Fixed-Effects Model")
   qqnorm(res2, main="Random-Effects Model")
   qqnorm(res3, main="Fixed-Effects with Moderators Model")
   qqnorm(res4, main="Mixed-Effects Model")

   par(opar)

})
