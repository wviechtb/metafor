### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:normal_qq_plots

context("Checking plots example: normal QQ plots")

test_that("plot can be drawn for rma().", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### set up 2x2 array for plotting
   par(mfrow=c(2,2))

   ### load BCG vaccine data
   data(dat.bcg, package="metafor")

   ### calculate log risk ratios and corresponding sampling variances
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

   ### draw plot with studentized residuals and labels
   qqnorm(res2, type="rstudent", label=TRUE)

})

test_that("plot can be drawn for rma.mh().", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)
   data(dat.bcg, package="metafor")
   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   qqnorm(res)
   qqnorm(res, type="rstudent", label=TRUE)
   par(opar)

})

test_that("plot can be drawn for rma.peto().", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)
   data(dat.bcg, package="metafor")
   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   qqnorm(res)
   qqnorm(res, type="rstudent", label=TRUE)
   par(opar)

})

test_that("plot cannot be drawn for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat)
   expect_error(qqnorm(res))

})
