### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:yusuf1985

context("Checking analysis example: yusuf1985")

source("tolerances.r") # read in tolerances

### create dataset for example
dat <- dat.yusuf1985
dat$grp_ratios <- round(dat$n1i / dat$n2i, 2)

test_that("log likelihood plot can be drawn.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)
   par(mfrow=c(1,2))
   expect_warning(llplot(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
                  subset=(table=="6"), drop00=FALSE, lwd=1, xlim=c(-5,5)))
   expect_warning(llplot(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
                  subset=(table=="6"), drop00=FALSE, lwd=1, xlim=c(-5,5), scale=FALSE))
   par(opar)

})

test_that("results are correct for the analysis using Peto's method.", {

   expect_warning(res <- rma.peto(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, subset=(table=="6")))
   out <- capture.output(print(res)) ### so that print.rma.peto() is run (at least once)
   out <- capture.output(print(summary(res))) ### so that print.rma.peto() is run (at least once) with showfit=TRUE

   sav <- predict(res, transf=exp)
   tmp <- c(sav$pred, sav$ci.lb, sav$ci.ub)

   ### compare with results on page 107
   expect_equivalent(tmp, c(.9332, .7385, 1.1792), tolerance=.tol[["pred"]])

})

test_that("results are correct for the analysis using the inverse-variance method.", {

   expect_warning(dat <- escalc(measure="PETO", ai=ai, n1i=n1i, ci=ci, n2i=n2i,
                         data=dat, subset=(table=="6"), add=0))

   expect_warning(res <- rma(yi, vi, data=dat, method="FE"))
   sav <- predict(res, transf=exp)
   tmp <- c(sav$pred, sav$ci.lb, sav$ci.ub)

   ### compare with results on page 107
   expect_equivalent(tmp, c(.9332, .7385, 1.1792), tolerance=.tol[["pred"]])

})
