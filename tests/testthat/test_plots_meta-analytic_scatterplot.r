### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/plots:meta_analytic_scatterplot

context("Checking plots example: meta-analytic scatterplot")

test_that("plot can be drawn.", {

   skip_on_cran()

   opar <- par()

   ### adjust margins so the space is better used
   par(mar=c(5,5,1,2))

   ### load BCG vaccine data
   data(dat.bcg)

   ### calculate (log) relative risks and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit mixed-effects model with absolute latitude as predictor
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   ### calculate predicted relative risks for 0 to 60 degrees absolute latitude
   preds <- predict(res, newmods=c(0:60), transf=exp)

   ### calculate point sizes by rescaling the standard errors
   wi    <- 1/sqrt(dat$vi)
   size  <- 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))

   ### plot the relative risks against absolute latitude
   plot(dat$ablat, exp(dat$yi), pch=19, cex=size,
        xlab="Absolute Latitude", ylab="Relative Risk",
        las=1, bty="l", log="y")

   ### add predicted values (and corresponding CI bounds)
   lines(0:60, preds$pred)
   lines(0:60, preds$ci.lb, lty="dashed")
   lines(0:60, preds$ci.ub, lty="dashed")

   ### dotted line at RR=1 (no difference between groups)
   abline(h=1, lty="dotted")

   ### labels some points in the plot
   ids <- c(4,7,12,13)
   pos <- c(3,3,1,1)
   text(dat$ablat[ids], exp(dat$yi)[ids], ids, cex=.9, pos=pos)

   par(opar)

})
