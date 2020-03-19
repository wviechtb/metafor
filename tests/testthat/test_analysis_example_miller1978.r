### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:miller1978

context("Checking analysis example: miller1978")

source("tolerances.r") # read in tolerances

### create dataset
dat <- data.frame(xi=c(3, 6, 10, 1), ni=c(11, 17, 21, 6))
dat$pi <- with(dat, xi/ni)
dat <- escalc(measure="PFT", xi=xi, ni=ni, data=dat)

test_that("calculations of escalc() for measure='PFT' are correct.", {

   ### compare with results on page 138
   expect_equivalent(dat$yi*2, c(1.1391, 1.2888, 1.5253, 0.9515), tolerance=.tol[["est"]]) ### need *2 factor due to difference in definition of measure
   expect_equivalent(dat$vi*4, c(0.0870, 0.0571, 0.0465, 0.1538), tolerance=.tol[["var"]])

})

test_that("results are correct for the fixed-effects model using unweighted estimation.", {

   res <- rma(yi, vi, method="FE", data=dat, weighted=FALSE)

   pred <- predict(res, transf=function(x) x*2)
   expect_equivalent(pred$pred, 1.2262, tolerance=.tol[["pred"]])

   pred <- predict(res, transf=transf.ipft.hm, targs=list(ni=dat$ni))
   expect_equivalent(pred$pred, 0.3164, tolerance=.tol[["pred"]])

})

test_that("results are correct for the fixed-effects model using weighted estimation.", {

   res <- rma(yi, vi, method="FE", data=dat)

   pred <- predict(res, transf=function(x) x*2)
   expect_equivalent(pred$pred, 1.3093, tolerance=.tol[["pred"]])

   pred <- predict(res, transf=transf.ipft.hm, targs=list(ni=dat$ni))
   expect_equivalent(pred$pred, 0.3595, tolerance=.tol[["pred"]])

})

test_that("results are correct when there are proportions of 0 and 1.", {

   ### create dataset
   dat <- data.frame(xi=c(0,10), ni=c(10,10))
   dat$pi <- with(dat, xi/ni)
   dat <- escalc(measure="PFT", xi=xi, ni=ni, data=dat, add=0)

   ### back-transformation of the individual outcomes
   expect_equivalent(transf.ipft(dat$yi, dat$ni), c(0,1))

})

test_that("back-transformations work as intended for individual studies and the model estimate.", {

   ### create dataset
   dat <- data.frame(xi = c( 0,  4,  9, 16, 20),
                     ni = c(10, 10, 15, 20, 20))
   dat$pi <- with(dat, xi/ni)
   dat <- escalc(measure="PFT", xi=xi, ni=ni, data=dat, add=0)

   ### back-transformation of the individual outcomes
   expect_equivalent(transf.ipft(dat$yi, dat$ni), c(0.0, 0.4, 0.6, 0.8, 1.0))

   ### back-transformation of the estimated average
   res <- rma(yi, vi, method="FE", data=dat)
   pred <- predict(res, transf=transf.ipft.hm, targs=list(ni=dat$ni))

   expect_equivalent(pred$pred,  0.6886, tolerance=.tol[["pred"]])
   expect_equivalent(pred$ci.lb, 0.5734, tolerance=.tol[["ci"]])
   expect_equivalent(pred$ci.ub, 0.7943, tolerance=.tol[["ci"]])

   ### calculate back-transformed CI bounds
   dat.back <- summary(dat, transf=transf.ipft, ni=dat$ni)

   skip_on_cran()

   ### create forest plot with CI bounds supplied and then add model estimate
   opar <- par(no.readonly=TRUE)
   forest(dat.back$yi, ci.lb=dat.back$ci.lb, ci.ub=dat.back$ci.ub, psize=1,
          xlim=c(-.5,1.8), alim=c(0,1), ylim=c(-1,8), refline=NA, digits=3,
          xlab="Proportion", header=c("Study", "Proportion [95% CI]"))
   addpoly(pred$pred, ci.lb=pred$ci.lb, ci.ub=pred$ci.ub, rows=-0.5, digits=3, mlab="FE Model", efac=1.3)
   abline(h=0.5)
   par(opar)

})
