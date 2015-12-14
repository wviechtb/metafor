### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:miller1978

context("Checking analysis example miller1978")

### create dataset
dat <- data.frame(xi=c(3, 6, 10, 1), ni=c(11, 17, 21, 6))
dat$pi <- with(dat, xi/ni)
dat <- escalc(measure="PFT", xi=xi, ni=ni, data=dat)

test_that("calculations of escalc() for measure='PFT' are correct.", {

   ### compare with results on page 138
   expect_equivalent(round(dat$yi*2,4), c(1.1391, 1.2888, 1.5253, 0.9515)) ### need *2 factor due to difference in definition of measure
   expect_equivalent(round(dat$vi*4,4), c(0.0870, 0.0571, 0.0465, 0.1538))

})

test_that("results are correct for the fixed-effects model using unweighted estimation.", {

   res <- rma(yi, vi, method="FE", data=dat, weighted=FALSE)

   pred <- predict(res, transf=function(x) x*2)
   expect_equivalent(round(pred$pred, 4), 1.2262)

   pred <- predict(res, transf=transf.ipft.hm, targs=list(ni=dat$ni))
   expect_equivalent(round(pred$pred, 4), 0.3164)

})

test_that("results are correct for the fixed-effects model using weighted estimation.", {

   res <- rma(yi, vi, method="FE", data=dat)

   pred <- predict(res, transf=function(x) x*2)
   expect_equivalent(round(pred$pred, 4), 1.3093)

   pred <- predict(res, transf=transf.ipft.hm, targs=list(ni=dat$ni))
   expect_equivalent(round(pred$pred, 4), 0.3595)

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

   expect_equivalent(round(pred$pred, 4), 0.6886)
   expect_equivalent(round(pred$ci.lb, 4), 0.5734)
   expect_equivalent(round(pred$ci.ub, 4), 0.7943)

   ### calculate back-transformed CI bounds
   dat.back <- summary(dat, transf=transf.ipft, ni=dat$ni)

   skip_on_cran()

   ### create forest plot with CI bounds supplied and then add model estimate
   opar <- par()
   forest(dat.back$yi, ci.lb=dat.back$ci.lb, ci.ub=dat.back$ci.ub, psize=1,
          xlim=c(-.5,1.8), alim=c(0,1), ylim=c(-1,8), refline=NA, digits=3, xlab="Proportion")
   addpoly(pred$pred, ci.lb=pred$ci.lb, ci.ub=pred$ci.ub, row=-0.5, digits=3, mlab="FE Model", efac=1.3)
   abline(h=0.5)
   text(-0.5, 7, "Study",               pos=4)
   text( 1.8, 7, "Proportion [95% CI]", pos=2)
   par(opar)

})
