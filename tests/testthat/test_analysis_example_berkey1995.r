### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:berkey1995

source("tolerances.r") # read in tolerances

context("Checking analysis example: berkey1995")

### calculate log ratio ratios and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### calculate "smoothed" sampling variances
dat$vi <- with(dat, sum(tneg/tpos)/(13*(tneg+tpos)) +
                    sum(cneg/cpos)/(13*(cneg+cpos)))

test_that("results are correct for the random-effects model.", {

   ### fit random-effects model using empirical Bayes method
   res.RE <- rma(yi, vi, data=dat, method="EB")
   out <- capture.output(print(res.RE)) ### so that print.rma.uni() is run (at least once)
   out <- capture.output(print(summary(res.RE))) ### so that print.summary.rma() is run (at least once)

   ### compare with results on page 408
   expect_equivalent(coef(res.RE), -0.5429, tolerance=.tol[["coef"]])
   expect_equivalent(res.RE$se, 0.1842, tolerance=.tol[["se"]])
   expect_equivalent(res.RE$tau2, 0.2682, tolerance=.tol[["var"]])

})

test_that("results are correct for the mixed-effects meta-regression model.", {

   ### fit random-effects model using empirical Bayes method
   res.RE <- rma(yi, vi, data=dat, method="EB")

   ### fit mixed-effects model with absolute latitude as moderator
   res.ME <- rma(yi, vi, mods=~I(ablat-33.46), data=dat, method="EB")
   out <- capture.output(print(res.ME))

   ### compare with results on page 408
   expect_equivalent(coef(res.ME), c(-0.6303, -0.0268), tolerance=.tol[["coef"]]) ### -0.6304 in article
   expect_equivalent(res.ME$se, c(0.1591, 0.0110), tolerance=.tol[["se"]])
   expect_equivalent(res.ME$tau2, 0.1572, tolerance=.tol[["var"]])
   expect_equivalent(anova(res.RE, res.ME)$R2, 41.3844, tolerance=.tol[["r2"]])

   ### predicted average risk ratios
   tmp <- predict(res.ME, newmods=c(33.46,42)-33.46, transf=exp, digits=2)

   ### compare with results on page 408
   expect_equivalent(tmp$pred, c(0.5324, 0.4236), tolerance=.tol[["pred"]])

})

test_that("results are correct for the fixed-effects meta-regression model.", {

   ### fit fixed-effects model with absolute latitude as moderator
   res.FE <- rma(yi, vi, mods=~I(ablat-33.46), data=dat, method="FE")

   ### compare with results on page 408
   expect_equivalent(coef(res.FE), c(-0.5949, -0.0282), tolerance=.tol[["coef"]]) ### -0.5950 in article
   expect_equivalent(res.FE$se, c(0.0696, 0.0040), tolerance=.tol[["se"]]) ### 0.0039 in article

   ### predicted risk ratios based on the fixed-effects model
   tmp <- predict(res.FE, newmods=c(33.46,42)-33.46, transf=exp, digits=2)

   ### compare with results on page 408
   expect_equivalent(tmp$pred, c(0.5516, 0.4336), tolerance=.tol[["pred"]])

})
