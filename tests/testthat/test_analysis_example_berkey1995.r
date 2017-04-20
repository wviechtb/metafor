### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:berkey1995

context("Checking analysis example: berkey1995")

### load BCG dataset
data(dat.bcg, package="metafor")

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
   expect_equivalent(round(coef(res.RE),4), -0.5429)
   expect_equivalent(round(res.RE$se,4), 0.1842)
   expect_equivalent(round(res.RE$tau2,3), 0.268)

})

test_that("results are correct for the mixed-effects meta-regression model.", {

   ### fit random-effects model using empirical Bayes method
   res.RE <- rma(yi, vi, data=dat, method="EB")

   ### fit mixed-effects model with absolute latitude as moderator
   res.ME <- rma(yi, vi, mods=~I(ablat-33.46), data=dat, method="EB")
   out <- capture.output(print(res.ME))

   ### compare with results on page 408
   expect_equivalent(round(coef(res.ME),4), c(-0.6303, -0.0268)) ### -0.6304 in article
   expect_equivalent(round(res.ME$se,4), c(0.1591, 0.0110))
   expect_equivalent(round(res.ME$tau2,3), 0.157)
   expect_equivalent(round(anova(res.RE, res.ME)$R2,0), 41)

   ### predicted average risk ratios
   tmp <- predict(res.ME, newmods=c(33.46,42)-33.46, transf=exp, digits=2)

   ### compare with results on page 408
   expect_equivalent(round(tmp$pred,2), c(0.53, 0.42))

})

test_that("results are correct for the fixed-effects meta-regression model.", {

   ### fit fixed-effects model with absolute latitude as moderator
   res.FE <- rma(yi, vi, mods=~I(ablat-33.46), data=dat, method="FE")

   ### compare with results on page 408
   expect_equivalent(round(coef(res.FE),4), c(-0.5949, -0.0282)) ### -0.5950 in article
   expect_equivalent(round(res.FE$se,4), c(0.0696, 0.0040)) ### 0.0039 in article

   ### predicted risk ratios based on the fixed-effects model
   tmp <- predict(res.FE, newmods=c(33.46,42)-33.46, transf=exp, digits=2)

   ### compare with results on page 408
   expect_equivalent(round(tmp$pred,2), c(0.55, 0.43))

})
