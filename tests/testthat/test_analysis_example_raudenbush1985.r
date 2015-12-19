### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:raudenbush1985

context("Checking analysis example raudenbush1985")

### load data
dat <- get(data(dat.raudenbush1985, package="metafor"))

test_that("results are correct for the random-effects model.", {

   ### random-effects model
   res <- rma(yi, vi, data=dat, digits=3)

   ### compare with results on pages 83, 85, and 86 (in text)
   expect_equivalent(round(res$tau2,3), 0.019)
   expect_equivalent(round(coef(res),3), 0.084)
   expect_equivalent(round(res$QE,2), 35.83) ### 35.85 in paper
   expect_equivalent(round(res$zval,2), 1.62)

   ### empirical Bayes estimates
   tmp <- blup(res)

   out <- capture.output(print(tmp)) ### so that print.list.rma() is run (at least once)

   ### compare with results in Figure 2
   expect_equivalent(round(tmp$pred,3), c(0.054, 0.101, -0.006, 0.214, 0.105, -0.008, 0.017, -0.029, 0.16, 0.249, 0.162, 0.11, 0.065, 0.11, -0.029, 0.026, 0.191, 0.074, 0.025))
   expect_equivalent(round(tmp$pi.lb,3), c(-0.132, -0.103, -0.223, -0.053, -0.162, -0.174, -0.148, -0.269, -0.054, 0, -0.097, -0.13, -0.192, -0.146, -0.241, -0.191, -0.008, -0.081, -0.195))
   expect_equivalent(round(tmp$pi.ub,3), c(0.241, 0.304, 0.21, 0.482, 0.372, 0.157, 0.183, 0.21, 0.375, 0.497, 0.421, 0.351, 0.321, 0.367, 0.183, 0.242, 0.389, 0.23, 0.245))

   skip_on_cran()

   ### profile tau^2
   opar <- par(no.readonly=TRUE)
   profile(res, xlim=c(0,.20), progbar=FALSE)
   par(opar)

})

test_that("results are correct for the mixed-effects model.", {

   ### recode weeks variable
   dat$weeks.c <- ifelse(dat$weeks > 3, 3, dat$weeks)

   ### mixed-effects model
   res <- rma(yi, vi, mods = ~ weeks.c, data=dat, digits=3)

   ### compare with results on pages 90 and 92 (in text)
   expect_equivalent(round(res$tau2,3), 0)
   expect_equivalent(round(coef(res),3), c(0.407, -0.157))
   expect_equivalent(round(res$QE,2), 16.57) ### 16.58 in paper
   expect_equivalent(round(res$zval,2), c(4.68, -4.39))

   ### empirical Bayes estimates
   tmp <- blup(res)

   ### (results for this not given in chapter)
   expect_equivalent(round(tmp$pred,3), c(0.093, -0.065, -0.065, 0.407, 0.407, -0.065, -0.065, -0.065, 0.407, 0.25, 0.407, 0.407, 0.25, 0.093, -0.065, -0.065, 0.25, 0.093, -0.065))
   expect_equivalent(round(tmp$pi.lb,3), c(0.02, -0.155, -0.155, 0.237, 0.237, -0.155, -0.155, -0.155, 0.237, 0.139, 0.237, 0.237, 0.139, 0.02, -0.155, -0.155, 0.139, 0.02, -0.155))
   expect_equivalent(round(tmp$pi.ub,3), c(0.166, 0.026, 0.026, 0.578, 0.578, 0.026, 0.026, 0.026, 0.578, 0.361, 0.578, 0.578, 0.361, 0.166, 0.026, 0.026, 0.361, 0.166, 0.026))

   ### predicted/fitted values
   tmp <- predict(res)

   ### (results for this not given in chapter)
   expect_equivalent(round(tmp$pred,3), c(0.093, -0.065, -0.065, 0.407, 0.407, -0.065, -0.065, -0.065, 0.407, 0.25, 0.407, 0.407, 0.25, 0.093, -0.065, -0.065, 0.25, 0.093, -0.065))
   expect_equivalent(round(tmp$ci.lb,3), c(0.02, -0.155, -0.155, 0.237, 0.237, -0.155, -0.155, -0.155, 0.237, 0.139, 0.237, 0.237, 0.139, 0.02, -0.155, -0.155, 0.139, 0.02, -0.155))
   expect_equivalent(round(tmp$ci.ub,3), c(0.166, 0.026, 0.026, 0.578, 0.578, 0.026, 0.026, 0.026, 0.578, 0.361, 0.578, 0.578, 0.361, 0.166, 0.026, 0.026, 0.361, 0.166, 0.026))

   skip_on_cran()

   ### profile tau^2
   opar <- par(no.readonly=TRUE)
   profile(res, xlim=c(0,.06), progbar=FALSE)
   par(opar)

})
