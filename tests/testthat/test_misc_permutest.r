### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: permutest() function")

### load data
dat <- get(data(dat.hine1989))

### calculate risk differences and corresponding sampling variances
dat <- escalc(measure="RD", n1i=n1i, n2i=n2i, ai=ai, ci=ci, data=dat)

test_that("permutest() gives correct results for a random-effects model", {

   ### fit random-effects model
   res <- rma(yi, vi, data=dat)

   ### exact permutation test
   sav <- permutest(res, progbar=FALSE)

   expect_equivalent(sav$pval, 0.0625)

   out <- capture.output(print(sav)) ### so that print.permutest.rma.uni() is run (at least once)

   tmp <- round(coef(sav), 4)
   expected <- structure(list(estimate = 0.0294, se = 0.0131, zval = 2.2531, pval = 0.0625, ci.lb = 0.0038, ci.ub = 0.0551),
                         .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub"), row.names = "intrcpt", class = "data.frame")

   expect_equivalent(tmp, expected)

})

test_that("permutest() gives correct results for a mixed-effects model", {

   skip_on_cran()

   ### add a fake moderator
   dat$mod <- c(3,1,2,2,4,5)

   ### fit mixed-effects model
   res <- rma(yi, vi, mods = ~ mod, data=dat)

   ### exact permutation test
   sav <- permutest(res, progbar=FALSE)

   expect_equivalent(round(sav$pval, 4), c(1, 0.0028))

})
