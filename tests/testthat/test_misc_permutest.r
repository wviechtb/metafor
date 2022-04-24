### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: permutest() function")

source("settings.r")

### load data
dat <- dat.hine1989

### calculate risk differences and corresponding sampling variances
dat <- escalc(measure="RD", n1i=n1i, n2i=n2i, ai=ai, ci=ci, data=dat)

test_that("permutest() gives correct results for a random-effects model.", {

   skip_on_cran()

   ### fit random-effects model
   res <- rma(yi, vi, data=dat)

   ### exact permutation test
   sav <- permutest(res, progbar=FALSE)

   expect_equivalent(sav$pval, 0.0625)

   out <- capture.output(print(sav)) ### so that print.permutest.rma.uni() is run (at least once)

   tmp <- coef(sav)
   expected <- structure(list(estimate = 0.029444, se = 0.013068, zval = 2.253107, pval = 0.0625, ci.lb = 0.003831, ci.ub = 0.055058),
                         .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub"), row.names = "intrcpt", class = "data.frame")

   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

   ### approximate permutation test
   set.seed(1234)
   sav <- permutest(res, iter=50, progbar=FALSE, control=list(p2defn="px2"))
   expect_equivalent(sav$pval, 0.08)
   set.seed(1234)
   sav <- permutest(res, iter=50, progbar=FALSE, control=list(p2defn="px2", stat="coef"))
   expect_equivalent(sav$pval, 0.08)

})

test_that("permutest() gives correct results for a mixed-effects model.", {

   skip_on_cran()

   ### add a fake moderator
   dat$mod <- c(3,1,2,2,4,5)

   ### fit mixed-effects model
   res <- rma(yi, vi, mods = ~ mod, data=dat)

   ### exact permutation test
   sav <- permutest(res, progbar=FALSE)

   expect_equivalent(sav$pval, c(1, 0.0028), tolerance=.tol[["pval"]])

   ### approximate permutation test
   set.seed(1234)
   sav <- permutest(res, iter=50, progbar=FALSE, control=list(p2defn="px2"))
   expect_equivalent(sav$pval, c(.04, .04))
   sav <- permutest(res, iter=50, progbar=FALSE, control=list(p2defn="px2", stat="coef"))
   expect_equivalent(sav$pval, c(.04, .04))

})

test_that("permutest() gives correct results for example in Follmann & Proschan (1999).", {

   skip_on_cran()

   ### data in Table 1
   dat <- read.table(header=TRUE, text = "
   ai n1i ci n2i
   173 5331 210 5296
   157 1906 193 1900
   131 4541 121 4516
   56 2051 84 2030
   52 424 65 422
   36 1149 42 1129
   62 6582 20 1663
   2 88 2 30")

   dat <- escalc(measure="PETO", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat)

   res <- rma(yi, vi, data=dat, method="DL")
   sav <- permutest(res, permci=TRUE, progbar=FALSE, control=list(stat="coef"))

   expect_equivalent(sav$pval, 10/256)
   expect_equivalent(sav$ci.lb, -0.3677, tolerance=.tol[["ci"]])
   expect_equivalent(sav$ci.ub, -0.0020, tolerance=.tol[["ci"]])

})

rm(list=ls())
