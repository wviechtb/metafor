### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:stijnen2010

context("Checking analysis example: stijnen2010")

source("settings.r")

### load data
dat <- dat.nielweise2007

test_that("results for the normal-normal model are correct (measure=='PLO')", {

   res <- rma(measure="PLO", xi=ci, ni=n2i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(coef(res), -3.3018, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2378, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.6629, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(tmp$pred, 0.0355, tolerance=.tol[["pred"]]) ### 0.035 in paper
   expect_equivalent(tmp$ci.lb, 0.0226, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.0554, tolerance=.tol[["ci"]]) ### 0.056 in paper

   res <- rma(measure="PLO", xi=ai, ni=n1i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(coef(res), -4.2604, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2589, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.3928, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(tmp$pred, 0.0139, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.0084, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.0229, tolerance=.tol[["ci"]])

})

test_that("results for the binomial-normal normal are correct (measure=='PLO')", {

   skip_on_cran()

   res <- rma.glmm(measure="PLO", xi=ci, ni=n2i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(coef(res), -3.4964, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2570, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.8124, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(tmp$pred, 0.0294, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.0180, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.0478, tolerance=.tol[["ci"]])

   res <- rma.glmm(measure="PLO", xi=ai, ni=n1i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(coef(res), -4.8121, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.3555, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.8265, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(tmp$pred, 0.0081, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.0040, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.0161, tolerance=.tol[["ci"]])

})

test_that("results for the normal-normal model are correct (measure=='OR')", {

   expect_warning(res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, drop00=TRUE))

   ### compare with results on page 3052 (Table III)
   expect_equivalent(coef(res), -0.9804, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2435, tolerance=.tol[["se"]]) ### 0.244 in paper
   expect_equivalent(sqrt(res$tau2), 0.1886, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 0.3752, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.2328, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.6046, tolerance=.tol[["ci"]]) ### 0.62 in paper

})

test_that("results for the conditional logistic model with exact likelihood are correct (measure=='OR')", {

   skip_on_cran()

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL"))
   out <- capture.output(print(res)) ### so that print.rma.glmm() is run (at least once)

   ### compare with results on page 3052 (Table III)
   expect_equivalent(coef(res), -1.3532, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.3511, tolerance=.tol[["se"]])
   expect_equivalent(sqrt(res$tau2), 0.8327, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 0.2584, tolerance=.tol[["pred"]]) ### 0.25 in paper
   expect_equivalent(tmp$ci.lb, 0.1299, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.5142, tolerance=.tol[["ci"]])

})

test_that("results for the conditional logistic model with approximate likelihood are correct (measure=='OR')", {

   skip_on_cran()

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.AL"))

   ### compare with results on page 3052 (Table III)
   expect_equivalent(coef(res), -1.3027, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.3386, tolerance=.tol[["se"]])
   expect_equivalent(sqrt(res$tau2), 0.7750, tolerance=.tol[["var"]]) ### 0.77 in paper
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 0.2718, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.1400, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.5279, tolerance=.tol[["ci"]])

})

############################################################################

### load data
dat <- dat.nielweise2008

### incidence rates reflect the expected number of events per 1000 days
dat$t1i <- dat$t1i/1000
dat$t2i <- dat$t2i/1000

test_that("results for the normal-normal model are correct (measure=='IRLN')", {

   res <- rma(measure="IRLN", xi=x2i, ti=t2i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(coef(res), 1.4676, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2425, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.3699, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 4.3389, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 2.6973, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 6.9795, tolerance=.tol[["ci"]]) ### 6.99 in paper

   res <- rma(measure="IRLN", xi=x1i, ti=t1i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(coef(res), 0.9808, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.3259, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.6393, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 2.6667, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 1.4078, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 5.0513, tolerance=.tol[["ci"]])

})

test_that("results for the Poisson-normal model are correct (measure=='IRLN')", {

   skip_on_cran()

   res <- rma.glmm(measure="IRLN", xi=x2i, ti=t2i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(coef(res), 1.4007, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2310, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.3165, tolerance=.tol[["var"]]) ### 0.316 in paper
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 4.0580, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 2.5803, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 6.3819, tolerance=.tol[["ci"]])

   res <- rma.glmm(measure="IRLN", xi=x1i, ti=t1i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(coef(res), 0.8494, tolerance=.tol[["coef"]]) ### 0.850 in paper
   expect_equivalent(res$se, 0.3303, tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.6543, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 2.3383, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 1.2240, tolerance=.tol[["ci"]]) ### 1.23 in paper
   expect_equivalent(tmp$ci.ub, 4.4670, tolerance=.tol[["ci"]])

})

test_that("results for the normal-normal model are correct (measure=='IRR')", {

   res <- rma(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat)

   ### compare with results on page 3055 (Table VIII)
   expect_equivalent(coef(res), -0.3963, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2268, tolerance=.tol[["se"]]) ### 0.223 in paper
   expect_equivalent(sqrt(res$tau2), 0.3060, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 0.6728, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.4314, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 1.0494, tolerance=.tol[["ci"]]) ### 1.04 in paper

})

test_that("results for the Poisson-normal model are correct (measure=='IRR')", {

   skip_on_cran()

   res <- rma.glmm(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat, model="CM.EL")

   ### compare with results on page 3055 (Table VIII)
   expect_equivalent(coef(res), -0.4762, tolerance=.tol[["coef"]])
   expect_equivalent(res$se, 0.2377, tolerance=.tol[["se"]])
   expect_equivalent(sqrt(res$tau2), 0.3501, tolerance=.tol[["var"]])
   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred, 0.6211, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.3898, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.9897, tolerance=.tol[["ci"]])

})

rm(list=ls())
