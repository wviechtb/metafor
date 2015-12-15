### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:stijnen2010

context("Checking analysis example: stijnen2010")

### load data
dat <- get(data(dat.nielweise2007, package="metafor"))

test_that("results for the normal-normal model are correct (measure=='PLO')", {

   res <- rma(measure="PLO", xi=ci, ni=n2i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(round(coef(res), digits=2), -3.30)
   expect_equivalent(round(res$se, digits=2), 0.24)
   expect_equivalent(round(res$tau2, digits=3), 0.663)
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(round(tmp$pred, digits=3), 0.036) ### 0.035 in paper
   expect_equivalent(round(tmp$ci.lb, digits=3), 0.023)
   expect_equivalent(round(tmp$ci.ub, digits=3), 0.055) ### 0.056 in paper

   res <- rma(measure="PLO", xi=ai, ni=n1i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(round(coef(res), digits=2), -4.26)
   expect_equivalent(round(res$se, digits=2), 0.26)
   expect_equivalent(round(res$tau2, digits=3), 0.393)
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(round(tmp$pred, digits=3), 0.014)
   expect_equivalent(round(tmp$ci.lb, digits=3), 0.008)
   expect_equivalent(round(tmp$ci.ub, digits=3), 0.023)

})

test_that("results for the binomial-normal normal are correct (measure=='PLO')", {

   res <- rma.glmm(measure="PLO", xi=ci, ni=n2i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(round(coef(res), digits=2), -3.50)
   expect_equivalent(round(res$se, digits=2), 0.26)
   expect_equivalent(round(res$tau2, digits=2), 0.81)
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(round(tmp$pred, digits=3), 0.029)
   expect_equivalent(round(tmp$ci.lb, digits=3), 0.018)
   expect_equivalent(round(tmp$ci.ub, digits=3), 0.048)

   res <- rma.glmm(measure="PLO", xi=ai, ni=n1i, data=dat)

   ### compare with results on page 3050 (Table II)
   expect_equivalent(round(coef(res), digits=2), -4.81)
   expect_equivalent(round(res$se, digits=2), 0.36)
   expect_equivalent(round(res$tau2, digits=2), 0.83)
   tmp <- predict(res, transf=transf.ilogit)
   expect_equivalent(round(tmp$pred, digits=3), 0.008)
   expect_equivalent(round(tmp$ci.lb, digits=3), 0.004)
   expect_equivalent(round(tmp$ci.ub, digits=3), 0.016)

})

test_that("results for the normal-normal model are correct (measure=='OR')", {

   res <- rma(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, drop00=TRUE)

   ### compare with results on page 3052 (Table III)
   expect_equivalent(round(coef(res), digits=3), -0.980)
   expect_equivalent(round(res$se, digits=3), 0.243) ### 0.244 in paper
   expect_equivalent(round(sqrt(res$tau2), digits=2), 0.19)
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 0.38)
   expect_equivalent(round(tmp$ci.lb, digits=2), 0.23)
   expect_equivalent(round(tmp$ci.ub, digits=2), 0.60) ### 0.62 in paper

})

test_that("results for the conditional logistic model with exact likelihood are correct (measure=='OR')", {

   skip_on_cran()

   res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.EL")
   res ### so that print.rma.glmm() is run (at least once)

   ### compare with results on page 3052 (Table III)
   expect_equivalent(round(coef(res), digits=3), -1.353)
   expect_equivalent(round(res$se, digits=3), 0.351)
   expect_equivalent(round(sqrt(res$tau2), digits=2), 0.83)
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 0.26) ### 0.25 in paper
   expect_equivalent(round(tmp$ci.lb, digits=2), 0.13)
   expect_equivalent(round(tmp$ci.ub, digits=2), 0.51)

})

test_that("results for the conditional logistic model with approximate likelihood are correct (measure=='OR')", {

   res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="CM.AL")

   ### compare with results on page 3052 (Table III)
   expect_equivalent(round(coef(res), digits=3), -1.303)
   expect_equivalent(round(res$se, digits=3), 0.339)
   expect_equivalent(round(sqrt(res$tau2), digits=2), 0.78) ### 0.77 in paper
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 0.27)
   expect_equivalent(round(tmp$ci.lb, digits=2), 0.14)
   expect_equivalent(round(tmp$ci.ub, digits=2), 0.53)

})

############################################################################

### load data
dat <- get(data(dat.nielweise2008, package="metafor"))

### incidence rates reflect the expected number of events per 1000 days
dat$t1i <- dat$t1i/1000
dat$t2i <- dat$t2i/1000

test_that("results for the normal-normal model are correct (measure=='IRLN')", {

   res <- rma(measure="IRLN", xi=x2i, ti=t2i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(round(coef(res), digits=3), 1.468)
   expect_equivalent(round(res$se, digits=3), 0.243)
   expect_equivalent(round(res$tau2, digits=3), 0.370)
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 4.34)
   expect_equivalent(round(tmp$ci.lb, digits=2), 2.70)
   expect_equivalent(round(tmp$ci.ub, digits=2), 6.98) ### 6.99 in paper

   res <- rma(measure="IRLN", xi=x1i, ti=t1i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(round(coef(res), digits=3), 0.981)
   expect_equivalent(round(res$se, digits=3), 0.326)
   expect_equivalent(round(res$tau2, digits=3), 0.639)
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 2.67)
   expect_equivalent(round(tmp$ci.lb, digits=2), 1.41)
   expect_equivalent(round(tmp$ci.ub, digits=2), 5.05)

})

test_that("results for the Poisson-normal model are correct (measure=='IRLN')", {

   res <- rma.glmm(measure="IRLN", xi=x2i, ti=t2i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(round(coef(res), digits=3), 1.401)
   expect_equivalent(round(res$se, digits=3), 0.231)
   expect_equivalent(round(res$tau2, digits=3), 0.317) ### 0.316 in paper
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 4.06)
   expect_equivalent(round(tmp$ci.lb, digits=2), 2.58)
   expect_equivalent(round(tmp$ci.ub, digits=2), 6.38)

   res <- rma.glmm(measure="IRLN", xi=x1i, ti=t1i, data=dat)

   ### compare with results on page 3054 (Table VII)
   expect_equivalent(round(coef(res), digits=3), 0.849) ### 0.850 in paper
   expect_equivalent(round(res$se, digits=3), 0.330)
   expect_equivalent(round(res$tau2, digits=3), 0.654)
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 2.34)
   expect_equivalent(round(tmp$ci.lb, digits=2), 1.22) ### 1.23 in paper
   expect_equivalent(round(tmp$ci.ub, digits=2), 4.47)

})

test_that("results for the normal-normal model are correct (measure=='IRR')", {

   res <- rma(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat)

   ### compare with results on page 3055 (Table VIII)
   expect_equivalent(round(coef(res), digits=3), -0.396)
   expect_equivalent(round(res$se, digits=3), 0.227) ### 0.223 in paper
   expect_equivalent(round(sqrt(res$tau2), digits=2), 0.31)
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 0.67)
   expect_equivalent(round(tmp$ci.lb, digits=2), 0.43)
   expect_equivalent(round(tmp$ci.ub, digits=2), 1.05) ### 1.04 in paper

})

test_that("results for the Poisson-normal model are correct (measure=='IRR')", {

   res <- rma.glmm(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat, model="CM.EL")

   ### compare with results on page 3055 (Table VIII)
   expect_equivalent(round(coef(res), digits=3), -0.476)
   expect_equivalent(round(res$se, digits=3), 0.238)
   expect_equivalent(round(sqrt(res$tau2), digits=2), 0.35)
   tmp <- predict(res, transf=exp)
   expect_equivalent(round(tmp$pred, digits=2), 0.62)
   expect_equivalent(round(tmp$ci.lb, digits=2), 0.39)
   expect_equivalent(round(tmp$ci.ub, digits=2), 0.99)

})
