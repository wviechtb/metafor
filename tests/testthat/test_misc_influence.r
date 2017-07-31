### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: influence() and related functions")

test_that("influence() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   inf <- influence(res)
   inf$inf <- round(inf$inf[1,,drop=FALSE], 4)
   inf$dfbs <- round(inf$dfbs[1,,drop=FALSE], 4)
   inf$tau2 <- round(inf$tau2, 4)
   inf$QE <- round(inf$QE, 4)
   inf$is.infl <- inf$is.infl[1]
   inf$not.na <- inf$not.na[1]

   sav <- structure(list(inf = structure(list(rstudent = -0.2181, dffits = -0.0407, cook.d = 0.0017, cov.r = 1.1164, tau2.del = 0.3362, QE.del = 151.5826, hat = 0.0506, weight = 5.0595), .Names = c("rstudent", "dffits", "cook.d", "cov.r", "tau2.del", "QE.del", "hat", "weight"), row.names = 1L, class = "data.frame"), dfbs = structure(list(intrcpt = -0.0403), .Names = "intrcpt", row.names = 1L, class = "data.frame"), tau2 = 0.3132, QE = 152.233, ids = 1:13, not.na = TRUE, is.infl = FALSE, k = 13L, p = 1L, digits = 4), .Names = c("inf", "dfbs", "tau2", "QE", "ids", "not.na", "is.infl", "k", "p", "digits"), class = "infl.rma.uni")

   expect_equivalent(sav, inf)

})

test_that("leave1out() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   inf <- leave1out(res)
   inf <- lapply(inf[1:12], function(x) round(x[1], 4))
   class(inf) <- "list.rma"

   sav <- structure(list(estimate = -0.7071, se = 0.19, zval = -3.7223, pval = 2e-04, ci.lb = -1.0794, ci.ub = -0.3348, Q = 151.5826, Qp = 0, tau2 = 0.3362, I2 = 93.2259, H2 = 14.7622, slab = 1, digits = 4, transf = 0), .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub", "Q", "Qp", "tau2", "I2", "H2", "slab", "digits", "transf"), class = "list.rma")

   expect_equivalent(sav, inf)

})

test_that("leave1out() works for rma.mh().", {

   data(dat.bcg, package="metafor")
   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   inf <- leave1out(res)
   inf <- lapply(inf[1:12], function(x) round(x[1], 4))
   class(inf) <- "list.rma"

   sav <- structure(list(estimate = -0.4514, se = 0.0394, zval = -11.4462, pval = 0, ci.lb = -0.5287, ci.ub = -0.3741, Q = 151.9153, Qp = 0, slab = 1, digits = 4, transf = 0), .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub", "Q", "Qp", "slab", "digits", "transf"), class = "list.rma")

   expect_equivalent(sav, inf)

})

test_that("leave1out() works for rma.peto().", {

   data(dat.bcg, package="metafor")
   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   inf <- leave1out(res)
   inf <- lapply(inf[1:12], function(x) round(x[1], 4))
   class(inf) <- "list.rma"

   sav <- structure(list(estimate = -0.4722, se = 0.0408, zval = -11.5791, pval = 0, ci.lb = -0.5521, ci.ub = -0.3923, Q = 167.2005, Qp = 0, slab = 1, digits = 4, transf = 0), .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub", "Q", "Qp", "slab", "digits", "transf"), class = "list.rma")

   expect_equivalent(sav, inf)

})

test_that("model.matrix() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   sav <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 44, 55, 42, 52, 13, 44, 19, 13, 27, 42, 18, 33, 33), .Dim = c(13L, 2L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"), c("intrcpt", "ablat")))

   expect_equivalent(sav, model.matrix(res))

})

test_that("hatvalues() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   expect_equivalent(round(hatvalues(res), 4), c(0.049, 0.1493, 0.0351, 0.3481, 0.2248, 0.2367, 0.064, 0.357, 0.0926, 0.1157, 0.2309, 0.0189, 0.0778))

   sav <- structure(c(0.049, 0.067, 0.0458, 0.0994, 0.1493, 0.0904, 0.0374, 0.0498, 0.0351), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"), c("1", "2", "3")))
   expect_equivalent(round(hatvalues(res, type="matrix")[1:3,1:3], 4), sav)

})

test_that("hatvalues() works for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, mods = ~ ablat, random = ~ 1 | trial, data=dat)

   expect_equivalent(round(hatvalues(res), 4), c(0.049, 0.1493, 0.0351, 0.3481, 0.2248, 0.2367, 0.064, 0.357, 0.0926, 0.1157, 0.2309, 0.0189, 0.0778))

   sav <- structure(c(0.049, 0.067, 0.0458, 0.0994, 0.1493, 0.0904, 0.0374, 0.0498, 0.0351), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"), c("1", "2", "3")))
   expect_equivalent(round(hatvalues(res, type="matrix")[1:3,1:3], 4), sav)

})

test_that("cooks.distance() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   expect_equivalent(round(cooks.distance(res), 4), c(0.0048, 0.0489, 0.0104, 0.2495, 0.0072, 0.2883, 0.3643, 0.2719, 0.02, 0.1645, 9e-04, 0.0403, 0.1433))

})

test_that("cooks.distance() works for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, mods = ~ ablat, random = ~ 1 | trial, data=dat)

   expect_equivalent(round(cooks.distance(res), 4), c(0.0048, 0.0489, 0.0104, 0.2495, 0.0072, 0.2883, 0.3643, 0.2719, 0.02, 0.1645, 9e-04, 0.0404, 0.1434))
   expect_equivalent(round(cooks.distance(res, cluster=dat$alloc), 4), c(2.4372, 0.2591, 0.1533))
   expect_equivalent(round(cooks.distance(res, cluster=dat$alloc, reestimate=FALSE), 4), c(2.2194, 0.3199, 0.2421))

})
