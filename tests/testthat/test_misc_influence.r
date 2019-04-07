### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: influence() and related functions")

test_that("influence() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   sav <- influence(res)
   sav$inf  <- sav$inf[1]
   sav$dfbs <- sav$dfbs[1]
   sav$inf  <- lapply(sav$inf,  function(x) if (class(x) == "numeric") round(x, 4) else x)
   sav$dfbs <- lapply(sav$dfbs, function(x) if (class(x) == "numeric") round(x, 4) else x)
   sav$tau2 <- round(sav$tau2, 4)
   sav$QE <- round(sav$QE, 4)
   sav$is.infl <- sav$is.infl[1]
   sav$not.na <- sav$not.na[1]

   tmp <- structure(list(inf = list(rstudent = -0.2181, dffits = -0.0407, cook.d = 0.0017, cov.r = 1.1164, tau2.del = 0.3362, QE.del = 151.5826, hat = 0.0506, weight = 5.0595, inf = "", slab = 1L, digits = c( est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), dfbs = list(intrcpt = -0.0403, slab = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), ids = 1:13, not.na = TRUE, is.infl = FALSE, tau2 = 0.3132, QE = 152.233, k = 13L, p = 1L, m = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), class = "infl.rma.uni")

   expect_equivalent(sav, tmp)

})

test_that("leave1out() works for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   inf <- leave1out(res)
   inf <- lapply(unclass(inf)[1:12], function(x) round(x[1], 4))
   class(inf) <- "list.rma"

   sav <- structure(list(estimate = -0.7071, se = 0.19, zval = -3.7223, pval = 2e-04, ci.lb = -1.0794, ci.ub = -0.3348, Q = 151.5826, Qp = 0, tau2 = 0.3362, I2 = 93.2259, H2 = 14.7622, slab = 1), class = "list.rma")

   expect_equivalent(sav, inf)

})

test_that("leave1out() works for rma.mh().", {

   data(dat.bcg, package="metafor")
   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   inf <- leave1out(res)
   inf <- lapply(unclass(inf)[1:9], function(x) round(x[1], 4))
   class(inf) <- "list.rma"

   sav <- structure(list(estimate = -0.4514, se = 0.0394, zval = -11.4462, pval = 0, ci.lb = -0.5287, ci.ub = -0.3741, Q = 151.9153, Qp = 0, slab = 1), class = "list.rma")

   expect_equivalent(sav, inf)

})

test_that("leave1out() works for rma.peto().", {

   data(dat.bcg, package="metafor")
   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   inf <- leave1out(res)
   inf <- lapply(unclass(inf)[1:9], function(x) round(x[1], 4))
   class(inf) <- "list.rma"

   sav <- structure(list(estimate = -0.4722, se = 0.0408, zval = -11.5791, pval = 0, ci.lb = -0.5521, ci.ub = -0.3923, Q = 167.2005, Qp = 0, slab = 1), class = "list.rma")

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

   expect_equivalent(round(cooks.distance(res), 4), c(0.0048, 0.0489, 0.0104, 0.2495, 0.0072, 0.2883, 0.3643, 0.2719, 0.02, 0.1645, 0.0009, 0.0403, 0.1433))

})

test_that("cooks.distance() works for rma.mv().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, mods = ~ ablat, random = ~ 1 | trial, data=dat)

   expect_equivalent(round(cooks.distance(res), 4), c(0.0048, 0.0489, 0.0104, 0.2495, 0.0072, 0.2883, 0.3643, 0.2719, 0.02, 0.1645, 0.0009, 0.0404, 0.1434))
   expect_equivalent(round(cooks.distance(res, cluster=dat$alloc), 4), c(0.2591, 2.4372, 0.1533))
   expect_equivalent(round(cooks.distance(res, cluster=dat$alloc, reestimate=FALSE), 4), c(0.3199, 2.2194, 0.2421))

})

test_that("influence() correctly works with 'na.omit' and 'na.pass'.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, slab=paste0("Trial ", dat.bcg$trial))

   dat$yi[2] <- NA
   dat$vi[3] <- NA
   dat$ablat[5] <- NA
   dat$trial12 <- ifelse(dat$trial == 12, 1, 0)

   options(na.action="na.omit")

   expect_warning(res <- rma(yi, vi, mods = ~ ablat + trial12, data=dat))
   sav <- influence(res)

   expect_equivalent(length(sav$inf$rstudent), 10)
   expect_equivalent(sum(is.na(sav$inf$rstudent)), 1)
   expect_equivalent(sum(is.na(sav$inf$hat)), 0)
   expect_equivalent(sum(is.na(sav$dfbs$intrcpt)), 1)

   options(na.action="na.pass")

   expect_warning(res <- rma(yi, vi, mods = ~ ablat + trial12, data=dat))
   sav <- influence(res)

   expect_equivalent(length(sav$inf$rstudent), 13)
   expect_equivalent(sum(is.na(sav$inf$rstudent)), 4)
   expect_equivalent(sum(is.na(sav$inf$hat)), 3)
   expect_equivalent(sum(is.na(sav$dfbs$intrcpt)), 4)

})
