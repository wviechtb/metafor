### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking misc: influence() and related functions")

source("settings.r")

test_that("influence() works for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   sav <- influence(res)
   sav$inf  <- sav$inf[1]
   sav$dfbs <- sav$dfbs[1]
   sav$is.infl <- sav$is.infl[1]
   sav$not.na <- sav$not.na[1]

   tmp <- structure(list(inf = structure(list(rstudent = -0.218142474344442, dffits = -0.0407075604868486, cook.d = 0.00171654236729195, cov.r = 1.11644891104804, tau2.del = 0.336156745300306, QE.del = 151.582572747109, hat = 0.0505948307931551, weight = 5.05948307931551, inf = "", slab = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), class = "list.rma"), dfbs = structure(list(intrcpt = -0.0402659025974144, slab = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), class = "list.rma"), ids = 1:13, not.na = TRUE, is.infl = FALSE, tau2 = 0.313243325980895, QE = 152.233008082373, k = 13L, p = 1L, m = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), class = "infl.rma.uni")

   expect_equivalent(sav, tmp, tolerance=.tol[["inf"]])

})

test_that("leave1out() works for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   inf <- leave1out(res)
   inf <- inf[1]

   sav <- structure(list(estimate = -0.707083788031436, se = 0.189961024702717, zval = -3.72225717953459, pval = 0.000197449759023198, ci.lb = -1.07940055491509, ci.ub = -0.334767021147788, Q = 151.582572747109, Qp = 7.0778599767807e-27, tau2 = 0.336156745300306, I2 = 93.2259349111223, H2 = 14.762184698253, slab = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4), transf = FALSE), class = "list.rma")

   expect_equivalent(sav, inf, tolerance=.tol[["misc"]])

})

test_that("leave1out() works for rma.mh().", {

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   inf <- leave1out(res)
   inf <- inf[1]

   sav <- structure(list(estimate = -0.451379469928476, se = 0.0394350331703394, zval = -11.4461541842439, pval = 2.45810944109134e-30, ci.lb = -0.528670714671484, ci.ub = -0.374088225185468, Q = 151.915260738878, Qp = 6.05181927235005e-27, I2 = 92.7591211399706, H2 = 13.8104782489889, slab = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4), transf = FALSE), class = "list.rma")

   expect_equivalent(sav, inf, tolerance=.tol[["misc"]])

})

test_that("leave1out() works for rma.peto().", {

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   inf <- leave1out(res)
   inf <- inf[1]

   sav <- structure(list(estimate = -0.472177269248539, se = 0.0407784291562603, zval = -11.5790941195696, pval = 5.25989306490064e-31, ci.lb = -0.552101521740927, ci.ub = -0.39225301675615, Q = 167.200450619361, Qp = 4.44309617192221e-30, I2 = 93.4210703623987, H2 = 15.2000409653964, slab = 1L, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4), transf = FALSE), class = "list.rma")

   expect_equivalent(sav, inf, tolerance=.tol[["misc"]])

})

test_that("model.matrix() works for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   sav <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 44, 55, 42, 52, 13, 44, 19, 13, 27, 42, 18, 33, 33), .Dim = c(13L, 2L), .Dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"), c("intrcpt", "ablat")))

   expect_equivalent(sav, model.matrix(res))

})

test_that("hatvalues() works for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   expect_equivalent(hatvalues(res), c(0.049, 0.1493, 0.0351, 0.3481, 0.2248, 0.2367, 0.064, 0.357, 0.0926, 0.1157, 0.2309, 0.0189, 0.0778), tolerance=.tol[["inf"]])

   sav <- structure(c(0.049, 0.067, 0.0458, 0.0994, 0.1493, 0.0904, 0.0374, 0.0498, 0.0351), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"), c("1", "2", "3")))
   expect_equivalent(hatvalues(res, type="matrix")[1:3,1:3], sav, tolerance=.tol[["inf"]])

})

test_that("hatvalues() works for rma.mv().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, mods = ~ ablat, random = ~ 1 | trial, data=dat, sparse=.sparse)

   expect_equivalent(hatvalues(res), c(0.049, 0.1493, 0.0351, 0.3481, 0.2248, 0.2367, 0.064, 0.357, 0.0926, 0.1157, 0.2309, 0.0189, 0.0778), tolerance=.tol[["inf"]])

   sav <- structure(c(0.049, 0.067, 0.0458, 0.0994, 0.1493, 0.0904, 0.0374, 0.0498, 0.0351), .Dim = c(3L, 3L), .Dimnames = list(c("1", "2", "3"), c("1", "2", "3")))
   expect_equivalent(hatvalues(res, type="matrix")[1:3,1:3], sav, tolerance=.tol[["inf"]])

})

test_that("cooks.distance() works for rma().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)

   expect_equivalent(cooks.distance(res), c(0.0048, 0.0489, 0.0104, 0.2495, 0.0072, 0.2883, 0.3643, 0.2719, 0.02, 0.1645, 0.0009, 0.0403, 0.1433), tolerance=.tol[["inf"]])

})

test_that("cooks.distance() works for rma.mv().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma.mv(yi, vi, mods = ~ ablat, random = ~ 1 | trial, data=dat, sparse=.sparse)

   expect_equivalent(cooks.distance(res), c(0.0048, 0.0489, 0.0104, 0.2495, 0.0072, 0.2883, 0.3643, 0.2719, 0.02, 0.1645, 0.0009, 0.0404, 0.1434), tolerance=.tol[["inf"]])
   expect_equivalent(cooks.distance(res, cluster=alloc), c(0.2591, 2.4372, 0.1533), tolerance=.tol[["inf"]])
   expect_equivalent(cooks.distance(res, cluster=alloc, reestimate=FALSE), c(0.3199, 2.2194, 0.2421), tolerance=.tol[["inf"]])

})

test_that("influence() correctly works with 'na.omit' and 'na.pass'.", {

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

   options(na.action="na.omit")

})

test_that("'infonly' argument works correctly with influence().", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, slab=paste0("Trial ", dat.bcg$trial))
   res <- rma(yi, vi, data=dat, method="EE")
   inf <- influence(res)
   tmp <- capture.output(sav <- print(inf))
   expect_equivalent(length(sav$rstudent), 13)
   tmp <- capture.output(sav <- print(inf, infonly=TRUE))
   expect_equivalent(length(sav$rstudent), 3)

})

rm(list=ls())
