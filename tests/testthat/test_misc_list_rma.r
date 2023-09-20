### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: head.list.rma() and tail.list.rma() functions")

source("settings.r")

test_that("head.list.rma() works correctly.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   res <- head(rstandard(res), 4)

   sav <- structure(list(resid = c(-0.1748, -0.8709, -0.6335, -0.727), se = c(0.7788, 0.6896, 0.8344, 0.5486), z = c(-0.2244, -1.2629, -0.7592, -1.3253), slab = 1:4, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), class = "list.rma")

   expect_equivalent(res, sav, tolerance=.tol[["misc"]])

})

test_that("tail.list.rma() works correctly.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   res <- tail(rstandard(res), 4)

   sav <- structure(list(resid = c(-0.6568, 0.3752, 1.1604, 0.6972), se = c(0.5949, 0.5416, 0.9019, 0.5936), z = c(-1.104, 0.6927, 1.2867, 1.1746 ), slab = 10:13, digits = c(est = 4, se = 4, test = 4, pval = 4, ci = 4, var = 4, sevar = 4, fit = 4, het = 4)), class = "list.rma")

   expect_equivalent(res, sav, tolerance=.tol[["misc"]])

})

test_that("as.data.frame.list.rma() works correctly.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)
   res <- predict(res)
   res <- as.data.frame(res)
   res <- res[1:3,1:2]

   sav <- structure(list(pred = c(-1.02900878645837, -1.34912705666653, -0.97080546460234), se = c(0.140375124151501, 0.201103941277043, 0.131456743392091)), row.names = c(NA, 3L), class = "data.frame")

   expect_equivalent(res, sav, tolerance=.tol[["misc"]])

})

test_that("as.matrix.list.rma() works correctly.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)
   res <- predict(res)
   res <- as.matrix(res)
   res <- res[1:3,1:2]

   sav <- structure(c(-1.02900878645837, -1.34912705666653, -0.97080546460234, 0.140375124151501, 0.201103941277043, 0.131456743392091), dim = 3:2, dimnames = list(c("1", "2", "3"), c("pred", "se")))

   expect_equivalent(res, sav, tolerance=.tol[["misc"]])

})

rm(list=ls())
