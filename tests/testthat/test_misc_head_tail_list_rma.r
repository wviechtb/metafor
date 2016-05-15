### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: head.list.rma() and tail.list.rma() functions")

test_that("head.list.rma() works correctly.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   res <- head(rstandard(res), 4)
   res <- lapply(res, round, 4)
   class(res) <- "list.rma"

   sav <- structure(list(resid = c(-0.1748, -0.8709, -0.6335, -0.727), se = c(0.7788, 0.6896, 0.8344, 0.5486), z = c(-0.2244, -1.2629, -0.7592, -1.3253), slab = c(1, 2, 3, 4), digits = 4), .Names = c("resid", "se", "z", "slab", "digits"), class = "list.rma")

   expect_equivalent(res, sav)

})

test_that("tail.list.rma() works correctly.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat)
   res <- tail(rstandard(res), 4)
   res <- lapply(res, round, 4)
   class(res) <- "list.rma"

   sav <- structure(list(resid = c(-0.6568, 0.3752, 1.1604, 0.6972), se = c(0.5949, 0.5416, 0.9019, 0.5936), z = c(-1.104, 0.6927, 1.2867, 1.1746), slab = c(10, 11, 12, 13), digits = 4), .Names = c("resid", "se", "z", "slab", "digits"), class = "list.rma")

   expect_equivalent(res, sav)

})
