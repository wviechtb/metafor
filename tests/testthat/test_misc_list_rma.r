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

test_that("as.data.frame.list.rma() works correctly.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)
   res <- predict(res)
   res <- as.data.frame(res)
   res <- round(res[1:3,1:2], 4)

   sav <- structure(list(pred = c(-1.029, -1.3491, -0.9708), se = c(0.1404, 0.2011, 0.1315)), .Names = c("pred", "se"), row.names = c(NA, 3L), class = "data.frame")

   expect_equivalent(res, sav)

})

test_that("as.matrix.list.rma() works correctly.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat, data=dat)
   res <- predict(res)
   res <- as.matrix(res)
   res <- round(res[1:3,1:2], 4)

   sav <- structure(c(-1.029, -1.3491, -0.9708, 0.1404, 0.2011, 0.1315), .Dim = c(3L, 2L), .Dimnames = list(c("1", "2", "3"), c("pred", "se")))

   expect_equivalent(res, sav)

})
