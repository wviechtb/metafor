### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: predict() function")

test_that("predict() correctly matches named vectors in 'newmods'", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi ~ ablat + alloc, vi, data=dat)
   pred1 <- predict(res, newmods = c(30, 0, 1))
   pred2 <- predict(res, newmods = c(abl = 30, ran = 0, sys = 1))
   pred3 <- predict(res, newmods = c(abl = 30, sys = 1, ran = 0))
   pred4 <- predict(res, newmods = c(ran = 0, abl = 30, sys = 1))
   pred5 <- predict(res, newmods = c(sys = 1, abl = 30, ran = 0))
   pred6 <- predict(res, newmods = c(ran = 0, sys = 1, abl = 30))
   pred7 <- predict(res, newmods = c(sys = 1, ran = 0, abl = 30))
   expect_equivalent(pred1, pred2)
   expect_equivalent(pred1, pred3)
   expect_equivalent(pred1, pred4)
   expect_equivalent(pred1, pred5)
   expect_equivalent(pred1, pred6)
   expect_equivalent(pred1, pred7)

})
