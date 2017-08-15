### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: predict() function")

test_that("predict() correctly matches named vectors in 'newmods'", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   dat$alloc[dat$alloc == "systematic"] <- "system"

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

   expect_error(predict(res, newmods = c(30, 0))) # not the right length
   expect_error(predict(res, newmods = c(abl = 30, alloc = 0, sys = 1))) # alloc matches up equally to allocrandom and allocsystem
   expect_error(predict(res, newmods = c(abl = 30, ran = 0, year = 1970))) # year not in the model
   expect_error(predict(res, newmods = c(abl = 30, ran = 0, sys = 1, ran = 1))) # ran used twice

   res <- rma(yi ~ ablat * year, vi, data=dat)
   pred1 <- predict(res, newmods = c(30, 1970, 30*1970))
   pred2 <- predict(res, newmods = c('ablat' = 30, 'year' = 1970, 'ablat:year' = 30*1970))
   pred3 <- predict(res, newmods = c('ablat:year' = 30*1970, 'year' = 1970, 'ablat' = 30))
   pred4 <- predict(res, newmods = c('ab' = 30, 'ye' = 1970, 'ablat:' = 30*1970))
   pred5 <- predict(res, newmods = c('ablat:' = 30*1970, 'ye' = 1970, 'ab' = 30))
   expect_equivalent(pred1, pred2)
   expect_equivalent(pred1, pred3)
   expect_equivalent(pred1, pred4)
   expect_equivalent(pred1, pred5)

})

test_that("predict() gives correct results when vcov=TRUE", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, data=dat)
   sav <- predict(res, vcov=TRUE)
   expect_equivalent(round(sav$pred$se, 4), round(c(sqrt(sav$vcov)), 4))

   res <- rma(yi, vi, mods = ~ ablat, data=dat)
   sav <- predict(res, vcov=TRUE)
   expect_equivalent(round(sav$pred$se, 4), round(c(sqrt(diag(sav$vcov))), 4))

}
