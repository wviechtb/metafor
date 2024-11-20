### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: predict() function")

source("settings.r")

test_that("predict() correctly matches named vectors in 'newmods'", {

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
   expect_error(predict(res, newmods = c(abl = 30, random = 0))) # not the right length
   expect_error(predict(res, newmods = c(abl = 30, alloc = 0, sys = 1))) # alloc matches up equally to allocrandom and allocsystem
   expect_error(predict(res, newmods = c(abl = 30, ran = 0, year = 1970))) # year not in the model
   expect_error(predict(res, newmods = c(abl = 30, ran = 0, sys = 1, ran = 1))) # ran used twice
   expect_error(predict(res, newmods = c(abl = 30, ran = 0, sys = 1, rand = 1))) # same issue

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

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, data=dat)
   sav <- predict(res, vcov=TRUE)
   expect_equivalent(sav$pred$se, c(sqrt(sav$vcov)), tolerance=.tol[["se"]])

   res <- rma(yi, vi, mods = ~ ablat, data=dat)
   sav <- predict(res, vcov=TRUE)
   expect_equivalent(sav$pred$se, c(sqrt(diag(sav$vcov))), tolerance=.tol[["se"]])

})

test_that("predict() correctly handles in/exclusion of the intercept term", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   #########################################################################

   # single quantitative predictor model with intercept included
   res <- rma(yi ~ ablat, vi, data=dat)

   # predicted average effect at ablat=0,10,...,60
   pred1 <- predict(res, newmods=seq(0,60,by=10))
   pred2 <- predict(res, newmods=cbind(1,seq(0,60,by=10)))
   expect_equivalent(pred1, pred2)

   # exclude the intercept from the prediction (i.e., assume it is 0)
   pred1 <- predict(res, newmods=seq(0,60,by=10), intercept=FALSE)
   pred2 <- predict(res, newmods=cbind(0,seq(0,60,by=10)))
   expect_equivalent(pred1, pred2)

   expect_warning(pred2 <- predict(res, newmods=cbind(0,seq(0,60,by=10)), intercept=FALSE))

   #########################################################################

   # single quantitative predictor model with intercept excluded
   res <- rma(yi ~ 0 + ablat, vi, data=dat)

   # predicted average effect at ablat=0,10,...,60
   pred1 <- predict(res, newmods=seq(0,60,by=10))
   pred2 <- predict(res, newmods=cbind(seq(0,60,by=10)))
   expect_equivalent(pred1, pred2)

   #########################################################################

   # multiple predictors one of which is categorical with intercept included/excluded
   res1 <- rma(yi ~ 1 + ablat + alloc, vi, data=dat)
   res0 <- rma(yi ~ 0 + ablat + alloc, vi, data=dat)

   # predicted average effect at ablat=20 for alloc='random'
   pred1 <- predict(res1, newmods=c(20,1,0))
   pred0 <- predict(res0, newmods=c(20,0,1,0))
   expect_equivalent(pred1, pred0)
   pred2 <- predict(res1, newmods=cbind(1,20,1,0))
   expect_equivalent(pred1, pred2)
   pred2 <- predict(res0, newmods=cbind(20,0,1,0))
   expect_equivalent(pred1, pred2)

   pred1 <- predict(res1, newmods=cbind(20,1,0))
   pred0 <- predict(res0, newmods=cbind(20,0,1,0))
   expect_equivalent(pred1, pred0)

   # predicted average effect at ablat=0,10,...,60 for alloc='random'
   pred1 <- predict(res1, newmods=cbind(seq(0,60,by=10),1,0))
   pred0 <- predict(res0, newmods=cbind(seq(0,60,by=10),0,1,0))
   expect_equivalent(pred1, pred0)
   pred2 <- predict(res1, newmods=cbind(1,seq(0,60,by=10),1,0))
   expect_equivalent(pred1, pred2)

   # contrast between alloc='random' and alloc='systematic' holding ablat constant
   pred1 <- predict(res1, newmods=c(0,1,-1), intercept=FALSE)
   pred0 <- predict(res0, newmods=c(0,0,1,-1))
   expect_equivalent(pred1, pred0)
   pred2 <- predict(res1, newmods=cbind(0,0,1,-1))
   expect_equivalent(pred1, pred2)
   pred2 <- predict(res0, newmods=cbind(0,0,1,-1))
   expect_equivalent(pred1, pred2)

   #########################################################################

})

test_that("predict() works correctly with adjusted level", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi ~ ablat, vi, data=dat)
   pred1 <- predict(res, newmods=seq(0,60,by=10), level=90)
   res <- rma(yi ~ ablat, vi, data=dat, level=90)
   pred2 <- predict(res, newmods=seq(0,60,by=10))
   expect_equivalent(pred1, pred2)

   res <- rma(yi ~ ablat, vi, data=dat)
   res <- robust(res, cluster=trial)
   pred1 <- predict(res, newmods=seq(0,60,by=10), level=90)
   res <- rma(yi ~ ablat, vi, data=dat, level=90)
   res <- robust(res, cluster=trial)
   pred2 <- predict(res, newmods=seq(0,60,by=10))
   expect_equivalent(pred1, pred2)

   res <- rma(yi ~ ablat, vi, data=dat)
   res <- robust(res, cluster=trial, clubSandwich=TRUE)
   pred1 <- predict(res, newmods=seq(0,60,by=10), level=90)
   res <- rma(yi ~ ablat, vi, data=dat, level=90)
   res <- robust(res, cluster=trial, clubSandwich=TRUE)
   pred2 <- predict(res, newmods=seq(0,60,by=10))
   expect_equivalent(pred1, pred2)

})

rm(list=ls())
