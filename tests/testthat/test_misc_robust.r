### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: robust() function")

source("settings.r")

test_that("robust() works correctly for 'rma' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, data=dat)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.032106, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.98776, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, adjust=FALSE)
   expect_equivalent(c(vcov(sav)), 0.029636, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -4.150592, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, clubSandwich=TRUE)
   expect_equivalent(c(vcov(sav)), 0.03229357, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 11.04125, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.97616, tolerance=.tol[["test"]])

   res <- rma(yi, vi, weights=1, data=dat)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.037028, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.848996, tolerance=.tol[["test"]])

})

test_that("robust() works correctly for 'rma' objects with moderators.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, mods = ~ ablat + year, data=dat)

   sav <- robust(res, cluster=trial)
   expect_equivalent(sav$se, c(23.910483, 0.007857, 0.012079), tolerance=.tol[["se"]])
   expect_equivalent(sav$dfs, 10, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, c(-0.148282, -3.564978, 0.157928), tolerance=.tol[["test"]])
   expect_equivalent(sav$QM, 11.8546, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, clubSandwich=TRUE)
   expect_equivalent(sav$se, c(33.655367, 0.011994, 0.016963), tolerance=.tol[["se"]])
   expect_equivalent(sav$dfs, c(2.724625, 2.112895, 2.745919), tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, c(-0.105347, -2.335398, 0.112456), tolerance=.tol[["test"]])
   expect_equivalent(sav$QM, 6.708996, tolerance=.tol[["test"]])
   expect_equivalent(sav$QMdf, c(2, 2.528214), tolerance=.tol[["misc"]])
   expect_equivalent(sav$QMp, 0.097479, tolerance=.tol[["pval"]])

   tmp <- anova(sav)
   expect_equivalent(tmp$QM, 6.708996, tolerance=.tol[["test"]])
   expect_equivalent(tmp$QMdf, c(2, 2.528214), tolerance=.tol[["misc"]])
   expect_equivalent(tmp$QMp, 0.097479, tolerance=.tol[["pval"]])

   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)
   sav <- robust(res, cluster=trial, clubSandwich=TRUE)
   tmp <- anova(sav, X=rbind(c(0,10,1,0),c(0,50,1,0)))
   expect_equivalent(tmp$se, c(0.210162, 0.321173), tolerance=.tol[["se"]])
   expect_equivalent(tmp$ddf, c(1.929902, 3.251262), tolerance=.tol[["misc"]])
   expect_equivalent(tmp$zval, c(-2.570637, -5.079127), tolerance=.tol[["test"]])
   expect_equivalent(tmp$QM, 9.914783, tolerance=.tol[["test"]])
   expect_equivalent(tmp$QMdf, c(2, 2.569003), tolerance=.tol[["misc"]])
   expect_equivalent(tmp$QMp, 0.06194173, tolerance=.tol[["pval"]])

   sav1 <- robust(res, cluster=trial)
   tmp1 <- anova(sav1, X=rbind(c(0,10,1,0),c(0,50,1,0)))
   sav2 <- robust(res, cluster=trial, clubSandwich=TRUE, vcov="CR1p", coef_test="naive-tp", wald_test="Naive-Fp")
   tmp2 <- anova(sav, X=rbind(c(0,10,1,0),c(0,50,1,0)))
   expect_equivalent(tmp1$se, tmp2$se, tolerance=.tol[["se"]])
   expect_equivalent(tmp1$ddf, tmp2$ddf, tolerance=.tol[["misc"]])
   expect_equivalent(tmp1$zval, tmp2$zval, tolerance=.tol[["test"]])
   expect_equivalent(tmp1$QM, tmp2$QM, tolerance=.tol[["test"]])
   expect_equivalent(tmp1$QMdf, tmp2$QMdf, tolerance=.tol[["misc"]])
   expect_equivalent(tmp1$QMp, tmp2$QMp, tolerance=.tol[["pval"]])

})
test_that("robust() works correctly for 'rma.mv' objects.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, sparse=sparse)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.032106, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.98776, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, adjust=FALSE)
   expect_equivalent(c(vcov(sav)), 0.029636, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -4.150592, tolerance=.tol[["test"]])

   sav <- robust(res, cluster=trial, clubSandwich=TRUE)
   expect_equivalent(c(vcov(sav)), 0.03229357, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 11.04125, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.97616, tolerance=.tol[["test"]])

   res <- rma.mv(yi, vi, W=1, random = ~ 1 | trial, data=dat, sparse=sparse)

   sav <- robust(res, cluster=trial)
   expect_equivalent(c(vcov(sav)), 0.037028, tolerance=.tol[["var"]])
   expect_equivalent(sav$dfs, 12, tolerance=.tol[["misc"]])
   expect_equivalent(sav$zval, -3.848996, tolerance=.tol[["test"]])

})

rm(list=ls())
