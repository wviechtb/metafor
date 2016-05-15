### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: confint() function")

test_that("confint() works correctly for rma().", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, data=dat, method="DL")
   sav <- confint(res, fixed=TRUE, transf=exp)

   expect_equivalent(round(c(sav$fixed), 4), c(0.4896, 0.3449, 0.6950))
   expect_equivalent(round(sav$random[1,], 4), c(0.3088, 0.1197, 1.1115))
   expect_equivalent(round(sav$random[3,], 4), c(92.1173, 81.9177, 97.6781))
   expect_equivalent(round(sav$random[4,], 4), c(12.6861, 5.5303, 43.0680))

})

test_that("confint() works correctly for rma.mh().", {

   data(dat.bcg, package="metafor")
   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   sav <- confint(res, transf=exp)

   expect_equivalent(round(c(sav$fixed), 4), c(0.6353, 0.5881, 0.6862))

})

test_that("confint() works correctly for rma.peto().", {

   data(dat.bcg, package="metafor")
   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   sav <- confint(res, transf=exp)

   expect_equivalent(round(c(sav$fixed), 4), c(0.6222, 0.5746, 0.6738))

})
