### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking tip: rma() results match up with those from lm()")

source("settings.r")

test_that("results for rma() and lm() match.", {

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)

   res1 <- rma(yi, 0, data=dat)
   res2 <- lm(yi ~ 1, data=dat)

   ### coefficients should be the same
   expect_equivalent(coef(res1), coef(res2))

   ### standard errors should be the same
   expect_equivalent(se(res1), se(res2))

})

test_that("results for rma.mv() and lm() match.", {

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.molloy2014)
   dat$id <- 1:nrow(dat)

   res1 <- rma.mv(yi, 0, random = ~ 1 | id, data=dat, sparse=.sparse)
   res2 <- lm(yi ~ 1, data=dat)

   ### coefficients should be the same
   expect_equivalent(coef(res1), coef(res2))

   ### standard errors should be the same
   expect_equivalent(se(res1), se(res2))

   ### get profile likelihood CI for sigma^2
   sav <- confint(res1)
   expect_equivalent(sav$random[1,2:3], c(.0111, .0474), tolerance=.tol[["var"]])

   ### fit with sparse=TRUE
   res1 <- rma.mv(yi, 0, random = ~ 1 | id, data=dat, sparse=TRUE)

   ### coefficients should be the same
   expect_equivalent(coef(res1), coef(res2))

   ### standard errors should be the same
   expect_equivalent(se(res1), se(res2))

   ### get profile likelihood CI for sigma^2
   sav <- confint(res1)
   expect_equivalent(sav$random[1,2:3], c(.0111, .0474), tolerance=.tol[["var"]])

})

rm(list=ls())
