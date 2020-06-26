### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: vif() function")

source("tolerances.r") # read in tolerances

test_that("vif() works correctly for 'rma.uni' objects.", {

   dat <- dat.bangertdrowns2004
   dat <- dat[!apply(dat[,c("length", "wic", "feedback", "info", "pers", "imag", "meta")], 1, anyNA),]
   res <- rma(yi, vi, mods = ~ length + wic + feedback + info + pers + imag + meta, data=dat)

   sav <- vif(res)
   out <- capture.output(print(sav))

   vifs <- c(length = 1.53710262575577, wic = 1.38604929927746, feedback = 1.64904565071108, info = 1.83396138431786, pers = 5.67803138275492, imag = 1.1553714953831, meta = 4.53327503733189)
   expect_equivalent(sav$vif, vifs)

   sav <- vif(res, table=TRUE)
   out <- capture.output(print(sav))

   expect_equivalent(sav$vif$vif[-1], vifs)

   sav <- vif(res, btt=2:3)
   out <- capture.output(print(sav))

   gvif <- 2.06507966959426
   expect_equivalent(sav$gvif, gvif)

})
