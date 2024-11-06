### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: fsn() function")

source("settings.r")

test_that("confint() gives correct results for the 'expectancy data' in Becker (2005).", {

   sav <- fsn(yi, vi, data=dat.raudenbush1985)

   expect_equivalent(sav$fsnum, 26)
   ### note: Becker uses p-values based on t-tests, which yields N =~ 23

   out <- capture.output(print(sav)) # so that print.fsn() is run (at least once)

   ### use Fisher's test
   sav <- fsn(yi, vi, data=dat.raudenbush1985, pool="Fisher")
   expect_equivalent(sav$fsnum, 40)

   sav <- fsn(yi, data=dat.raudenbush1985, type="Orwin", target=.05)
   expect_equivalent(sav$fsnum, 44)

   out <- capture.output(print(sav)) # so that print.fsn() is run (at least once) with type="Orwin"

   sav <- fsn(yi, vi, data=dat.raudenbush1985, type="Orwin", target=.05)
   expect_equivalent(sav$fsnum, 4)

   sav <- fsn(yi, vi, data=dat.raudenbush1985, type="Rosenberg")
   expect_equivalent(sav$fsnum, 0)

   out <- capture.output(print(sav)) # so that print.fsn() is run (at least once) with type="Rosenberg"

   skip_on_cran()

   sav <- fsn(yi, vi, data=dat.raudenbush1985, type="General")
   expect_equivalent(sav$fsnum, 0)

   sav <- fsn(yi, vi, data=dat.raudenbush1985, type="General", exact=TRUE)
   expect_equivalent(sav$fsnum, 0)

   out <- capture.output(print(sav)) # so that print.fsn() is run (at least once) with type="General"

   res <- rma(yi, vi, data=dat.raudenbush1985)
   sav <- fsn(res, target=.05)
   expect_equivalent(sav$fsnum, 12)

})

test_that("confint() gives correct results for the 'passive smoking data' in Becker (2005).", {

   sav <- fsn(yi, vi, data=dat.hackshaw1998)

   expect_equivalent(sav$fsnum, 393)
   ### note: Becker finds N =~ 398 (due to rounding)

   sav <- fsn(yi, data=dat.hackshaw1998, type="Orwin", target=.049)
   expect_equivalent(sav$fsnum, 186)

   sav <- fsn(yi, vi, data=dat.hackshaw1998, type="Orwin", target=.049)
   expect_equivalent(sav$fsnum, 104) # not 103 as fsn() always rounds up

   sav <- fsn(yi, vi, data=dat.hackshaw1998, type="Rosenberg")
   expect_equivalent(sav$fsnum, 202)

   skip_on_cran()

   sav <- fsn(yi, vi, data=dat.hackshaw1998, type="General")
   expect_equivalent(sav$fsnum, 112)

   sav <- fsn(yi, vi, data=dat.hackshaw1998, type="General", exact=TRUE)
   expect_equivalent(sav$fsnum, 119)

})

test_that("confint() gives correct results for the 'interview data' in Becker (2005).", {

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.mcdaniel1994)

   sav <- fsn(yi, vi, data=dat)

   expect_equivalent(sav$fsnum, 50364)
   ### note: Becker uses p-values based on t-tests, which yields N =~ 51226

   sav <- fsn(yi, data=dat, type="Orwin", target=.15)
   expect_equivalent(sav$fsnum, 129)

   sav <- fsn(yi, vi, data=dat, type="Orwin", target=.15)
   expect_equivalent(sav$fsnum, 65) # not 64 as fsn() always rounds up

   sav <- fsn(yi, vi, data=dat, type="Rosenberg")
   expect_equivalent(sav$fsnum, 45528)

   skip_on_cran()

   sav <- fsn(yi, vi, data=dat, type="General")
   expect_equivalent(sav$fsnum, 6068)

   sav <- fsn(yi, vi, data=dat, type="General", exact=TRUE)
   expect_equivalent(sav$fsnum, 6068)

   res <- rma(yi, vi, data=dat)

   sav <- fsn(res)
   expect_equivalent(sav$fsnum, 6068)

})

rm(list=ls())
