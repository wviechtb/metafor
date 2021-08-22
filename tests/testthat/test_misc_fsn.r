### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: fsn() function")

source("tolerances.r") # read in tolerances

test_that("confint() gives correct results for the 'expectancy data' in Becker (2005).", {

   data(dat.raudenbush1985)
   sav <- fsn(yi, vi, data=dat.raudenbush1985)

   expect_equivalent(sav$fsnum, 26)
   ### note: Becker uses p-values based on t-tests, which yields N =~ 23

   out <- capture.output(print(sav)) ### so that print.fsn() is run (at least once)

   sav <- fsn(yi, data=dat.raudenbush1985, type="Orwin", target=.05)
   expect_equivalent(sav$fsnum, 44)
   ### note: Becker finds N = 4, but uses the FE model estimate with 1/vi weights for
   ### the average effect size, but Orwin's methods is based on units weighting

   out <- capture.output(print(sav)) ### so that print.fsn() is run (at least once) with type="Orwin"

   sav <- fsn(yi, vi, data=dat.raudenbush1985, type="Rosenberg")
   expect_equivalent(sav$fsnum, 0)

   out <- capture.output(print(sav)) ### so that print.fsn() is run (at least once) with type="Rosenberg"

})

test_that("confint() gives correct results for the 'passive smoking data' in Becker (2005).", {

   data(dat.hackshaw1998)
   sav <- fsn(yi, vi, data=dat.hackshaw1998)

   expect_equivalent(sav$fsnum, 393)
   ### note: Becker finds N =~ 398 (due to rounding)

   sav <- fsn(yi, data=dat.hackshaw1998, type="Orwin", target=.049)
   expect_equivalent(sav$fsnum, 186)
   ### note: Becker finds N = 103, but uses the FE model estimate with 1/vi weights for
   ### the average effect size, but Orwin's methods is based on units weighting

   sav <- fsn(yi, vi, data=dat.hackshaw1998, type="Rosenberg")
   expect_equivalent(sav$fsnum, 202)

})

test_that("confint() gives correct results for the 'interview data' in Becker (2005).", {

   data(dat.mcdaniel1994)
   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat.mcdaniel1994)
   sav <- fsn(yi, vi, data=dat)

   expect_equivalent(sav$fsnum, 50364)
   ### note: Becker uses p-values based on t-tests, which yields N =~ 51226

   sav <- fsn(yi, data=dat, type="Orwin", target=.15)
   expect_equivalent(sav$fsnum, 129)
   ### note: Becker finds N = 64, but uses the FE model estimate with 1/vi weights for
   ### the average effect size, but Orwin's methods is based on units weighting

   sav <- fsn(yi, vi, data=dat, type="Rosenberg")
   expect_equivalent(sav$fsnum, 45528)

})
