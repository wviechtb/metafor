### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: escalc() function")

test_that("escalc() works correctly for measure='RR'", {

   data(dat.bcg, package="metafor")

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(round(dat$yi[1],4), -0.8893)
   expect_equivalent(round(dat$vi[1],4),  0.3256)

   ### need to add lots of stuff here ...

})

test_that("escalc.formula() works correctly for measure='RR'", {

   data(dat.bcg, package="metafor")

   dat.long <- to.long(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
                       data=dat.bcg, append=FALSE, vlong=TRUE)

   dat <- escalc(measure="RR", outcome ~ group | study, weights=freq, data=dat.long)

   expect_equivalent(round(dat$yi[1],4), -0.8893)
   expect_equivalent(round(dat$vi[1],4),  0.3256)

})

test_that("escalc() works correctly for measure='PHI/YUQ/YUY/RTET/PBIT/OR2D/OR2DN'", {

   ### see Table 13.4 (p. 242) in the Handbook of Research Synthesis and Meta-Analysis

   dat <- escalc(measure="PHI", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(round(sav$yi, 3), 0.131)
   expect_equivalent(round(sav$sei, 3), 0.079)

   dat <- escalc(measure="YUQ", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(round(sav$yi, 3), 0.385)
   expect_equivalent(round(sav$sei, 3), 0.190)

   dat <- escalc(measure="YUY", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(round(sav$yi, 3), 0.200)
   expect_equivalent(round(sav$sei, 3), 0.107)

   dat <- escalc(measure="RTET", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(round(sav$yi, 3), 0.260)
   expect_equivalent(round(sav$sei, 3), 0.142)

   dat <- escalc(measure="PBIT", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(round(sav$yi, 3), 0.440)
   expect_equivalent(round(sav$sei, 3), 0.246)

   dat <- escalc(measure="OR2D", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(round(sav$yi, 3), 0.447)
   expect_equivalent(round(sav$sei, 3), 0.246)

   dat <- escalc(measure="OR2DN", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(round(sav$yi, 3), 0.491)
   expect_equivalent(round(sav$sei, 3), 0.270)

})
