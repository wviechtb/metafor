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

test_that("escalc() works correctly for measure='SMD/SMDH/ROM'", {

   dat <- get(data(dat.normand1999, package="metafor"))

   sav <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(c(round(sav$yi, 4)), c(-0.3552, -0.3479, -2.3176, -1.8880))
   expect_equivalent(c(round(sav$vi, 4)), c( 0.0131,  0.0645,  0.0458,  0.1606))

   sav <- escalc(measure="SMDH", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(c(round(sav$yi, 4)), c(-0.3553, -0.3465, -2.3018, -1.8880))
   expect_equivalent(c(round(sav$vi, 4)), c( 0.0131,  0.0657,  0.0509,  0.1874))

   sav <- escalc(measure="ROM", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(c(round(sav$yi, 4)), c(-0.3102, -0.0715, -0.6202, -0.7303))
   expect_equivalent(c(round(sav$vi, 4)), c( 0.0094,  0.0028,  0.0018,  0.0119))

})

test_that("escalc() works correctly for measure='RPB/RBIS'", {

   x <- c(20, 31, 18, 22, 30, 16, 28, 24, 23, 27,  1,  4,  8, 15,  9, 11, 11,  6,  8,  4)
   y <- c(3,   3,  4,  5,  6,  4,  7,  6,  5,  4,  3,  5,  1,  5,  2,  4,  6,  4,  2,  4)
   xb <- ifelse(x > median(x), 1, 0)

   sav <- escalc(measure="RPB", m1i=mean(y[xb==1]), sd1i=sd(y[xb==1]), n1i=sum(xb==1), m2i=mean(y[xb==0]), sd2i=sd(y[xb==0]), n2i=sum(xb==0))
   expect_equivalent(c(round(sav$yi, 4)), 0.3685)
   expect_equivalent(c(round(sav$vi, 4)), 0.0384)

   sav <- escalc(measure="RBIS", m1i=mean(y[xb==1]), sd1i=sd(y[xb==1]), n1i=sum(xb==1), m2i=mean(y[xb==0]), sd2i=sd(y[xb==0]), n2i=sum(xb==0))
   expect_equivalent(c(round(sav$yi, 4)), 0.4619)
   expect_equivalent(c(round(sav$vi, 4)), 0.0570)

})

test_that("escalc() works correctly for measure='D2ORL/D2ORN'", {

   dat <- get(data(dat.gibson2002, package="metafor"))

   sav <- escalc(measure="D2ORL", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(c(round(sav$yi, 4)), c(-0.4315, -0.9285,  0.5932, -0.1890))
   expect_equivalent(c(round(sav$vi, 4)), c( 0.1276,  0.0493,  0.3204,  0.0690))

   sav <- escalc(measure="D2ORN", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(c(round(sav$yi, 4)), c(-0.3925, -0.8447,  0.5397, -0.1719))
   expect_equivalent(c(round(sav$vi, 4)), c( 0.1056,  0.0408,  0.2651,  0.0571))

})

test_that("escalc() works correctly for measure='COR/UCOR/ZCOR'", {

   dat <- get(data(dat.mcdaniel1994, package="metafor"))

   sav <- escalc(measure="COR", ri=ri, ni=ni, data=dat, subset=c(1,13,33,102))
   expect_equivalent(c(round(sav$yi, 4)), c(0.0000, 0.6200, 0.9900, -0.1300))
   expect_equivalent(c(round(sav$vi, 4)), c(0.0082, 0.0271, 0.0001,  0.0242))

   sav <- escalc(measure="UCOR", ri=ri, ni=ni, data=dat, subset=c(1,13,33,102))
   expect_equivalent(c(round(sav$yi, 4)), c(0.0000, 0.6363, 0.9925, -0.1317))
   expect_equivalent(c(round(sav$vi, 4)), c(0.0082, 0.0253, 0.0000,  0.0241))

   sav <- escalc(measure="UCOR", ri=ri, ni=ni, data=dat, vtype="UB", subset=c(1,13,33,102))
   expect_equivalent(c(round(sav$vi, 4)), c(0.0084, 0.0283, 0.0000,  0.0261))

   sav <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat, subset=c(1,13,33,102))
   expect_equivalent(c(round(sav$yi, 4)), c(0.0000, 0.7250, 2.6467, -0.1307))
   expect_equivalent(c(round(sav$vi, 4)), c(0.0083, 0.0833, 0.3333,  0.0263))

})
