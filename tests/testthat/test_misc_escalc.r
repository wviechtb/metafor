### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: escalc() function")

source("settings.r")

test_that("escalc() works correctly for measure='RR'", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(dat$yi[1], -0.8893, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi[1],  0.3256, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='PHI/YUQ/YUY/RTET/PBIT/OR2D/OR2DN'", {

   ### see Table 13.4 (p. 242) in the Handbook of Research Synthesis and Meta-Analysis

   dat <- escalc(measure="PHI", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(sav$yi, 0.1309, tolerance=.tol[["est"]])
   expect_equivalent(sav$sei, 0.0789, tolerance=.tol[["var"]])

   dat <- escalc(measure="YUQ", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(sav$yi, 0.3846, tolerance=.tol[["est"]])
   expect_equivalent(sav$sei, 0.1901, tolerance=.tol[["var"]])

   dat <- escalc(measure="YUY", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(sav$yi, 0.2000, tolerance=.tol[["est"]])
   expect_equivalent(sav$sei, 0.1071, tolerance=.tol[["var"]])

   dat <- escalc(measure="RTET", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(sav$yi, 0.2603, tolerance=.tol[["est"]])
   expect_equivalent(sav$sei, 0.1423, tolerance=.tol[["var"]])

   dat <- escalc(measure="PBIT", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(sav$yi, 0.4399, tolerance=.tol[["est"]])
   expect_equivalent(sav$sei, 0.2456, tolerance=.tol[["var"]])

   dat <- escalc(measure="OR2D", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(sav$yi, 0.4471, tolerance=.tol[["est"]])
   expect_equivalent(sav$sei, 0.2460, tolerance=.tol[["var"]])

   dat <- escalc(measure="OR2DN", ai=135, bi=15, ci=40, di=10)
   sav <- summary(dat)
   expect_equivalent(sav$yi, 0.4915, tolerance=.tol[["est"]])
   expect_equivalent(sav$sei, 0.2704, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='SMD/SMDH/ROM'", {

   dat <- dat.normand1999

   sav <- escalc(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(sav$yi, c(-0.3552, -0.3479, -2.3176, -1.8880), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c( 0.0131,  0.0645,  0.0458,  0.1606), tolerance=.tol[["var"]])

   sav <- escalc(measure="SMDH", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(sav$yi, c(-0.3553, -0.3465, -2.3018, -1.8880), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c( 0.0131,  0.0657,  0.0509,  0.1874), tolerance=.tol[["var"]])

   sav <- escalc(measure="ROM", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(sav$yi, c(-0.3102, -0.0715, -0.6202, -0.7303), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c( 0.0094,  0.0028,  0.0018,  0.0119), tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='CVR/VR'", {

   dat <- dat.normand1999
   dat <- escalc(measure="CVR", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1)

   expect_equivalent(dat$yi[1], 0.0014, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi[1], 0.0159, tolerance=.tol[["var"]])

   dat <- dat.normand1999
   dat <- escalc(measure="VR", sd1i=sd1i, n1i=n1i, sd2i=sd2i, n2i=n2i, data=dat, subset=1)

   expect_equivalent(dat$yi[1], -0.3087, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi[1], 0.0065, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='RPB/RBIS'", {

   x <- c(20, 31, 18, 22, 30, 16, 28, 24, 23, 27,  1,  4,  8, 15,  9, 11, 11,  6,  8,  4)
   y <- c(3,   3,  4,  5,  6,  4,  7,  6,  5,  4,  3,  5,  1,  5,  2,  4,  6,  4,  2,  4)
   xb <- ifelse(x > median(x), 1, 0)

   sav <- escalc(measure="RPB", m1i=mean(y[xb==1]), sd1i=sd(y[xb==1]), n1i=sum(xb==1), m2i=mean(y[xb==0]), sd2i=sd(y[xb==0]), n2i=sum(xb==0))
   expect_equivalent(sav$yi, 0.3685, tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, 0.0384, tolerance=.tol[["var"]])

   sav <- escalc(measure="RBIS", m1i=mean(y[xb==1]), sd1i=sd(y[xb==1]), n1i=sum(xb==1), m2i=mean(y[xb==0]), sd2i=sd(y[xb==0]), n2i=sum(xb==0))
   expect_equivalent(sav$yi, 0.4619, tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, 0.0570, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='D2ORL/D2ORN'", {

   dat <- dat.gibson2002

   sav <- escalc(measure="D2ORL", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(sav$yi, c(-0.4315, -0.9285,  0.5932, -0.1890), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c( 0.1276,  0.0493,  0.3204,  0.0690), tolerance=.tol[["var"]])

   sav <- escalc(measure="D2ORN", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   expect_equivalent(sav$yi, c(-0.3925, -0.8447,  0.5397, -0.1719), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c( 0.1056,  0.0408,  0.2651,  0.0571), tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='COR/UCOR/ZCOR'", {

   dat <- dat.mcdaniel1994

   sav <- escalc(measure="COR", ri=ri, ni=ni, data=dat, subset=c(1,13,33,102))
   expect_equivalent(sav$yi, c(0.0000, 0.6200, 0.9900, -0.1300), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c(0.0082, 0.0271, 0.0001,  0.0242), tolerance=.tol[["var"]])

   sav <- escalc(measure="UCOR", ri=ri, ni=ni, data=dat, subset=c(1,13,33,102))
   expect_equivalent(sav$yi, c(0.0000, 0.6363, 0.9925, -0.1317), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c(0.0082, 0.0253, 0.0000,  0.0241), tolerance=.tol[["var"]])

   sav <- escalc(measure="UCOR", ri=ri, ni=ni, data=dat, vtype="UB", subset=c(1,13,33,102))
   expect_equivalent(sav$vi, c(0.0084, 0.0283, 0.0000,  0.0261), tolerance=.tol[["var"]])

   sav <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat, subset=c(1,13,33,102))
   expect_equivalent(sav$yi, c(0.0000, 0.7250, 2.6467, -0.1307), tolerance=.tol[["est"]])
   expect_equivalent(sav$vi, c(0.0083, 0.0833, 0.3333,  0.0263), tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='PCOR/ZPCOR/SPCOR'", {

   ### data from Aloe and Thompson (2013)
   dat <- data.frame(ti = c(4.61, 6.19, 4.07, -0.77, 1.16),
                     ni = c(218, 232, 156, 382, 259),
                     mi = c(4, 7, 6, 19, 15),
                     r2i = c(.240, .455, .500, .327, .117))

   dat <- escalc(measure="PCOR", ti=ti, ni=ni, mi=mi, data=dat)
   expect_equivalent(dat$yi[1], 0.3012, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi[1], 0.0039, tolerance=.tol[["var"]])

   dat <- escalc(measure="ZPCOR", ti=ti, ni=ni, mi=mi, data=dat)
   expect_equivalent(dat$yi[1], 0.3108, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi[1], 0.0047, tolerance=.tol[["var"]])

   dat <- escalc(measure="SPCOR", ti=ti, ni=ni, mi=mi, r2i=r2i, data=dat)
   expect_equivalent(dat$yi[1], 0.2754, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi[1], 0.0033, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='MC/SMCRH'", {

   dat <- escalc(measure="MC", m1i=26, m2i=22, sd1i=sqrt(30), sd2i=sqrt(20), ni=60, ri=0.7)
   expect_equivalent(dat$yi, 4.0000, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.2618, tolerance=.tol[["var"]])

   dat <- escalc(measure="SMCRH", m1i=26, m2i=22, sd1i=sqrt(30), sd2i=sqrt(20), ni=60, ri=0.7)
   expect_equivalent(dat$yi, 0.7210, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0129, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='PAS'", {

   dat <- escalc(measure="PAS", xi=10, ni=20)
   expect_equivalent(dat$yi, 0.7854, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0125, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='IRS/IRFT'", {

   dat <- escalc(measure="IRS", xi=10, ti=20)
   expect_equivalent(dat$yi, 0.7071, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0125, tolerance=.tol[["var"]])

   dat <- escalc(measure="IRFT", xi=10, ti=20)
   expect_equivalent(dat$yi, 0.7244, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0125, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='ROMC'", {

   dat <- escalc(measure="ROMC", m1i=26, m2i=22, sd1i=sqrt(30), sd2i=sqrt(20), ni=60, ri=0.7)
   expect_equivalent(dat$yi, 0.1671, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0004, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='MPRD'", {

   dat <- escalc(measure="MPRD", ai=20, bi=10, ci=5, di=20)
   expect_equivalent(dat$yi, 0.0909, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0048, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='MPRR'", {

   dat <- escalc(measure="MPRR", ai=20, bi=10, ci=5, di=20)
   expect_equivalent(dat$yi, 0.1823, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0200, tolerance=.tol[["var"]])


})

test_that("escalc() works correctly for measure='MPOR'", {

   dat <- escalc(measure="MPOR", ai=20, bi=10, ci=5, di=20)
   expect_equivalent(dat$yi, 0.3646, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0782, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='MPORC'", {

   dat <- escalc(measure="MPORC", ai=20, bi=10, ci=5, di=20)
   expect_equivalent(dat$yi, 0.6931, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.3000, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='MPPETO'", {

   dat <- escalc(measure="MPPETO", ai=20, bi=10, ci=5, di=20)
   expect_equivalent(dat$yi, 0.6667, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.2667, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='IRSD'", {

   dat <- escalc(measure="IRSD", x1i=10, x2i=6, t1i=20, t2i=20)
   expect_equivalent(dat$yi, 0.1594, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0250, tolerance=.tol[["var"]])

})

test_that("escalc() works correctly for measure='MNLN/CVLN/SDLN'", {

   dat <- escalc(measure="MNLN", mi=10, sdi=2, ni=20)
   expect_equivalent(dat$yi, 2.3026, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0020, tolerance=.tol[["var"]])

   dat <- escalc(measure="CVLN", mi=10, sdi=2, ni=20)
   expect_equivalent(dat$yi, -1.5831, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0283, tolerance=.tol[["var"]])

   dat <- escalc(measure="SDLN", sdi=2, ni=20)
   expect_equivalent(dat$yi, 0.7195, tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, 0.0263, tolerance=.tol[["var"]])

})

test_that("'var.names' argument works correctly for 'escalc' objects.", {

   dat <- dat.bcg
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(author, ", ", year))
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y2","v2"), slab=paste0(author, ", ", year))
   expect_identical(tail(names(dat), 4), c("y1","v1","y2","v2"))
   expect_identical(attributes(dat)$yi.names, c("y2","y1"))
   expect_identical(attributes(dat)$vi.names, c("v2","v1"))
   expect_identical(attr(dat$y1, "measure"), "RR")
   expect_identical(attr(dat$y2, "measure"), "OR")

})

test_that("`[`, cbind(), and rbind() work correctly for 'escalc' objects.", {

   dat <- dat.bcg
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(author, ", ", year))
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y2","v2"), slab=paste0(author, ", ", year))
   dat <- cbind(dat[,1:9], dat[,c(12:13,10:11)])
   expect_identical(tail(names(dat), 4), c("y2","v2","y1","v1"))
   expect_identical(attributes(dat)$yi.names, c("y2","y1"))
   expect_identical(attributes(dat)$vi.names, c("v2","v1"))
   expect_identical(attr(dat$y1, "measure"), "RR")
   expect_identical(attr(dat$y2, "measure"), "OR")

   dat <- rbind(dat[13,], dat[1:12,])
   expect_equivalent(attr(dat$y2, "ni"), rowSums(dat[,c("tpos", "tneg", "cpos", "cneg")]))
   expect_identical(attr(dat$y2, "slab"), paste0(dat$author, ", ", dat$year))

   dat <- dat.bcg
   dat1 <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(author, ", ", year))
   dat2 <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(author, ", ", year))
   dat1 <- dat1[1:4,]
   dat2 <- dat2[4:1,]
   dat <- rbind(dat1, dat2)
   expect_equivalent(attr(dat$y1, "ni"), rowSums(dat[,c("tpos", "tneg", "cpos", "cneg")]))
   attr(dat1$y1, "ni") <- NULL
   dat <- rbind(dat1, dat2)
   expect_null(attr(dat$y1, "ni"))

})

test_that("summary() of 'escalc' objects works correctly with the 'out.names' argument.", {

   dat <- dat.bcg
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(author, ", ", year))
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y2","v2"), slab=paste0(author, ", ", year))
   dat <- summary(dat, var.names=c("y1","v1"), out.names=c("sei1","zi1","pval1","ci.lb1","ci.ub1"))
   dat <- summary(dat, var.names=c("y2","v2"), out.names=c("sei2","zi2","pval2","ci.lb2","ci.ub2"))
   expect_equivalent(with(dat, c(zi1[1], sei1[1], ci.lb1[1], ci.ub1[1])), c(-1.5586, 0.5706, -2.0077, 0.2290), tolerance=.tol[["est"]])
   expect_equivalent(with(dat, c(zi2[1], sei2[1], ci.lb2[1], ci.ub2[1])), c(-1.5708, 0.5976, -2.1100, 0.2326), tolerance=.tol[["est"]])

   dat <- dat[,1:11]
   expect_identical(attr(dat, "yi.names"), "y1")
   expect_identical(attr(dat, "vi.names"), "v1")

})

test_that("'subset' and 'include' arguments work correctly in 'escalc'.", {

   all <- dat.bcg
   all$tpos[1] <- NA

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=all, subset=1:4)
   expect_equivalent(c(dat$yi), c(NA, -1.5854, -1.3481, -1.4416), tolerance=.tol[["est"]])

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=all, subset=1:4, include=1:3)
   expect_equivalent(c(dat$yi), c(NA, -1.5854, -1.3481, NA), tolerance=.tol[["est"]])
   expect_identical(attributes(dat$yi)$ni, c(NA, 609L, 451L, NA))

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=all, subset=1:4, include=1:3, add.measure=TRUE)
   expect_identical(dat$measure, c("", "RR", "RR", ""))
   attributes(dat$yi)$ni[3] <- 1L
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, include=3:4, add.measure=TRUE)
   expect_equivalent(c(dat$yi), c(NA, -1.5854, -1.3863, -1.4564), tolerance=.tol[["est"]])
   expect_identical(dat$measure, c("", "RR", "OR", "OR"))
   expect_identical(attributes(dat$yi)$ni, c(NA, 609L, 451L, 26465L))

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=all, subset=1:4, include=1:3, add.measure=TRUE)
   attributes(dat$yi)$ni[3] <- 1L
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, include=3:4, replace=FALSE, add.measure=TRUE)
   expect_equivalent(c(dat$yi), c(NA, -1.5854, -1.3481, -1.4564), tolerance=.tol[["est"]])
   expect_identical(dat$measure, c("", "RR", "RR", "OR"))
   expect_identical(attributes(dat$yi)$ni, c(NA, 609L, 1L, 26465L))

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=all, subset=1:4, include=1:3, append=FALSE, add.measure=TRUE)
   expect_equivalent(c(dat$yi), c(NA, -1.5854, -1.3481, NA), tolerance=.tol[["est"]])
   expect_identical(dat$measure, c("", "RR", "RR", ""))

})

rm(list=ls())
