### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: escalc() function")

test_that("escalc() works correctly for measure='RR'", {

   data(dat.bcg, package="metafor")

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   expect_equivalent(round(dat$yi[1],4), -0.8893)
   expect_equivalent(round(dat$vi[1],4),  0.3256)

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

test_that("escalc() works correctly for measure='PCOR/ZPCOR/SPCOR'", {

   ### data from Aloe and Thompson (2013)
   dat <- data.frame(ti = c(4.61, 6.19, 4.07, -0.77, 1.16),
                     ni = c(218, 232, 156, 382, 259),
                     mi = c(4, 7, 6, 19, 15),
                     r2i = c(.240, .455, .500, .327, .117))

   dat <- escalc(measure="PCOR", ti=ti, ni=ni, mi=mi, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.3012, 0.0039))

   dat <- escalc(measure="ZPCOR", ti=ti, ni=ni, mi=mi, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.3108, 0.0047))

   dat <- escalc(measure="SPCOR", ti=ti, ni=ni, mi=mi, r2i=r2i, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.2754, 0.0033))

})

test_that("escalc() works correctly for measure='SMCRH'", {

   dat <- escalc(measure="SMCRH", m1i=26, m2i=22, sd1i=sqrt(30), sd2i=sqrt(20), ni=60, ri=0.7)
   expect_equivalent(round(c(dat$yi, dat$vi), 4), c(0.7210, 0.0129))

})

test_that("escalc() works correctly for measure='ROMC'", {

   dat <- escalc(measure="ROMC", m1i=26, m2i=22, sd1i=sqrt(30), sd2i=sqrt(20), ni=60, ri=0.7)
   expect_equivalent(round(c(dat$yi, dat$vi), 4), c(0.1671, 0.0004))

})

test_that("escalc() with formula works correctly for measure='IRR'", {

   dat <- get(data(dat.nielweise2008, package="metafor"))
   dat <- to.long(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat)
   dat <- escalc(measure="IRR", events/ptime ~ group | factor(study), data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(-0.0605, 0.2338))

})

test_that("escalc() with formula works correctly for measure='SMD'", {

   dat <- get(data(dat.normand1999, package="metafor"))
   dat <- to.long(measure="SMD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat)
   dat <- escalc(measure="SMD", mean/sd ~ group | factor(study), weights=n, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(-0.3552, 0.0131))

})

test_that("escalc() with formula works correctly for measure='COR'", {

   dat <- get(data(dat.mcdaniel1994, package="metafor"))
   dat <- to.long(measure="COR", ri=ri, ni=ni, data=dat)
   dat <- escalc(measure="COR", ri ~ 1 | factor(study), weights=n, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.0000, 0.0082))

})

test_that("escalc() with formula works correctly for measure='PR'", {

   dat <- get(data(dat.debruin2009, package="metafor"))
   dat <- to.long(measure="PR", xi=xi, ni=ni, data=dat, vlong=TRUE)
   dat <- escalc(measure="PR", outcome ~ 1 | factor(study), weights=freq, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.3793, 0.0081))

})

test_that("escalc() with formula works correctly for measure='IR'", {

   dat <- get(data(dat.hart1999, package="metafor"))
   dat <- to.long(measure="IR", xi=x1i, ti=t1i, data=dat)
   dat <- escalc(measure="IR", events/ptime ~ 1 | factor(trial), data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.0218, 0.0001))

})

test_that("escalc() with formula works correctly for measure='MN'", {

   dat <- get(data(dat.normand1999, package="metafor"))
   dat <- to.long(measure="MN", mi=m1i, sdi=sd1i, ni=n1i, data=dat)
   dat <- escalc(measure="MN", mean/sd ~ 1 | factor(study), weights=n, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(55.0000, 14.2516))

})

test_that("escalc() with formula works correctly for measure='ARAW/AHW/ABT'", {

   dat <- get(data(dat.bonett2010, package="metafor"))
   dat <- to.long(measure="ARAW", ai=ai, mi=mi, ni=ni, data=dat)
   dat <- escalc(measure="ARAW", alpha/m ~ 1 | factor(study), weights=n, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.9300, 0.0001))

   dat <- get(data(dat.bonett2010, package="metafor"))
   dat <- escalc(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(0.5879, 0.0004))

   dat <- escalc(measure="ABT", ai=ai, mi=mi, ni=ni, data=dat)
   expect_equivalent(round(c(dat$yi[1], dat$vi[1]), 4), c(2.6593, 0.0208))

})

test_that("'var.names' argument works correctly for 'escalc' objects.", {

   dat <- get(data(dat.bcg, package="metafor"))
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(dat$author, ", ", dat$year))
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y2","v2"), slab=paste0(dat$author, ", ", dat$year))
   expect_identical(tail(names(dat), 4), c("y1","v1","y2","v2"))
   expect_identical(attributes(dat)$yi.names, c("y2","y1"))
   expect_identical(attributes(dat)$vi.names, c("v2","v1"))
   expect_identical(attr(dat$y1, "measure"), "RR")
   expect_identical(attr(dat$y2, "measure"), "OR")

})

test_that("`[`, cbind(), and rbind() work correctly for 'escalc' objects.", {

   dat <- get(data(dat.bcg, package="metafor"))
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(dat$author, ", ", dat$year))
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y2","v2"), slab=paste0(dat$author, ", ", dat$year))
   dat <- cbind(dat[,1:9], dat[,c(12:13,10:11)])
   expect_identical(tail(names(dat), 4), c("y2","v2","y1","v1"))
   expect_identical(attributes(dat)$yi.names, c("y2","y1"))
   expect_identical(attributes(dat)$vi.names, c("v2","v1"))
   expect_identical(attr(dat$y1, "measure"), "RR")
   expect_identical(attr(dat$y2, "measure"), "OR")

   dat <- rbind(dat[13,], dat[1:12,])
   expect_equivalent(attr(dat$y2, "ni"), rowSums(dat[,c("tpos", "tneg", "cpos", "cneg")]))
   expect_identical(attr(dat$y2, "slab"), paste0(dat$author, ", ", dat$year))

   dat <- get(data(dat.bcg, package="metafor"))
   dat1 <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(dat$author, ", ", dat$year))
   dat2 <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(dat$author, ", ", dat$year))
   dat1 <- dat1[1:4,]
   dat2 <- dat2[4:1,]
   dat <- rbind(dat1, dat2)
   expect_equivalent(attr(dat$y1, "ni"), rowSums(dat[,c("tpos", "tneg", "cpos", "cneg")]))
   attr(dat1$y1, "ni") <- NULL
   dat <- rbind(dat1, dat2)
   expect_null(attr(dat$y1, "ni"))

})

test_that("summary() of 'escalc' objects works correctly with the 'out.names' argument.", {

   dat <- get(data(dat.bcg, package="metafor"))
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y1","v1"), slab=paste0(dat$author, ", ", dat$year))
   dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, var.names=c("y2","v2"), slab=paste0(dat$author, ", ", dat$year))
   dat <- summary(dat, var.names=c("y1","v1"), out.names=c("sei1","zi1","ci.lb1","ci.ub1"))
   dat <- summary(dat, var.names=c("y2","v2"), out.names=c("sei2","zi2","ci.lb2","ci.ub2"))
   expect_equivalent(round(with(dat, c(zi1[1], sei1[1], ci.lb1[1], ci.ub1[1])), 4), c(-1.5586, 0.5706, -2.0077, 0.2290))
   expect_equivalent(round(with(dat, c(zi2[1], sei2[1], ci.lb2[1], ci.ub2[1])), 4), c(-1.5708, 0.5976, -2.1100, 0.2326))

   dat <- dat[,1:11]
   expect_identical(attr(dat, "yi.names"), "y1")
   expect_identical(attr(dat, "vi.names"), "v1")

})
