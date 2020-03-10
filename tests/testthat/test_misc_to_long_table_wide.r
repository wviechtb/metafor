### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: to.long() function")

source("tolerances.r") # read in tolerances

test_that("to.long() works correctly for measure='MD'", {

   dat <- dat.normand1999

   sav <- to.long(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   sav <- sav[,c(1,10:13)]

   expected <- structure(list(study = c(1, 1, 2, 2, 3, 3, 4, 4), group = structure(c(1, 2, 1, 2, 1, 2, 1, 2), .Label = c("1", "2"), class = "factor"), mean = c(55, 75, 27, 29, 64, 119, 66, 137), sd = c(47, 64, 7, 4, 17, 29, 20, 48), n = c(155, 156, 31, 32, 75, 71, 18, 18)), .Names = c("study", "group", "mean", "sd", "n"), class = "data.frame", row.names = c(NA, 8))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='MD'", {

   dat <- dat.normand1999

   sav <- to.table(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)

   expected <- structure(c(55, 75, 47, 64, 155, 156, 27, 29, 7, 4, 31, 32, 64, 119, 17, 29, 75, 71, 66, 137, 20, 48, 18, 18), .Dim = 2:4, .Dimnames = list(c("Grp1", "Grp2"), c("Mean", "SD", "n"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='COR'", {

   dat <- dat.molloy2014

   sav <- to.long(measure="COR", ri=ri, ni=ni, data=dat, subset=1:4)
   sav <- sav[,c(11:13)]

   expected <- structure(list(study = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), r = c(0.187, 0.162, 0.34, 0.32), n = c(109, 749, 55, 107)), .Names = c("study", "r", "n"), class = "data.frame", row.names = c(NA, 4))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='COR'", {

   dat <- dat.molloy2014

   sav <- to.table(measure="COR", ri=ri, ni=ni, data=dat, subset=1:4)

   expected <- structure(c(0.187, 109, 0.162, 749, 0.34, 55, 0.32, 107), .Dim = c(1L, 2L, 4L), .Dimnames = list("Grp", c("r", "n"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='PR'", {

   dat <- dat.debruin2009

   sav <- to.long(measure="PR", xi=xi, ni=ni, data=dat, subset=1:4)
   sav <- sav[,c(11:13)]

   expected <- structure(list(study = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), out1 = c(11L, 24L, 179L, 82L), out2 = c(18L, 9L, 147L, 158L)), class = "data.frame", row.names = c(NA, 4L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='PR'", {

   dat <- dat.debruin2009

   sav <- to.table(measure="PR", xi=xi, ni=ni, data=dat, subset=1:4)

   expected <- structure(c(11, 18, 24, 9, 179, 147, 82, 158), .Dim = c(1, 2, 4), .Dimnames = list("Grp", c("Out1", "Out2"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='IR'", {

   dat <- dat.hart1999

   sav <- to.long(measure="IR", xi=x1i, ti=t1i, data=dat, subset=1:4)
   sav <- sav[,c(1,14:15)]

   expected <- structure(list(trial = 1:4, events = c(9, 8, 3, 6), ptime = c(413, 263, 487, 237)), .Names = c("trial", "events", "ptime"), class = "data.frame", row.names = c(NA, 4))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='IR'", {

   dat <- dat.hart1999

   sav <- to.table(measure="IR", xi=x1i, ti=t1i, data=dat, subset=1:4)

   expected <- structure(c(9, 413, 8, 263, 3, 487, 6, 237), .Dim = c(1, 2, 4), .Dimnames = list("Grp", c("Events", "Person-Time"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='MN'", {

   dat <- dat.normand1999

   sav <- to.long(measure="MN", mi=m1i, sdi=sd1i, ni=n1i, data=dat, subset=1:4)
   sav <- sav[,c(1,10:12)]

   expected <- structure(list(study = 1:4, mean = c(55, 27, 64, 66), sd = c(47, 7, 17, 20), n = c(155, 31, 75, 18)), .Names = c("study", "mean", "sd", "n"), class = "data.frame", row.names = c(NA, 4))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='MN'", {

   dat <- dat.normand1999

   sav <- to.table(measure="MN", mi=m1i, sdi=sd1i, ni=n1i, data=dat, subset=1:4)

   expected <- structure(c(55, 47, 155, 27, 7, 31, 64, 17, 75, 66, 20, 18), .Dim = c(1, 3, 4), .Dimnames = list("Grp", c("Mean", "SD", "n"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

### create dataset (from Morris, 2008)

datT <- data.frame(
m_pre   = c(30.6, 23.5, 0.5, 53.4, 35.6),
m_post  = c(38.5, 26.8, 0.7, 75.9, 36.0),
sd_pre  = c(15.0, 3.1, 0.1, 14.5, 4.7),
sd_post = c(11.6, 4.1, 0.1, 4.4, 4.6),
ni      = c(20, 50, 9, 10, 14),
ri      = c(.47, .64, .77, .89, .44))

test_that("to.long() works correctly for measure='SMCR'", {

   sav <- to.long(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ni=ni, ri=ri, data=datT, subset=2:4)
   sav <- sav[,c(7:12)]

   expected <- structure(list(study = structure(1:3, .Label = c("2", "3", "4"), class = "factor"), mean1 = c(26.8, 0.7, 75.9), mean2 = c(23.5, 0.5, 53.4), sd1 = c(3.1, 0.1, 14.5), n = c(50, 9, 10), r = c(0.64, 0.77, 0.89)), .Names = c("study", "mean1", "mean2", "sd1", "n", "r"), class = "data.frame", row.names = c(NA, 3))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='SMCR'", {

   sav <- to.table(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ni=ni, ri=ri, data=datT, subset=2:4)

   expected <- structure(c(26.8, 23.5, 3.1, 50, 0.64, 0.7, 0.5, 0.1, 9, 0.77, 75.9, 53.4, 14.5, 10, 0.89), .Dim = c(1, 5, 3), .Dimnames = list("Grp", c("Mean1", "Mean2", "SD1", "n", "r"), c("2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='ARAW'", {

   dat <- dat.bonett2010

   sav <- to.long(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat, subset=1:4)
   sav <- sav[,c(1,8:10)]

   expected <- structure(list(study = 1:4, alpha = c(0.93, 0.91, 0.94, 0.89), m = c(20, 20, 20, 20), n = c(103, 64, 118, 401)), .Names = c("study", "alpha", "m", "n"), class = "data.frame", row.names = c(NA, 4))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='ARAW'", {

   dat <- dat.bonett2010

   sav <- to.table(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat, subset=1:4)

   expected <- structure(c(0.93, 20, 103, 0.91, 20, 64, 0.94, 20, 118, 0.89, 20, 401), .Dim = c(1, 3, 4), .Dimnames = list("Grp", c("alpha", "m", "n"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.wide() works correctly.", {

   dat.l <- dat.hasselblad1998
   dat.c <- to.wide(dat.l, study="study", grp="trt", ref="no_contact", grpvars=6:7)

   expect_equivalent(dat.c$xi.1, c(363, 10, 23, 9, 237, 9, 16, 31, 26, 29, 12, 17, 77, 21, 107, 20, 3, 32, 8, 34, 9, 19, 143, 36, 73, 54))
   expect_equivalent(dat.c$xi.2, c(75, 9, 9, 2, 58, 0, 20, 3, 1, 11, 11, 6, 79, 18, 64, 12, 9, 7, 5, 20, 0, 8, 95, 15, 78, 69))
   expect_equivalent(dat.c$comp, c("in-no", "gr-no", "in-no", "in-no", "in-no", "in-no", "in-se", "in-no", "in-no", "gr-se", "in-se", "in-no", "se-no", "se-no", "in-no", "gr-in", "gr-in", "gr-se", "in-no", "in-no", "gr-no", "se-no", "in-no", "in-no", "in-no", "in-no"))
   expect_equivalent(dat.c$design, c("in-no", "gr-in-no", "gr-in-no", "in-no", "in-no", "in-no", "in-se", "in-no", "in-no", "gr-in-se", "gr-in-se", "in-no", "se-no", "se-no", "in-no", "gr-in", "gr-in", "gr-se", "in-no", "in-no", "gr-no", "se-no", "in-no", "in-no", "in-no", "in-no"))

   dat.l$trt <- factor(dat.l$trt, levels=c("no_contact", "ind_counseling", "grp_counseling", "self_help"))
   dat.c <- to.wide(dat.l, study="study", grp="trt", grpvars=5:7, postfix=c(".T",".C"), minlen=1)

   expect_equivalent(dat.c$xi.T, c(363, 23, 10, 9, 237, 9, 16, 31, 26, 12, 29, 17, 77, 21, 107, 12, 9, 32, 8, 34, 9, 19, 143, 36, 73, 54))
   expect_equivalent(dat.c$xi.C, c(75, 9, 9, 2, 58, 0, 20, 3, 1, 11, 11, 6, 79, 18, 64, 20, 3, 7, 5, 20, 0, 8, 95, 15, 78, 69))
   expect_equivalent(dat.c$comp, c("i-n", "i-n", "g-n", "i-n", "i-n", "i-n", "i-s", "i-n", "i-n", "i-s", "g-s", "i-n", "s-n", "s-n", "i-n", "i-g", "i-g", "g-s", "i-n", "i-n", "g-n", "s-n", "i-n", "i-n", "i-n", "i-n"))
   expect_equivalent(dat.c$design, c("i-n", "i-g-n", "i-g-n", "i-n", "i-n", "i-n", "i-s", "i-n", "i-n", "i-g-s", "i-g-s", "i-n", "s-n", "s-n", "i-n", "i-g", "i-g", "g-s", "i-n", "i-n", "g-n", "s-n", "i-n", "i-n", "i-n", "i-n"))

})
