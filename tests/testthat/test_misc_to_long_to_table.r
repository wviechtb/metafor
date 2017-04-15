### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: to.long() function")

test_that("to.long() works correctly for measure='MD'", {

   dat <- get(data(dat.normand1999, package="metafor"))

   sav <- to.long(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)
   sav <- sav[,c(1,10:13)]

   expected <- structure(list(study = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L), group = structure(c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), .Label = c("1", "2"), class = "factor"), mean = c(55, 75, 27, 29, 64, 119, 66, 137), sd = c(47, 64, 7, 4, 17, 29, 20, 48), n = c(155, 156, 31, 32, 75, 71, 18, 18)), .Names = c("study", "group", "mean", "sd", "n"), class = "data.frame", row.names = c(NA, 8L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='MD'", {

   dat <- get(data(dat.normand1999, package="metafor"))

   sav <- to.table(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, subset=1:4)

   expected <- structure(c(55L, 75L, 47L, 64L, 155L, 156L, 27L, 29L, 7L, 4L, 31L, 32L, 64L, 119L, 17L, 29L, 75L, 71L, 66L, 137L, 20L, 48L, 18L, 18L), .Dim = 2:4, .Dimnames = list(c("Grp1", "Grp2"), c("Mean", "SD", "n"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='COR'", {

   dat <- get(data(dat.molloy2014, package="metafor"))

   sav <- to.long(measure="COR", ri=ri, ni=ni, data=dat, subset=1:4)
   sav <- sav[,c(11:13)]

   expected <- structure(list(study = structure(1:4, .Label = c("1", "2", "3", "4"), class = "factor"), r = c(0.187, 0.162, 0.34, 0.32), n = c(109, 749, 55, 107)), .Names = c("study", "r", "n"), class = "data.frame", row.names = c(NA, 4L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='COR'", {

   dat <- get(data(dat.molloy2014, package="metafor"))

   sav <- to.table(measure="COR", ri=ri, ni=ni, data=dat, subset=1:4)

   expected <- structure(c(0.187, 109, 0.162, 749, 0.34, 55, 0.32, 107), .Dim = c(1L, 2L, 4L), .Dimnames = list("Grp", c("r", "n"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='PR'", {

   dat <- get(data(dat.debruin2009, package="metafor"))

   sav <- to.long(measure="PR", xi=xi, ni=ni, data=dat, subset=1:4)
   sav <- sav[,c(11:13)]

   expected <- structure(list(study = 1:4, out1 = c(11L, 24L, 179L, 82L), out2 = c(18L, 9L, 147L, 158L)), .Names = c("study", "out1", "out2"), class = "data.frame", row.names = c(NA, 4L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='PR'", {

   dat <- get(data(dat.debruin2009, package="metafor"))

   sav <- to.table(measure="PR", xi=xi, ni=ni, data=dat, subset=1:4)

   expected <- structure(c(11L, 18L, 24L, 9L, 179L, 147L, 82L, 158L), .Dim = c(1L, 2L, 4L), .Dimnames = list("Grp", c("Out1", "Out2"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='IR'", {

   dat <- get(data(dat.hart1999, package="metafor"))

   sav <- to.long(measure="IR", xi=x1i, ti=t1i, data=dat, subset=1:4)
   sav <- sav[,c(1,14:15)]

   expected <- structure(list(trial = 1:4, events = c(9L, 8L, 3L, 6L), ptime = c(413L, 263L, 487L, 237L)), .Names = c("trial", "events", "ptime"), class = "data.frame", row.names = c(NA, 4L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='IR'", {

   dat <- get(data(dat.hart1999, package="metafor"))

   sav <- to.table(measure="IR", xi=x1i, ti=t1i, data=dat, subset=1:4)

   expected <- structure(c(9L, 413L, 8L, 263L, 3L, 487L, 6L, 237L), .Dim = c(1L, 2L, 4L), .Dimnames = list("Grp", c("Events", "Person-Time"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='MN'", {

   dat <- get(data(dat.normand1999, package="metafor"))

   sav <- to.long(measure="MN", mi=m1i, sdi=sd1i, ni=n1i, data=dat, subset=1:4)
   sav <- sav[,c(1,10:12)]

   expected <- structure(list(study = 1:4, mean = c(55L, 27L, 64L, 66L), sd = c(47L, 7L, 17L, 20L), n = c(155L, 31L, 75L, 18L)), .Names = c("study", "mean", "sd", "n"), class = "data.frame", row.names = c(NA, 4L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='MN'", {

   dat <- get(data(dat.normand1999, package="metafor"))

   sav <- to.table(measure="MN", mi=m1i, sdi=sd1i, ni=n1i, data=dat, subset=1:4)

   expected <- structure(c(55L, 47L, 155L, 27L, 7L, 31L, 64L, 17L, 75L, 66L, 20L, 18L), .Dim = c(1L, 3L, 4L), .Dimnames = list("Grp", c("Mean", "SD", "n"), c("1", "2", "3", "4")))
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

   expected <- structure(list(study = structure(1:3, .Label = c("2", "3", "4"), class = "factor"), mean1 = c(26.8, 0.7, 75.9), mean2 = c(23.5, 0.5, 53.4), sd1 = c(3.1, 0.1, 14.5), n = c(50, 9, 10), r = c(0.64, 0.77, 0.89)), .Names = c("study", "mean1", "mean2", "sd1", "n", "r"), class = "data.frame", row.names = c(NA, 3L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='SMCR'", {

   sav <- to.table(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ni=ni, ri=ri, data=datT, subset=2:4)

   expected <- structure(c(26.8, 23.5, 3.1, 50, 0.64, 0.7, 0.5, 0.1, 9, 0.77, 75.9, 53.4, 14.5, 10, 0.89), .Dim = c(1L, 5L, 3L), .Dimnames = list("Grp", c("Mean1", "Mean2", "SD1", "n", "r"), c("2", "3", "4")))
   expect_equivalent(sav, expected)

})

test_that("to.long() works correctly for measure='ARAW'", {

   dat <- get(data(dat.bonett2010, package="metafor"))

   sav <- to.long(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat, subset=1:4)
   sav <- sav[,c(1,8:10)]

   expected <- structure(list(study = 1:4, alpha = c(0.93, 0.91, 0.94, 0.89), m = c(20, 20, 20, 20), n = c(103, 64, 118, 401)), .Names = c("study", "alpha", "m", "n"), class = "data.frame", row.names = c(NA, 4L))
   expect_equivalent(sav, expected)

})

test_that("to.table() works correctly for measure='ARAW'", {

   dat <- get(data(dat.bonett2010, package="metafor"))

   sav <- to.table(measure="AHW", ai=ai, mi=mi, ni=ni, data=dat, subset=1:4)

   expected <- structure(c(0.93, 20, 103, 0.91, 20, 64, 0.94, 20, 118, 0.89, 20, 401), .Dim = c(1L, 3L, 4L), .Dimnames = list("Grp", c("alpha", "m", "n"), c("1", "2", "3", "4")))
   expect_equivalent(sav, expected)

})
