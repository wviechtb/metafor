### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: handling of NAs")

source("settings.r")

dat <- data.frame(yi = c(NA, 1, 3, 2, 5, 4, 6), vi = c(1, NA, 1, 1, 1, 1, 1), xi = c(0, 1, NA, 3, 4, 5, 6))

test_that("NAs are correctly handled by various method functions for rma.uni() intercept-only models.", {

   expect_warning(res <- rma(yi, vi, data=dat))

   expect_equivalent(res$k, 5)

   options(na.action = "na.omit")

   expect_equivalent(fitted(res), c(4, 4, 4, 4, 4))
   expect_equivalent(resid(res), c(-1, -2, 1, 0, 2))
   expect_equivalent(predict(res)$pred, 4)
   expect_equivalent(blup(res)$pred, c(3.4, 2.8, 4.6, 4.0, 5.2))
   expect_equivalent(cooks.distance(res), c(0.125, 0.5, 0.125, 0, 0.5))
   expect_equivalent(dfbetas(res)[[1]], c(-0.3273, -0.8660, 0.3273, 0, 0.8660), tolerance=.tol[["inf"]])
   expect_equivalent(hatvalues(res), c(0.2, 0.2, 0.2, 0.2, 0.2))
   expect_equivalent(leave1out(res)$estimate, c(4.25, 4.5, 3.75, 4, 3.5))
   expect_equivalent(ranef(res)$pred, c(-0.6, -1.2, 0.6, 0, 1.2))
   expect_equivalent(rstandard(res)$resid, c(-1, -2, 1, 0, 2))
   expect_equivalent(rstudent(res)$resid, c(-1.25, -2.5, 1.25, 0, 2.5))
   expect_equivalent(length(simulate(res, seed=1234)[[1]]), 5)
   expect_equivalent(weights(res), c(20, 20, 20, 20, 20))

   options(na.action = "na.pass")

   # note: all of these are of the same length as the original data (except for predict(), which gives a single value for intercept-only models)

   expect_equivalent(fitted(res), c(4, 4, 4, 4, 4, 4, 4)) # note: can compute fitted value even for the study with missing yi and the study with missing vi
   expect_equivalent(resid(res), c(NA, -3, -1, -2, 1, 0, 2)) # note: can compute residual value even for the study with missing vi
   expect_equivalent(predict(res)$pred, 4)
   expect_equivalent(blup(res)$pred, c(NA, NA, 3.4, 2.8, 4.6, 4.0, 5.2))
   expect_equivalent(cooks.distance(res), c(NA, NA, 0.125, 0.5, 0.125, 0, 0.5))
   expect_equivalent(dfbetas(res)[[1]], c(NA, NA, -0.3273, -0.8660, 0.3273, 0, 0.8660), tolerance=.tol[["inf"]])
   expect_equivalent(hatvalues(res), c(NA, NA, 0.2, 0.2, 0.2, 0.2, 0.2))
   expect_equivalent(leave1out(res)$estimate, c(NA, NA, 4.25, 4.5, 3.75, 4, 3.5))
   expect_equivalent(ranef(res)$pred, c(NA, NA, -0.6, -1.2, 0.6, 0, 1.2))
   expect_equivalent(rstandard(res)$resid, c(NA, NA, -1, -2, 1, 0, 2))
   expect_equivalent(rstudent(res)$resid, c(NA, NA, -1.25, -2.5, 1.25, 0, 2.5))
   expect_equivalent(length(simulate(res, seed=1234)[[1]]), 7)
   expect_equivalent(weights(res), c(NA, NA, 20, 20, 20, 20, 20))

   options(na.action = "na.exclude")

   # note: all of these are of the same length as the original data, but are NA for studies 1 and 2

   expect_equivalent(fitted(res), c(NA, NA, 4, 4, 4, 4, 4)) # note: all of these are of the same length as the original data, but are NA for studies 1 and 2
   expect_equivalent(resid(res), c(NA, NA, -1, -2, 1, 0, 2))
   expect_equivalent(predict(res)$pred, 4)
   expect_equivalent(blup(res)$pred, c(NA, NA, 3.4, 2.8, 4.6, 4.0, 5.2))
   expect_equivalent(cooks.distance(res), c(NA, NA, 0.125, 0.5, 0.125, 0, 0.5))
   expect_equivalent(dfbetas(res)[[1]], c(NA, NA, -0.3273, -0.8660, 0.3273, 0, 0.8660), tolerance=.tol[["inf"]])
   expect_equivalent(hatvalues(res), c(NA, NA, 0.2, 0.2, 0.2, 0.2, 0.2))
   expect_equivalent(leave1out(res)$estimate, c(NA, NA, 4.25, 4.5, 3.75, 4, 3.5))
   expect_equivalent(ranef(res)$pred, c(NA, NA, -0.6, -1.2, 0.6, 0, 1.2))
   expect_equivalent(rstandard(res)$resid, c(NA, NA, -1, -2, 1, 0, 2))
   expect_equivalent(rstudent(res)$resid, c(NA, NA, -1.25, -2.5, 1.25, 0, 2.5))
   expect_equivalent(length(simulate(res, seed=1234)[[1]]), 7)
   expect_equivalent(weights(res), c(NA, NA, 20, 20, 20, 20, 20))

   options(na.action = "na.omit")

})

test_that("NAs are correctly handled by various method functions for rma.uni() meta-regression models.", {

   expect_warning(res <- rma(yi, vi, mods = ~ xi, data=dat))

   expect_equivalent(res$k, 4)

   options(na.action = "na.omit")

   expect_equivalent(fitted(res), c(2.6, 3.7, 4.8, 5.9))
   expect_equivalent(resid(res), c(-0.6, 1.3, -0.8, 0.1))
   expect_equivalent(predict(res)$pred, c(2.6, 3.7, 4.8, 5.9))
   expect_equivalent(blup(res)$pred, c(2.4444, 4.0370, 4.5926, 5.9259), tolerance=.tol[["pred"]])
   expect_equivalent(cooks.distance(res), c(2.0741, 0.7664, 0.2902, 0.0576), tolerance=.tol[["inf"]])
   expect_equivalent(dfbetas(res)[[2]], c(1.0954, -0.4153, -0.1912, 0.1369), tolerance=.tol[["inf"]])
   expect_equivalent(hatvalues(res), c(0.7, 0.3, 0.3, 0.7))
   expect_equivalent(ranef(res)$pred, c(-0.1556, 0.3370, -0.2074, 0.0259), tolerance=.tol[["pred"]])
   expect_equivalent(rstandard(res)$resid, c(-0.6, 1.3, -0.8, 0.1))
   expect_equivalent(rstudent(res)$resid, c(-2, 1.8571, -1.1429, 0.3333), tolerance=.tol[["pred"]])
   expect_equivalent(length(simulate(res, seed=1234)[[1]]), 4)
   expect_equivalent(weights(res), c(25, 25, 25, 25))

   options(na.action = "na.pass")

   # note: all of these are of the same length as the original data

   expect_equivalent(fitted(res), c(-0.7, 0.4, NA, 2.6, 3.7, 4.8, 5.9)) # note: can compute fitted value even for the study with missing yi and the study with missing vi
   expect_equivalent(resid(res), c(NA, 0.6, NA, -0.6, 1.3, -0.8, 0.1)) # note: can compute residual value even for the study with missing vi
   expect_equivalent(predict(res)$pred, c(-0.7, 0.4, NA, 2.6, 3.7, 4.8, 5.9))
   expect_equivalent(blup(res)$pred, c(NA, NA, NA, 2.4444, 4.0370, 4.5926, 5.9259), tolerance=.tol[["pred"]])
   expect_equivalent(cooks.distance(res), c(NA, NA, NA, 2.0741, 0.7664, 0.2902, 0.0576), tolerance=.tol[["inf"]])
   expect_equivalent(dfbetas(res)[[2]], c(NA, NA, NA, 1.0954, -0.4153, -0.1912, 0.1369), tolerance=.tol[["inf"]])
   expect_equivalent(hatvalues(res), c(NA, NA, NA, 0.7, 0.3, 0.3, 0.7))
   expect_equivalent(ranef(res)$pred, c(NA, NA, NA, -0.1556, 0.3370, -0.2074, 0.0259), tolerance=.tol[["pred"]])
   expect_equivalent(rstandard(res)$resid, c(NA, NA, NA, -0.6, 1.3, -0.8, 0.1))
   expect_equivalent(rstudent(res)$resid, c(NA, NA, NA, -2, 1.8571, -1.1429, 0.3333), tolerance=.tol[["pred"]])
   expect_equivalent(length(simulate(res, seed=1234)[[1]]), 7)
   expect_equivalent(weights(res), c(NA, NA, NA, 25, 25, 25, 25))

   options(na.action = "na.exclude")

   # note: all of these are of the same length as the original data, but are NA for studies 1, 2, and 3

   expect_equivalent(fitted(res), c(NA, NA, NA, 2.6, 3.7, 4.8, 5.9))
   expect_equivalent(resid(res), c(NA, NA, NA, -0.6, 1.3, -0.8, 0.1))
   expect_equivalent(predict(res)$pred, c(NA, NA, NA, 2.6, 3.7, 4.8, 5.9))
   expect_equivalent(blup(res)$pred, c(NA, NA, NA, 2.4444, 4.0370, 4.5926, 5.9259), tolerance=.tol[["pred"]])
   expect_equivalent(cooks.distance(res), c(NA, NA, NA, 2.0741, 0.7664, 0.2902, 0.0576), tolerance=.tol[["inf"]])
   expect_equivalent(dfbetas(res)[[2]], c(NA, NA, NA, 1.0954, -0.4153, -0.1912, 0.1369), tolerance=.tol[["inf"]])
   expect_equivalent(hatvalues(res), c(NA, NA, NA, 0.7, 0.3, 0.3, 0.7))
   expect_equivalent(ranef(res)$pred, c(NA, NA, NA, -0.1556, 0.3370, -0.2074, 0.0259), tolerance=.tol[["pred"]])
   expect_equivalent(rstandard(res)$resid, c(NA, NA, NA, -0.6, 1.3, -0.8, 0.1))
   expect_equivalent(rstudent(res)$resid, c(NA, NA, NA, -2, 1.8571, -1.1429, 0.3333), tolerance=.tol[["pred"]])
   expect_equivalent(length(simulate(res, seed=1234)[[1]]), 7)
   expect_equivalent(weights(res), c(NA, NA, NA, 25, 25, 25, 25))

   options(na.action = "na.omit")

})

test_that("NAs are correctly handled by rma.mv() intercept-only models.", {

   dat <- dat.konstantopoulos2011

   res1 <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=sparse)
   res2 <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat, sparse=sparse)
   expect_equivalent(logLik(res1), logLik(res2), tolerance=.tol[["fit"]])

   dat$yi[1:2] <- NA

   expect_warning(res1 <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=sparse))
   expect_warning(res2 <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat, sparse=sparse))
   expect_equivalent(logLik(res1), logLik(res2), tolerance=.tol[["fit"]])

   dat$yi[1:4] <- NA # entire district 11 is missing

   expect_warning(res1 <- rma.mv(yi, vi, random = ~ 1 | district/study, data=dat, sparse=sparse))
   expect_warning(res2 <- rma.mv(yi, vi, random = ~ factor(study) | district, data=dat, sparse=sparse))
   expect_equivalent(logLik(res1), logLik(res2), tolerance=.tol[["fit"]])

})

rm(list=ls())
