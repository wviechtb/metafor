### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: anova() function")

test_that("anova() works correctly for comparing nested models.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res1 <- rma(yi, vi, data=dat, method="ML")
   res2 <- rma(yi ~ ablat, vi, data=dat, method="ML")
   sav <- anova(res1, res2)
   out <- capture.output(print(sav))

   expect_equivalent(round(c(sav$LRT), 4), 9.9588)

   res1 <- rma(yi, vi, data=dat, method="REML")
   res2 <- rma(yi ~ ablat, vi, data=dat, method="REML")
   expect_warning(sav <- anova(res1, res2))

   expect_equivalent(round(c(sav$LRT), 4), 8.2301)

})

test_that("anova() works correctly when using the 'btt' argument.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)
   sav <- anova(res, btt=3:4)
   out <- capture.output(print(sav))

   expect_equivalent(round(c(sav$QM), 4), 1.2850)
   expect_equivalent(round(c(sav$QMp), 2), 0.53)

   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat, test="knha")
   sav <- anova(res, btt=3:4)
   out <- capture.output(print(sav))

   expect_equivalent(round(c(sav$QM), 4), 0.6007)
   expect_equivalent(round(c(sav$QMp), 2), 0.57)

})

test_that("anova() works correctly when using the 'L' argument.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat)
   sav <- anova(res, L=rbind(c(1, 10, 0, 0), c(1, 30, 0, 0), c(1, 50, 0, 0)))
   out <- capture.output(print(sav))

   expect_equivalent(round(c(sav$zval), 4), c(0.0588, -1.7964, -3.1210))

   res <- rma(yi, vi, mods = ~ ablat + alloc, data=dat, test="knha")
   sav <- anova(res, L=rbind(c(1, 10, 0, 0), c(1, 10, 1, 0), c(1, 10, 0, 1)))
   out <- capture.output(print(sav))

   expect_equivalent(round(c(sav$zval), 4), c(0.0568, -0.8252, 0.2517))
   expect_equivalent(round(c(sav$QM), 4), 0.4230)
   expect_equivalent(round(c(sav$QMp), 2), 0.74)

})
