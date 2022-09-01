### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: emmprep() function")

source("settings.r")

test_that("emmprep() gives correct results for an intercept-only model.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, data=dat)

   sav <- capture.output(emmprep(res, verbose=TRUE))
   sav <- emmprep(res)

   skip_on_cran()

   tmp <- emmeans::emmeans(sav, specs="1", type="response")
   tmp <- as.data.frame(tmp)

   expect_equivalent(tmp$response,  0.4894209, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$asymp.LCL, 0.3440743, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$asymp.UCL, 0.6961661, tolerance=.tol[["ci"]])

})

test_that("emmprep() gives correct results for a meta-regression model.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   dat$yi[1] <- NA
   res <- suppressWarnings(rma(yi, vi, mods = ~ ablat + alloc, data=dat, subset=-2, test="knha"))

   sav <- emmprep(res)

   skip_on_cran()

   tmp <- emmeans::emmeans(sav, specs="1", type="response")
   tmp <- as.data.frame(tmp)

   expect_equivalent(tmp$response,  0.5395324, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$lower.CL,  0.3564229, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$upper.CL,  0.8167130, tolerance=.tol[["ci"]])

   sav <- emmprep(res, data=dat[-c(1,2),], df=7, sigma=sqrt(res$tau2), tran="log")
   tmp <- as.data.frame(tmp)

   expect_equivalent(tmp$response,  0.5395324, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$lower.CL,  0.3564229, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$upper.CL,  0.8167130, tolerance=.tol[["ci"]])

})

test_that("emmprep() gives correct results for the r-to-z transformation.", {

   dat <- dat.mcdaniel1994
   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat)

   res <- suppressWarnings(rma(yi, vi, mods = ~ factor(type), data=dat, test="knha"))

   sav <- emmprep(res)

   skip_on_cran()

   tmp <- emmeans::emmeans(sav, specs="1", type="response")
   tmp <- as.data.frame(tmp)

   expect_equivalent(tmp$response,  0.2218468, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$lower.CL,  0.1680606, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$upper.CL,  0.2743160, tolerance=.tol[["ci"]])

})

rm(list=ls())
