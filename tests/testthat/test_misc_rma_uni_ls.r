### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma() function with location-scale models")

test_that("location-scale model results are correct for in intercept-only model", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi, vi, data=dat, test="t")
   res2 <- rma(yi, vi, scale = ~ 1, data=dat, test="t", control=list(optimizer="optim"))
   res3 <- suppressWarnings(rma(yi, vi, scale = ~ 1, link="identity", data=dat, test="t", control=list(optimizer="optim", optmethod="Nelder-Mead")))
   expect_equivalent(round(res1$tau2, 3), round(as.vector(exp(coef(res2)$alpha)), 3))
   expect_equivalent(round(res1$tau2, 3), round(as.vector(coef(res3)$alpha), 3))

})

test_that("location-scale model results are correct for a categorical predictor", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res1 <- rma(yi ~ alloc, vi, scale = ~ alloc - 1, data=dat)
   res2 <- rma(yi ~ alloc, vi, scale = ~ alloc - 1, link = "identity", data=dat)
   res3 <- rma.mv(yi ~ alloc, vi, random = ~ alloc | trial, struct="DIAG", data=dat)
   expect_equivalent(round(as.vector(exp(coef(res1)$alpha)), 3), round(as.vector(coef(res2)$alpha), 3))
   expect_equivalent(round(as.vector(exp(coef(res1)$alpha)), 3), round(res3$tau2, 3))

})

test_that("location-scale model results are correct for a continuous predictor", {

   data(dat.laopaiboon2015, package="metafor")
   dat <- escalc(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.laopaiboon2015)

   dat$ni <- dat$n1i + dat$n2i
   dat$ni[dat$study == "Whitlock"] <- dat$ni[dat$study == "Whitlock"] + 2

   res <- suppressWarnings(rma(yi, vi, scale = ~ I(1/ni) - 1, link="identity", data=dat, method="ML"))
   expect_equivalent(round(as.vector(coef(res)$alpha), 3), 79.108)
   expect_equivalent(round(exp(c(res$beta, res$ci.lb, res$ci.ub)), 3), c(0.854, 0.548, 1.330))

   res <- rma(yi, vi, scale = ~ I(1/ni), link="identity", data=dat, method="ML")
   expect_equivalent(round(as.vector(coef(res)$alpha), 3), c(0.275, 31.513))
   expect_equivalent(round(exp(c(res$beta, res$ci.lb, res$ci.ub)), 3), c(1.016, 0.621, 1.662))

   res <- rma(yi, vi, scale = ~ I(1/ni) - 1, data=dat)
   expect_equivalent(round(as.vector(coef(res)$alpha), 3), -34.519)
   expect_equivalent(round(exp(c(res$beta, res$ci.lb, res$ci.ub)), 3), c(1.125, 0.638, 1.984))

   res <- rma(yi, vi, scale = ~ I(1/ni), data=dat)
   expect_equivalent(round(as.vector(coef(res)$alpha), 3), c(-0.887, 42.407))
   expect_equivalent(round(exp(c(res$beta, res$ci.lb, res$ci.ub)), 3), c(1.047, 0.624, 1.758))

})
