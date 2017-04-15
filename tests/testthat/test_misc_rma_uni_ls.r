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

   sav <- coef(summary(res))
   sav <- lapply(sav, round, 3)

   expected <- structure(list(beta = structure(list(estimate = 0.046, se = 0.264, zval = 0.175, pval = 0.861, ci.lb = -0.471, ci.ub = 0.564),
                  .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub"), row.names = "intrcpt", class = "data.frame"),
                  alpha = structure(list(estimate = c(-0.887, 42.407), se = c(1.239, 118.693), zval = c(-0.716, 0.357), pval = c(0.474, 0.721), ci.lb = c(-3.316, -190.228), ci.ub = c(1.542, 275.041)),
                  .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub"), row.names = c("intrcpt", "I(1/ni)"), class = "data.frame")),
                  .Names = c("beta", "alpha"))

   expect_equivalent(sav, expected)

   sav <- model.matrix(res)$scale
   expect_equivalent(sav, cbind(1, 1/dat$ni))

   sav <- fitted(res)$scale
   expect_equivalent(round(sav,3), c(-0.479, -0.588, -0.831, -0.711, -0.494, -0.254, -0.661, -0.458, -0.542, -0.039, -0.039, -0.13, -0.405, -0.764, -0.357))

})
