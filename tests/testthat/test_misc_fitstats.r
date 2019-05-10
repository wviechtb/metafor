### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: computations of fit statistics")

source("tolerances.r") # read in tolerances

test_that("fit statistics are correct for rma.uni().", {

   ### load BCG dataset
   data(dat.bcg, package="metafor")

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### fit random- and mixed-effects models (with ML estimation)
   res1 <- rma(yi, vi, data=dat, method="ML")
   res2 <- rma(yi ~ ablat, vi, data=dat, method="ML")

   tmp <- c(logLik(res1))
   expect_equivalent(tmp, -12.6651, tolerance=.tol[["fit"]])
   expect_equivalent(tmp, sum(dnorm(dat$yi, coef(res1), sqrt(dat$vi+res1$tau2), log=TRUE)), tolerance=.tol[["fit"]])

   tmp <- deviance(res1)
   expect_equivalent(tmp, 37.1160, tolerance=.tol[["fit"]])
   expect_equivalent(tmp, -2 * (sum(dnorm(dat$yi, coef(res1), sqrt(dat$vi+res1$tau2), log=TRUE)) - sum(dnorm(dat$yi, dat$yi, sqrt(dat$vi), log=TRUE))), tolerance=.tol[["fit"]])

   tmp <- AIC(res1)
   expect_equivalent(tmp, 29.3302, tolerance=.tol[["fit"]])
   expect_equivalent(tmp, -2 * sum(dnorm(dat$yi, coef(res1), sqrt(dat$vi+res1$tau2), log=TRUE)) + 2*2, tolerance=.tol[["fit"]])

   tmp <- AIC(res1, res2)
   expect_equivalent(tmp, structure(list(df = c(2, 3), AIC = c(29.3302, 21.3713)), .Names = c("df", "AIC"), row.names = c("res1", "res2"), class = "data.frame"), tolerance=.tol[["fit"]])

   tmp <- BIC(res1)
   expect_equivalent(tmp, 30.4601, tolerance=.tol[["fit"]])
   expect_equivalent(tmp, -2 * sum(dnorm(dat$yi, coef(res1), sqrt(dat$vi+res1$tau2), log=TRUE)) + 2*log(res1$k), tolerance=.tol[["fit"]])

   tmp <- BIC(res1, res2)
   expect_equivalent(tmp, structure(list(df = c(2, 3), BIC = c(30.4601, 23.0662)), .Names = c("df", "BIC"), row.names = c("res1", "res2"), class = "data.frame"), tolerance=.tol[["fit"]])

   tmp <- c(fitstats(res1))
   expect_equivalent(tmp, c(-12.6651, 37.1160, 29.3302, 30.4601, 30.5302), tolerance=.tol[["fit"]])

   tmp <- fitstats(res1, res2)
   expect_equivalent(tmp, structure(list(res1 = c(-12.6651, 37.116, 29.3302, 30.4601, 30.5302), res2 = c(-7.6857, 27.1572, 21.3713, 23.0662, 24.038)), .Names = c("res1", "res2"), row.names = c("logLik:", "deviance:", "AIC:", "BIC:", "AICc:"), class = "data.frame"), tolerance=.tol[["fit"]])

   tmp <- nobs(res1)
   expect_equivalent(tmp, 13)

   tmp <- df.residual(res1)
   expect_equivalent(tmp, 12)

})

test_that("fit statistics are correct for rma.mv().", {

   ### load data
   dat <- dat.berkey1998

   ### construct variance-covariance matrix of the observed outcomes
   V <- bldiag(lapply(split(dat[,c("v1i", "v2i")], dat$trial), as.matrix))

   ### multiple outcomes random-effects model (with ML estimation)
   res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")

   tmp <- c(logLik(res))
   expect_equivalent(tmp, 5.8407, tolerance=.tol[["fit"]])

   tmp <- deviance(res)
   expect_equivalent(tmp, 25.6621, tolerance=.tol[["fit"]])

   tmp <- AIC(res)
   expect_equivalent(tmp, -1.6813, tolerance=.tol[["fit"]])
   expect_equivalent(tmp, -2 * c(logLik(res)) + 2*5, tolerance=.tol[["fit"]])

   tmp <- BIC(res)
   expect_equivalent(tmp, -0.1684, tolerance=.tol[["fit"]])
   expect_equivalent(tmp, -2 * c(logLik(res)) + 5*log(res$k), tolerance=.tol[["fit"]])

   tmp <- c(fitstats(res))
   expect_equivalent(tmp, c(5.8407, 25.6621, -1.6813, -0.1684, 13.3187), tolerance=.tol[["fit"]])

   tmp <- nobs(res)
   expect_equivalent(tmp, 10)

   tmp <- df.residual(res)
   expect_equivalent(tmp, 8)

})
