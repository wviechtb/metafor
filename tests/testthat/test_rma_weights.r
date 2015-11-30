### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking that weights() functions provide the correct results")

test_that("weights are correct for rma() with method='FE'.", {

   data(dat.bcg, package="metafor")

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### weighted analysis
   res <- rma(yi, vi, data=dat, method="FE")

   ### weights should be the same as 1/vi (scaled to percentages)
   expect_equivalent(weights(res), (1/dat$vi)/sum(1/dat$vi) * 100)

   ### weighted analysis with user defined weights
   res <- rma(yi, vi, data=dat, method="FE", weights=1:13)

   ### weights should match (scaled to percentages)
   expect_equivalent(weights(res), (1:13)/sum(1:13) * 100)

   ### unweighted analysis
   res <- rma(yi, vi, data=dat, method="FE", weighted=FALSE)

   ### weights should be the same as 1/k (scaled to percentages)
   expect_equivalent(weights(res), rep(1/res$k, res$k) * 100)

   ### unweighted analysis (but user has specified weights nevertheless)
   res <- rma(yi, vi, data=dat, method="FE", weighted=FALSE, weights=1:13)

   ### weights should be the same as 1/k (scaled to percentages)
   expect_equivalent(weights(res), rep(1/res$k, res$k) * 100)

})

test_that("weights are correct for rma() with method='DL'.", {

   data(dat.bcg, package="metafor")

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### weighted analysis
   res <- rma(yi, vi, data=dat, method="DL")

   ### weights should be the same as 1/(vi+tau2) (scaled to percentages)
   expect_equivalent(weights(res), (1/(dat$vi+res$tau2)/sum(1/(dat$vi+res$tau2)) * 100))

   ### weighted analysis with user defined weights
   res <- rma(yi, vi, data=dat, method="DL", weights=1:13)

   ### weights should match (scaled to percentages)
   expect_equivalent(weights(res), (1:13)/sum(1:13) * 100)

   ### unweighted analysis
   res <- rma(yi, vi, data=dat, method="DL", weighted=FALSE)

   ### weights should be the same as 1/k (scaled to percentages)
   expect_equivalent(weights(res), rep(1/res$k, res$k) * 100)

   ### unweighted analysis (but user has specified weights nevertheless)
   res <- rma(yi, vi, data=dat, method="FE", weighted=FALSE, weights=1:13)

   ### weights should be the same as 1/k (scaled to percentages)
   expect_equivalent(weights(res), rep(1/res$k, res$k) * 100)

})
