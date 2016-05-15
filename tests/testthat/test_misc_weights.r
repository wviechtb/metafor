### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: weights() function")

test_that("weights are correct for rma() with method='FE'.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### weighted analysis
   res <- rma(yi, vi, data=dat, method="FE")

   ### weights should be the same as 1/vi (scaled to percentages)
   expect_equivalent(weights(res), (1/dat$vi)/sum(1/dat$vi) * 100)

   ### weights should be the same as 1/vi
   expect_equivalent(diag(weights(res, type="matrix")), 1/dat$vi)

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

   ### weights should be the same as 1/(vi+tau2)
   expect_equivalent(diag(weights(res, type="matrix")), 1/(dat$vi+res$tau2))

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

test_that("weights are correct for rma.mv() with method='REML'.", {

   data(dat.bcg, package="metafor")
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   ### weighted analysis
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat)

   ### weights should be the same as 1/(vi+sigma2) (scaled to percentages)
   expect_equivalent(weights(res), (1/(dat$vi+res$sigma2)/sum(1/(dat$vi+res$sigma2)) * 100))

   ### weights should be the same as 1/(vi+sigma2)
   expect_equivalent(diag(weights(res, type="matrix")), 1/(dat$vi+res$sigma2))

   ### weighted analysis with user defined weights
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, W=1:13)

   ### weights should match (scaled to percentages)
   expect_equivalent(weights(res), (1:13)/sum(1:13) * 100)

   ### unweighted analysis
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, W=1)

   ### weights should be the same as 1/k (scaled to percentages)
   expect_equivalent(weights(res), rep(1/res$k, res$k) * 100)

})

test_that("weights are correct for rma.mh() with measure='RD/RR/OR'.", {

   dat <- get(data(dat.bcg, package="metafor"))

   res <- rma.mh(measure="RD", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)
   sav <- weights(res)
   expect_equivalent(round(coef(res), 6), round(sum(res$yi * sav/100), 6))
   tmp <- diag(weights(res, type="matrix"))
   expect_equivalent(sav, tmp/sum(tmp)*100)

   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)
   sav <- weights(res)
   expect_equivalent(round(exp(coef(res)), 6), round(sum(exp(res$yi) * sav/100), 6))
   tmp <- diag(weights(res, type="matrix"))
   expect_equivalent(sav, tmp/sum(tmp)*100)

   res <- rma.mh(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)
   sav <- weights(res)
   expect_equivalent(round(exp(coef(res)), 6), round(sum(exp(res$yi) * sav/100), 6))
   tmp <- diag(weights(res, type="matrix"))
   expect_equivalent(sav, tmp/sum(tmp)*100)

})

test_that("weights are correct for rma.mh() with measure='IRD/IRR'.", {

   dat <- get(data(dat.nielweise2008, package="metafor"))

   res <- rma.mh(measure="IRD", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat)
   sav <- weights(res)
   expect_equivalent(round(coef(res), 6), round(sum(res$yi * sav/100), 6))
   tmp <- diag(weights(res, type="matrix"))
   expect_equivalent(sav, tmp/sum(tmp)*100)

   res <- rma.mh(measure="IRR", x1i=x1i, t1i=t1i, x2i=x2i, t2i=t2i, data=dat)
   sav <- weights(res)
   expect_equivalent(round(exp(coef(res)), 6), round(sum(exp(res$yi) * sav/100), 6))
   tmp <- diag(weights(res, type="matrix"))
   expect_equivalent(sav, tmp/sum(tmp)*100)

})

test_that("weights are correct for rma.peto().", {

   dat <- get(data(dat.bcg, package="metafor"))

   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)
   sav <- weights(res)
   expect_equivalent(round(coef(res), 6), round(sum(res$yi * sav/100), 6))
   tmp <- diag(weights(res, type="matrix"))
   expect_equivalent(sav, tmp/sum(tmp)*100)

})
