### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma.mv() function")

source("settings.r")

dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

test_that("rma.mv() correctly handles a formula for the 'yi' argument", {

   res1 <- rma.mv(yi ~ ablat, vi, random = ~ 1 | trial, data=dat, sparse=sparse)
   res2 <- rma.mv(yi, vi, mods = ~ ablat, random = ~ 1 | trial, data=dat, sparse=sparse)
   expect_equivalent(coef(res1), coef(res2), tolerance=.tol[["coef"]])

})

test_that("rma.mv() works correctly when using user-defined weights", {

   res <- rma.mv(yi, vi, W=1, random = ~ 1 | trial, data=dat, sparse=sparse)
   expect_equivalent(coef(res), mean(dat$yi), tolerance=.tol[["coef"]])
   expect_equivalent(c(vcov(res)), 0.0358, tolerance=.tol[["var"]])

})

test_that("rma.mv() correctly handles negative sampling variances", {

   dat$vi[1] <- -.01
   expect_warning(res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, sparse=sparse))
   expect_equivalent(coef(res), -0.7220, tolerance=.tol[["coef"]])
   expect_equivalent(c(vcov(res)), 0.0293, tolerance=.tol[["var"]])

})

test_that("rma.mv() correctly handles a missing value", {

   dat$vi[1] <- NA
   expect_warning(res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, sparse=sparse))
   expect_equivalent(coef(res), -0.7071, tolerance=.tol[["coef"]])
   expect_equivalent(c(vcov(res)), 0.0361, tolerance=.tol[["var"]])

})

test_that("rma.mv() correctly handles the R argument", {

   P <- structure(c(1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                    0.000, 1.000, 0.621, 0.621, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128,
                    0.000, 0.621, 1.000, 0.642, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128,
                    0.000, 0.621, 0.642, 1.000, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128, 0.128,
                    0.000, 0.128, 0.128, 0.128, 1.000, 0.266, 0.266, 0.221, 0.221, 0.221, 0.157, 0.157, 0.157, 0.157, 0.157,
                    0.000, 0.128, 0.128, 0.128, 0.266, 1.000, 0.467, 0.221, 0.221, 0.221, 0.157, 0.157, 0.157, 0.157, 0.157,
                    0.000, 0.128, 0.128, 0.128, 0.266, 0.467, 1.000, 0.221, 0.221, 0.221, 0.157, 0.157, 0.157, 0.157, 0.157,
                    0.000, 0.128, 0.128, 0.128, 0.221, 0.221, 0.221, 1.000, 0.605, 0.296, 0.157, 0.157, 0.157, 0.157, 0.157,
                    0.000, 0.128, 0.128, 0.128, 0.221, 0.221, 0.221, 0.605, 1.000, 0.296, 0.157, 0.157, 0.157, 0.157, 0.157,
                    0.000, 0.128, 0.128, 0.128, 0.221, 0.221, 0.221, 0.296, 0.296, 1.000, 0.157, 0.157, 0.157, 0.157, 0.157,
                    0.000, 0.128, 0.128, 0.128, 0.157, 0.157, 0.157, 0.157, 0.157, 0.157, 1.000, 0.773, 0.390, 0.390, 0.390,
                    0.000, 0.128, 0.128, 0.128, 0.157, 0.157, 0.157, 0.157, 0.157, 0.157, 0.773, 1.000, 0.390, 0.390, 0.390,
                    0.000, 0.128, 0.128, 0.128, 0.157, 0.157, 0.157, 0.157, 0.157, 0.157, 0.390, 0.390, 1.000, 0.606, 0.606,
                    0.000, 0.128, 0.128, 0.128, 0.157, 0.157, 0.157, 0.157, 0.157, 0.157, 0.390, 0.390, 0.606, 1.000, 0.697,
                    0.000, 0.128, 0.128, 0.128, 0.157, 0.157, 0.157, 0.157, 0.157, 0.157, 0.390, 0.390, 0.606, 0.697, 1.000),
                    .Dim = c(15L, 15L), .Dimnames = list(c("S11", "S15", "S06", "S10", "S08", "S02", "S07", "S14", "S09", "S01", "S12", "S05", "S13", "S04", "S03"),
                    c("S11", "S15", "S06", "S10", "S08", "S02", "S07", "S14", "S09", "S01", "S12", "S05", "S13", "S04", "S03")))

   dat <- structure(list(study = 1:44, species = c("S01", "S01", "S02", "S02", "S02", "S02", "S03", "S03", "S03", "S03", "S04", "S04", "S04", "S04", "S05", "S05", "S05", "S06", "S06", "S06", "S06", "S07", "S07", "S08", "S08", "S08", "S09", "S09", "S10", "S10", "S10", "S11", "S11", "S11", "S11", "S12", "S12", "S13", "S13", "S13", "S14", "S14", "S15", "S15"), phylogeny = c("S01", "S01", "S02", "S02", "S02", "S02", "S03", "S03", "S03", "S03", "S04", "S04", "S04", "S04", "S05", "S05", "S05", "S06", "S06", "S06", "S06", "S07", "S07", "S08", "S08", "S08", "S09", "S09", "S10", "S10", "S10", "S11", "S11", "S11", "S11", "S12", "S12", "S13", "S13", "S13", "S14", "S14", "S15", "S15"),
                         yi = c(1.91, 1.67, -0.92, -0.1, -0.58, -1.29, 0.04, -1.33, 0.02, -1, 0.2, 1.75, -0.75, 1.36, 1.24, 0.64, 0.52, 1.93, 1.11, 1.12, 1.17, 0.25, 1.95, -0.06, -0.79, 0.39, 1.61, 1.96, 0.93, 0.5, 0.73, -0.7, 0.11, 0.84, 1.83, -0.59, 0.19, 0.14, 0.74, 0.55, 0.34, -1.16, 1.93, 1.85),
                         vi = c(0.213, 0.387, 0.381, 0.467, 0.132, 0.603, 0.374, 0.2, 0.119, 0.092, 0.139, 0.449, 0.412, 0.398, 0.25, 0.168, 0.303, 0.125, 0.164, 0.229, 0.482, 0.059, 0.421, 0.111, 0.373, 0.032, 0.062, 0.126, 0.066, 0.155, 0.229, 0.276, 0.039, 0.409, 0.312, 0.304, 0.601, 0.096, 0.216, 0.181, 0.537, 0.16, 0.303, 0.281)),
                         .Names = c("study", "species", "phylogeny", "yi", "vi"), row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44"), class = "data.frame")

   res <- rma.mv(yi, vi, random = list(~ 1 | study, ~ 1 | species, ~ 1 | phylogeny), R = list(phylogeny=P), data=dat, sparse=sparse)

   expect_equivalent(coef(res), .5504, tolerance=.tol[["coef"]])
   expect_equivalent(res$sigma2, c(0.1763, 0.5125, 0.1062), tolerance=.tol[["var"]])
   expect_equivalent(c(logLik(res)), -54.6272, tolerance=.tol[["fit"]])

})

dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

test_that("rma.mv() correctly computes the Hessian", {

   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, cvvc=TRUE, sparse=sparse)
   expect_equivalent(c(sqrt(res$vvc)), 0.1678, tolerance=.tol[["se"]])

})

test_that("rma.mv() works correctly with test='t'", {

   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, test="t", sparse=sparse)
   expect_equivalent(res$pval, 0.0018, tolerance=.tol[["pval"]])

})

test_that("rma.mv() works correctly with different optimizers", {

   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="BFGS"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="L-BFGS-B"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="Nelder-Mead"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3133, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="nlminb"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="uobyqa"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="newuoa"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="bobyqa"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="nloptr"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="nlm"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="hjk"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="nmk"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3131, tolerance=.tol[["var"]])
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, control=list(optimizer="ucminf"), sparse=sparse)
   expect_equivalent(res$sigma2, 0.3132, tolerance=.tol[["var"]])

})

test_that("rma.mv() correctly handles 'beta' argument", {

   dat <- dat.berkey1998
   V <- vcalc(vi=1, cluster=author, rvars=c(v1i, v2i), data=dat)
   res.unc <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")
   res.01 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML", beta=c(0,0))
   res.02 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML", beta=c(NA,0))
   res.03 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML", beta=c(0,NA))

   fstats <- fitstats(res.01, res.02, res.03, res.unc)
   expect_equivalent(unlist(fstats[1,]), c(-2.464111, -0.691524, 1.010033, 5.840657), tolerance=.tol[["fit"]])

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
   res <- rma(yi, vi, scale = ~ 1, data=dat, optbeta=TRUE, beta=0)
   ll1 <- logLik(res)
   res <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, beta=0)
   ll2 <- logLik(res)

   expect_equivalent(ll1, ll2, tolerance=.tol[["fit"]])

})

rm(list=ls())
