### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

source("settings.r")

context("Checking misc: cumul() functions")

test_that("cumul() works correctly for 'rma.uni' object.", {

   ### calculate log risk ratios and corresponding sampling variances
   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
                 data=dat.bcg, slab=paste0(author, ", ", year))

   ### fit random-effects model
   res <- rma(yi, vi, data=dat)

   ### cumulative meta-analysis (in the order of publication year)
   out <- cumul(res, order=year)

   expect_equivalent(out$estimate, c(-0.889311, -1.325003, -0.97208, -1.001094, -1.101424, -0.973464, -0.901251, -0.788566, -0.865607, -0.785211, -0.708206, -0.794768, -0.714532), tolerance=.tol[["est"]])

   ### with transformation
   out <- cumul(res, order=year, transf=exp)

   expect_equivalent(out$estimate, c(0.410939, 0.265802, 0.378296, 0.367477, 0.332398, 0.377772, 0.406061, 0.454496, 0.420796, 0.456024, 0.492527, 0.451686, 0.489421), tolerance=.tol[["est"]])

   ### add studies with the same publication year simultaneously
   out <- cumul(res, order=year, transf=exp, collapse=TRUE)

   expect_equivalent(out$estimate, c(0.410939, 0.265802, 0.378296, 0.367477, 0.332398, 0.377772, 0.406061, 0.420796, 0.456024, 0.492527, 0.451686, 0.489421), tolerance=.tol[["est"]])

})

test_that("cumul() works correctly for 'rma.mh' object.", {

   ### meta-analysis of the (log) risk ratios using the Mantel-Haenszel method
   res <- rma.mh(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg,
                 data=dat.bcg, slab=paste0(author, ", ", year))

   ### cumulative meta-analysis (in the order of publication year)
   out <- cumul(res, order=year)

   expect_equivalent(out$estimate, c(-0.889311, -1.351739, -0.827295, -0.837892, -0.894003, -0.850672, -0.835785, -0.774442, -0.789605, -0.666015, -0.635076, -0.775798, -0.45371), tolerance=.tol[["est"]])

   ### with transformation
   out <- cumul(res, order=year, transf=exp)

   expect_equivalent(out$estimate, c(0.410939, 0.25879, 0.437231, 0.432621, 0.409015, 0.427128, 0.433534, 0.460961, 0.454024, 0.513752, 0.529895, 0.460336, 0.635267), tolerance=.tol[["est"]])

   ### add studies with the same publication year simultaneously
   out <- cumul(res, order=year, transf=exp, collapse=TRUE)

   expect_equivalent(out$estimate, c(0.410939, 0.25879, 0.437231, 0.432621, 0.409015, 0.427128, 0.433534, 0.454024, 0.513752, 0.529895, 0.460336, 0.635267), tolerance=.tol[["est"]])

})

test_that("cumul() works correctly for 'rma.peto' object.", {

   ### meta-analysis of the (log) odds ratios using Peto's method
   res <- rma.peto(ai=tpos, bi=tneg, ci=cpos, di=cneg,
                   data=dat.bcg, slab=paste0(author, ", ", year))

   ### cumulative meta-analysis (in the order of publication year)
   out <- cumul(res, order=year)

   expect_equivalent(out$estimate, c(-0.860383, -1.240086, -0.957016, -0.964222, -1.000315, -0.940979, -0.924632, -0.850145, -0.871198, -0.726791, -0.691247, -0.816134, -0.474446), tolerance=.tol[["est"]])

   ### with transformation
   out <- cumul(res, order=year, transf=exp)

   expect_equivalent(out$estimate, c(0.423, 0.289359, 0.384037, 0.38128, 0.367764, 0.390246, 0.396677, 0.427353, 0.41845, 0.483458, 0.500951, 0.442138, 0.622229), tolerance=.tol[["est"]])

   ### add studies with the same publication year simultaneously
   out <- cumul(res, order=year, transf=exp, collapse=TRUE)

   expect_equivalent(out$estimate, c(0.423, 0.289359, 0.384037, 0.38128, 0.367764, 0.390246, 0.396677, 0.41845, 0.483458, 0.500951, 0.442138, 0.622229), tolerance=.tol[["est"]])

})

rm(list=ls())
