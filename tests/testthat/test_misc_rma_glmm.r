### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: rma.glmm() function")

dat <- get(data(dat.nielweise2007, package="metafor"))

test_that("rma.glmm() works correctly for 'UM.RS' model.", {

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.RS", method="FE"))
   out <- capture.output(print(res))

   expect_equivalent(round(coef(res), 3), -1.221)
   expect_equivalent(round(res$tau2, 3),   0)
   expect_equivalent(round(res$sigma2, 3), 0.616)

   expect_warning(res <- rma.glmm(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat, model="UM.RS", test="t"))
   out <- capture.output(print(res))

   expect_equivalent(round(coef(res), 3), -1.281)
   expect_equivalent(round(res$tau2, 3),   0.726)
   expect_equivalent(round(res$sigma2, 3), 0.521)

   ### check some (current) stop()'s

   expect_error(confint(res))
   expect_error(plot(res))
   expect_error(qqnorm(res))
   expect_error(weights(res))

})
