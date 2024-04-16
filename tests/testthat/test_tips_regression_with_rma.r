### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking tip: rma() results match up with those from lm()")

source("settings.r")

### this is essentially checking the equivalence of the results as explained here:
### https://www.metafor-project.org/doku.php/tips:regression_with_rma

test_that("results for rma() and lm() match for method='FE'.", {

   stackloss$vi <- 0

   res.lm  <- lm(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., data=stackloss)
   res.rma <- rma(stack.loss, vi, mods =  ~ Air.Flow + Water.Temp + Acid.Conc., data=stackloss, test="knha", control=list(REMLf=FALSE))

   ### log likelihood (REML) should be the same
   expect_equivalent(logLik(res.lm, REML=TRUE), logLik(res.rma), tolerance=.tol[["fit"]])

   ### coefficients should be the same
   expect_equivalent(coef(res.lm), coef(res.rma), tolerance=.tol[["coef"]])

   ### var-cov matrix should be the same
   expect_equivalent(matrix(vcov(res.lm), nrow=4, ncol=4), matrix(vcov(res.rma), nrow=4, ncol=4), tolerance=.tol[["var"]])

   ### fitted values should be the same
   expect_equivalent(fitted(res.lm), fitted(res.rma), tolerance=.tol[["pred"]])

   ### standardized residuals should be the same
   expect_equivalent(rstandard(res.lm), rstandard(res.rma)$z, tolerance=.tol[["test"]])

   ### studentized residuals should be the same
   expect_equivalent(rstudent(res.lm), rstudent(res.rma)$z, tolerance=.tol[["test"]])

   ### hat values should be the same
   expect_equivalent(hatvalues(res.lm), hatvalues(res.rma), tolerance=.tol[["inf"]])

   ### dffits should be the same
   expect_equivalent(dffits(res.lm), influence(res.rma)$inf$dffits, tolerance=.tol[["inf"]])

   ### covratios should be the same
   expect_equivalent(covratio(res.lm), influence(res.rma)$inf$cov.r, tolerance=.tol[["inf"]])

   ### dfbetas should be the same
   expect_equivalent(as.matrix(dfbetas(res.lm)), as.matrix(dfbetas(res.rma)), tolerance=.tol[["inf"]])

   ### Cook's distancs should differ by a factor of p
   expect_equivalent(cooks.distance(res.lm), cooks.distance(res.rma)/res.rma$p, tolerance=.tol[["inf"]])

})

rm(list=ls())
