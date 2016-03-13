### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking tip: rma() results match up with those from lm()")

### this is essentially checking the equivalence of the results as explained here:
### http://www.metafor-project.org/doku.php/tips:regression_with_rma

test_that("results for rma() and lm() match for method='FE'.", {

   stackloss$vi <- 0

   res.lm  <- lm(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc., data=stackloss)
   res.rma <- rma(stack.loss, vi, mods =  ~ Air.Flow + Water.Temp + Acid.Conc., data=stackloss, knha=TRUE, control=list(REMLf=FALSE))

   ### log likelihood (REML) should be the same
   expect_equivalent(logLik(res.lm, REML=TRUE), logLik(res.rma))

   ### coefficients should be the same
   expect_equivalent(coef(res.lm), coef(res.rma))

   ### var-cov matrix should be the same
   expect_equivalent(matrix(vcov(res.lm), nrow=4, ncol=4), matrix(vcov(res.rma), nrow=4, ncol=4))

   ### fitted values should be the same
   expect_equivalent(fitted(res.lm), fitted(res.rma))

   ### standardized residuals should be the same
   expect_equivalent(rstandard(res.lm), rstandard(res.rma)$z)

   ### studentized residuals should be the same
   expect_equivalent(rstudent(res.lm), rstudent(res.rma)$z)

   ### hat values should be the same
   expect_equivalent(hatvalues(res.lm), hatvalues(res.rma))

   ### dffits should be the same
   expect_equivalent(dffits(res.lm), influence(res.rma)$inf$dffits)

   ### covratios should be the same
   expect_equivalent(covratio(res.lm), influence(res.rma)$inf$cov.r)

   ### dfbetas should be the same
   expect_equivalent(as.matrix(dfbetas(res.lm)), as.matrix(dfbetas(res.rma)))

   ### Cook's distancs should differ by a factor of p
   expect_equivalent(cooks.distance(res.lm), cooks.distance(res.rma)/res.rma$p)

})
