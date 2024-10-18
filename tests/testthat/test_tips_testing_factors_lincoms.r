### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

### see: https://www.metafor-project.org/doku.php/tips:testing_factors_lincoms

context("Checking tip: testing factors and linear combinations of parameters")

source("settings.r")

dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat

dat$year  <- dat$year  - 1948
dat$ablat <- dat$ablat - 13

test_that("results are correct when testing factors.", {

   res <- rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat)
   sav <- anova(res, btt=2:3)
   expect_equivalent(sav$QM, 1.366284, tolerance=.tol[["test"]])

   ### use linearHypothesis() from 'car' package for the same purpose

   if (!require(car))
      stop("Cannot load 'car' package.")

   sav2 <- linearHypothesis(res, rbind(c(0,1,0,0,0),c(0,0,1,0,0)))
   expect_equivalent(sav$QM, sav2$Chisq[2], tolerance=.tol[["test"]])

   sav3 <- linearHypothesis(res, c("factor(alloc)random = 0", "factor(alloc)systematic = 0"))
   expect_equivalent(sav$QM, sav3$Chisq[2], tolerance=.tol[["test"]])

   ### use glht() from 'multcomp' package for the same purpose

   if (!require(multcomp))
      stop("Cannot load 'multcomp' package.")

   sav4 <- summary(glht(res, linfct=rbind(b1=c(0,1,0,0,0), b2=c(0,0,1,0,0))), test=Chisqtest())
   expect_equivalent(sav$QM, sav4$test$SSH[1,1], tolerance=.tol[["test"]])

   ### show that reference level is not relevant

   res2 <- rma(yi, vi, mods = ~ relevel(factor(alloc), ref="random") + year + ablat, data=dat)
   sav5 <- anova(res2, btt=2:3)
   expect_equivalent(sav$QM, sav5$QM, tolerance=.tol[["test"]])

   ### likelihood ratio test

   res0 <- rma(yi, vi, mods = ~ year + ablat, data=dat, method="ML")
   res1 <- rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat, method="ML")
   sav <- anova(res0, res1)
   expect_equivalent(sav$LRT, 1.451038, tolerance=.tol[["test"]])

   ### Knapp & Hartung method

   res <- rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat, test="knha")
   sav <- anova(res, btt=2:3)
   expect_equivalent(sav$QM, 0.6502793, tolerance=.tol[["test"]])

   ### use linearHypothesis() from 'car' package for the same purpose

   sav2 <- linearHypothesis(res, c("factor(alloc)random = 0", "factor(alloc)systematic = 0"), test="F")
   expect_equivalent(sav$QM, sav2$F[2], tolerance=.tol[["test"]])

   ### use glht() from 'multcomp' package for the same purpose

   sav3 <- summary(glht(res, linfct=rbind(b1=c(0,1,0,0,0), b2=c(0,0,1,0,0))), test=Ftest())
   expect_equivalent(sav$QM, sav3$test$fstat[1,1], tolerance=.tol[["test"]])

})

test_that("results are correct when testing linear combinations.", {

   res <- rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat)

   sav1 <- anova(res, X=c(0,1,-1,0,0))
   sav2 <- linearHypothesis(res, c(0,1,-1,0,0))
   expect_equivalent(sav1$QM, sav2$Chisq[2], tolerance=.tol[["test"]])

   res <- rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat, test="knha")
   sav1 <- anova(res, X=c(0,1,-1,0,0))
   sav2 <- linearHypothesis(res, c(0,1,-1,0,0), test="F")
   expect_equivalent(sav1$QM, sav2$F[2], tolerance=.tol[["test"]])

   sav1 <- anova(res, X=c(1,1,0,1970-1948,30-13))
   sav2 <- linearHypothesis(res, c(1,1,0,1970-1948,30-13), test="F")
   expect_equivalent(sav1$QM, sav2$F[2], tolerance=.tol[["test"]])

   tmp <- predict(res, newmods=c(1,0,1970-1948,30-13))
   expect_equivalent(sav1$QM, (tmp$pred/tmp$se)^2, tolerance=.tol[["test"]])

})

test_that("results are correct when testing all pairwise comparisons.", {

   res <- rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat)

   sav1 <- anova(res, X=rbind(c(0,1,0,0,0), c(0,0,1,0,0), c(0,-1,1,0,0)))
   sav2 <- summary(glht(res, linfct=rbind(c(0,1,0,0,0), c(0,0,1,0,0), c(0,-1,1,0,0))), test=adjusted("none"))
   expect_equivalent(sav1$zval, sav2$test$tstat, tolerance=.tol[["test"]])

   sav1 <- confint(glht(res, linfct=rbind(c(0,1,0,0,0), c(0,0,1,0,0), c(0,-1,1,0,0))), calpha=univariate_calpha())
   sav2 <- predict(res, newmods=rbind(c(1,0,0,0), c(0,1,0,0), c(-1,1,0,0)), intercept=FALSE)
   expect_equivalent(sav1$confint[,2], sav2$ci.lb, tolerance=.tol[["ci"]])

   ### same results but leaving out the intercept

   res <- rma(yi, vi, mods = ~ 0 + factor(alloc) + year + ablat, data=dat)

   sav1 <- anova(res, X=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0)))
   sav2 <- anova(res, X=pairmat(btt=1:3))
   expect_equivalent(sav1$zval, sav2$zval, tolerance=.tol[["test"]])

   sav3 <- anova(res, X=pairmat(btt=1:3), adjust="holm")
   expect_equivalent(sav3$pval, c(0.882646, 0.981965, 0.882646), tolerance=.tol[["pval"]])

   sav4 <- summary(glht(res, linfct=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0))), test=adjusted("none"))
   expect_equivalent(sav1$zval, sav4$test$tstat, tolerance=.tol[["test"]])

   sav1 <- confint(glht(res, linfct=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0))), calpha=univariate_calpha())
   sav2 <- predict(res, newmods=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0)))
   expect_equivalent(sav1$confint[,2], sav2$ci.lb, tolerance=.tol[["ci"]])

   sav1 <- anova(res, X=pairmat(btt=1:3))
   sav2 <- summary(glht(res, linfct=cbind(contrMat(c("alternate"=1,"random"=1,"systematic"=1), type="Tukey"), 0, 0)), test=adjusted("none"))
   expect_equivalent(sav1$zval, sav2$test$tstat, tolerance=.tol[["test"]])

   ### with Knapp & Hartung adjustment

   res <- rma(yi, vi, mods = ~ factor(alloc) + year + ablat, data=dat, test="knha")

   sav1 <- anova(res, X=rbind(c(0,1,0,0,0), c(0,0,1,0,0), c(0,-1,1,0,0)))
   sav2 <- anova(res, X=pairmat(btt=2:3, btt2=2:3))
   expect_equivalent(sav1$zval, sav2$zval, tolerance=.tol[["test"]])

   sav3 <- summary(glht(res, linfct=rbind(c(0,1,0,0,0), c(0,0,1,0,0), c(0,-1,1,0,0)), df=df.residual(res)), test=adjusted("none"))
   expect_equivalent(sav1$zval, sav3$test$tstat, tolerance=.tol[["test"]])

   sav1 <- confint(glht(res, linfct=rbind(c(0,1,0,0,0), c(0,0,1,0,0), c(0,-1,1,0,0)), df=df.residual(res)), calpha=univariate_calpha())
   sav2 <- predict(res, newmods=rbind(c(1,0,0,0), c(0,1,0,0), c(-1,1,0,0)), intercept=FALSE)
   expect_equivalent(sav1$confint[,2], sav2$ci.lb, tolerance=.tol[["ci"]])

   ### same results but leaving out the intercept

   res <- rma(yi, vi, mods = ~ 0 + factor(alloc) + year + ablat, data=dat, test="knha")

   sav1 <- anova(res, X=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0)))
   sav2 <- anova(res, X=pairmat(btt=1:3))
   expect_equivalent(sav1$zval, sav2$zval, tolerance=.tol[["test"]])

   sav3 <- summary(glht(res, linfct=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0)), df=df.residual(res)), test=adjusted("none"))
   expect_equivalent(sav1$zval, sav3$test$tstat, tolerance=.tol[["test"]])

   sav4 <- summary(glht(res, linfct=cbind(contrMat(c("alternate"=1,"random"=1,"systematic"=1), type="Tukey"), 0, 0), df=df.residual(res)), test=adjusted("none"))
   expect_equivalent(sav1$zval, sav4$test$tstat, tolerance=.tol[["test"]])

   sav1 <- confint(glht(res, linfct=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0)), df=df.residual(res)), calpha=univariate_calpha())
   sav2 <- predict(res, newmods=rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(0,-1,1,0,0)))
   expect_equivalent(sav1$confint[,2], sav2$ci.lb, tolerance=.tol[["ci"]])

   sav3 <- predict(res, newmods=pairmat(btt=1:3))
   expect_equivalent(sav2$ci.lb, sav3$ci.lb, tolerance=.tol[["ci"]])

})

rm(list=ls())
