### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

### see: https://www.metafor-project.org/doku.php/tips:multiple_imputation_with_mice_and_metafor

source("settings.r")

context("Checking tip: multiple imputation with the mice and metafor packages")

dat <- dat.bangertdrowns2004

### keep only variables needed in the analysis
dat <- dat[c("yi", "vi", "length", "wic", "feedback", "info", "pers", "imag", "meta")]

test_that("results are correct for package mice.", {

   skip_on_cran()

   if (!require(mice))
      stop("Cannot load 'mice' package.")

   ### turn dummy variables into proper factors
   ### (so logistic regression is used for imputing missing values on these moderators)
   dat$wic      <- factor(dat$wic)
   dat$feedback <- factor(dat$feedback)
   dat$info     <- factor(dat$info)
   dat$pers     <- factor(dat$pers)
   dat$imag     <- factor(dat$imag)
   dat$meta     <- factor(dat$meta)

   ### create default predictor matrix
   predMatrix <- make.predictorMatrix(dat)

   ### adjust predictor matrix
   predMatrix[,"vi"] <- 0 # don't use vi for imputing
   predMatrix["yi",] <- 0 # don't impute yi (since yi has no NAs, this is actually irrelevant here)
   predMatrix["vi",] <- 0 # don't impute vi (since vi has no NAs, this is actually irrelevant here)

   ### create imputation methods vector
   impMethod <- make.method(dat)

   ### generate imputed datasets
   imp <- mice(dat, print=FALSE, m=20, predictorMatrix=predMatrix, method=impMethod, seed=1234)

   ### fit model to each completed dataset
   fit <- with(imp, rma(yi, vi, mods = ~ length + wic + feedback + info + pers + imag + meta))

   ### pool results
   pool <- summary(pool(fit))

   expect_equivalent(pool$estimate, c(0.381857, 0.013467, -0.056704, 0.000742, -0.30994, -0.309623, 0.201466, 0.448205), tolerance=.tol[["coef"]])
   expect_equivalent(pool$std.error, c(0.241471, 0.008772, 0.12978, 0.122219, 0.227993, 0.197197, 0.21107, 0.17783), tolerance=.tol[["se"]])
   expect_equivalent(pool$statistic, c(1.581377, 1.53531, -0.436921, 0.006069, -1.359426, -1.570121, 0.954497, 2.520414), tolerance=.tol[["test"]])
   expect_equivalent(pool$p.value, c(0.122388, 0.133287, 0.664798, 0.995194, 0.182211, 0.125322, 0.345947, 0.016331), tolerance=.tol[["pval"]])

})

test_that("results are correct for package Amelia.", {

   skip_on_cran()

   if (!require(Amelia))
      stop("Cannot load 'Amelia' package.")

   set.seed(1234)
   invisible(capture.output(imp <- amelia(dat, m=20, idvars=2, noms=4:9, incheck=TRUE, p2s=0)))

   fit <- lapply(imp$imputations, function(x) if (length(x)==1L) NULL else rma(yi, vi, mods = ~ length + wic + feedback + info + pers + imag + meta, data=x))
   fit <- fit[!sapply(fit, is.null)]

   b  <- sapply(fit, function(x) coef(x))
   se <- sapply(fit, function(x) x$se)
   pool <- mi.meld(b, se, byrow=FALSE)
   pool <- data.frame(estimate=pool$q.mi[1,], se=pool$se.mi[1,])

   expect_equivalent(pool$estimate, c(0.364127, 0.013386, -0.062038, 0.011519, -0.295951, -0.265128, 0.182751, 0.403387), tolerance=.tol[["coef"]])
   expect_equivalent(pool$se, c(0.240302, 0.008748, 0.133176, 0.121411, 0.229097, 0.200866, 0.212578, 0.18000), tolerance=.tol[["se"]])

})

rm(list=ls())
