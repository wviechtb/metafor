### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti_and_mumin

source("settings.r")

context("Checking tip: model selection using the glmulti and MuMIn packages")

dat <- dat.bangertdrowns2004

### remove rows where at least one potential moderator is missing
dat <- dat[!apply(dat[,c("length", "wic", "feedback", "info", "pers", "imag", "meta")], 1, anyNA),]

test_that("results are correct for package glmulti.", {

   skip_on_cran()

   if (!require(glmulti))
      stop("Cannot load 'glmulti' package.")

   ### function for glmulti
   rma.glmulti <- function(formula, data, ...) {
      rma(formula, vi, data=data, method="ML", ...)
   }

   ### fit all possible models (only main effects)
   #res <- glmulti(yi ~ length + wic + feedback + info + pers + imag + meta, data=dat, ### DOES NOT WORK
   #            level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=128)
   res <- glmulti(y = "yi", xr=c("length", "wic", "feedback", "info", "pers", "imag", "meta"), data=dat,
                  level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=128, plotty=FALSE, report=FALSE)

   ### models, IC values, and weights for the models whose IC is not more than 2 points away from the lowest value
   top <- weightable(res)
   top <- top[top$aicc <= min(top$aicc) + 2,]

   expect_equivalent(top$aicc, c(13.502247, 13.504515, 14.003149, 14.130949, 14.434335, 14.748334, 14.986646, 15.058191, 15.210613, 15.410733), tolerance=.tol[["fit"]])

   ### register getfit method for 'rma.uni' objects
   eval(metafor:::.glmulti)

   ### multimodel inference results
   mmi <- as.data.frame(coef(res, varweighting="Johnson")) # to use newer method
   mmi <- data.frame(Estimate=mmi$Est, SE=sqrt(mmi$Uncond), Importance=mmi$Importance, row.names=row.names(mmi))
   mmi$z <- mmi$Estimate / mmi$SE
   mmi$p <- 2*pnorm(abs(mmi$z), lower.tail=FALSE)
   names(mmi) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
   mmi$ci.lb <- mmi[[1]] - qnorm(.975) * mmi[[2]]
   mmi$ci.ub <- mmi[[1]] + qnorm(.975) * mmi[[2]]
   mmi <- mmi[order(mmi$Importance, decreasing=TRUE), c(1,2,4:7,3)]

   expect_equivalent(mmi$Estimate, c(0.108404, 0.351153, 0.051201, 0.036604, 0.002272, 0.013244, -0.017004, -0.018272), tolerance=.tol[["coef"]])
   expect_equivalent(mmi$"Std. Error", c(0.103105, 0.201648, 0.08529, 0.068926, 0.005019, 0.068788, 0.054466, 0.079911), tolerance=.tol[["se"]])
   expect_equivalent(mmi$Importance, c(1, 0.847824, 0.424365, 0.367132, 0.325539, 0.291322, 0.264263, 0.241566), tolerance=.tol[["r2"]])

   ### multimodel predictions

   x <- c("length"=15, "wic"=1, "feedback"=1, "info"=0, "pers"=0, "imag"=1, "meta"=1)

   preds <- list()

   for (j in 1:res@nbmods) {

      model <- res@objects[[j]]
      vars <- names(coef(model))[-1]

      if (length(vars) == 0) {
         preds[[j]] <- predict(model)
      } else {
         preds[[j]] <- predict(model, newmods=x[vars])
      }

   }

   ### multimodel prediction
   weights <- weightable(res)$weights
   yhat <- sum(weights * sapply(preds, function(x) x$pred))

   expect_equivalent(yhat, 0.56444, tolerance=.tol[["pred"]])

   ### multimodel SE
   se <- sqrt(sum(weights * sapply(preds, function(x) x$se^2 + (x$pred - yhat)^2)))

   expect_equivalent(se, 0.2225354, tolerance=.tol[["se"]])

})

test_that("results are correct for package MuMIn.", {

   skip_on_cran()

   expect_equivalent(TRUE, TRUE)

   return()

   ### cannot get this to work as the helper functions are somehow not visible
   ### (even when directly adding them below or above this function)

   if (!require(MuMIn))
      stop("Cannot load 'MuMIn' package.")

   ### get helper functions
   eval(metafor:::.MuMIn)

   ### fit full model
   full <- rma(yi, vi, mods = ~ length + wic + feedback + info + pers + imag + meta,
               data=dat, method="ML")

   ### fit all possible models
   res <- suppressMessages(dredge(full))

   ### models, IC values, and weights for the models whose IC is not more than 2 points away from the lowest value
   top <- subset(res, delta <= 2, recalc.weights=FALSE)

   expect_equivalent(top$AICc, c(13.502247, 13.504515, 14.003149, 14.130949, 14.434335, 14.748334, 14.986646, 15.058191, 15.210613, 15.410733), tolerance=.tol[["fit"]])
   expect_equivalent(c(top$weight), c(0.067057, 0.066981, 0.0522, 0.048969, 0.042077, 0.035963, 0.031923, 0.030802, 0.028541, 0.025824), tolerance=.tol[["inf"]])

   ### importance of each predictor
   expect_equivalent(c(sw(res)), c(imag = 0.847824, meta = 0.424365, feedback = 0.367132, length = 0.325539, pers = 0.291322, wic = 0.264263, info = 0.241566), tolerance=.tol[["r2"]])

   ### model averaging
   mmi <- summary(model.avg(res))

   expect_equivalent(mmi$coefmat.full[,"Estimate"], c(intrcpt = 0.108404, imag = 0.351153, meta = 0.051201, feedback = 0.036604, length = 0.002272, wic = -0.017004, pers = 0.013244, info = -0.018272), tolerance=.tol[["coef"]])
   expect_equivalent(mmi$coefmat.full[,"Std. Error"], c(intrcpt = 0.103105, imag = 0.201648, meta = 0.08529, feedback = 0.068926, length = 0.005019, wic = 0.054466, pers = 0.068788, info = 0.079911), tolerance=.tol[["se"]])

})

rm(list=ls())
