### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: matreg() function")

source("settings.r")

test_that("matreg() works correctly for the 'mtcars' dataset.", {

   if (!require(MuMIn))
      stop("Cannot load 'MuMIn' package.")

   dat <- mtcars
   res1 <- lm(mpg ~ hp + wt + am, data=dat)
   S <- cov(dat)
   res2 <- matreg(y="mpg", x=c("hp","wt","am"), R=S, cov=TRUE, means=colMeans(dat), n=nrow(dat))
   expect_equivalent(coef(res1), coef(res2), tolerance=.tol[["coef"]])
   expect_equivalent(vcov(res1), vcov(res2), tolerance=.tol[["var"]])
   expect_equivalent(se(res1), se(res2), tolerance=.tol[["se"]])
   expect_equivalent(sigma(res1), sigma(res2), tolerance=.tol[["var"]])
   expect_equivalent(nobs(res1), nobs(res2), tolerance=.tol[["count"]])
   expect_equivalent(df.residual(res1), df.residual(res2), tolerance=.tol[["count"]])
   expect_equivalent(deviance(res1), deviance(res2), tolerance=.tol[["fit"]])
   expect_equivalent(logLik(res1), logLik(res2), tolerance=.tol[["fit"]])
   expect_equivalent(logLik(res1, REML=TRUE), logLik(res2, REML=TRUE), tolerance=.tol[["fit"]])
   expect_equivalent(AIC(res1), AIC(res2), tolerance=.tol[["fit"]])
   expect_equivalent(BIC(res1), BIC(res2), tolerance=.tol[["fit"]])
   expect_equivalent(AICc(res1), AICc(res2), tolerance=.tol[["fit"]])
   expect_equivalent(AICc(res2), AIC(res2, correct=TRUE), tolerance=.tol[["fit"]])
   expect_equivalent(confint(res1), confint(res2)$tab[,2:3], tolerance=.tol[["ci"]])

   res2 <- matreg(y="mpg", x=c("hp","wt","am"), R=S, cov=TRUE, n=nrow(dat))
   b1 <- coef(res1)
   b1[1] <- NA
   expect_equivalent(b1, coef(res2), tolerance=.tol[["coef"]])
   vb1 <- vcov(res1)
   vb1[1,] <- vb1[,1] <- NA
   expect_equivalent(vb1, vcov(res2), tolerance=.tol[["var"]])
   expect_equivalent(se(res1)[-1], se(res2)[-1], tolerance=.tol[["se"]])
   expect_equivalent(sigma(res1), sigma(res2), tolerance=.tol[["var"]])
   expect_equivalent(nobs(res1), nobs(res2), tolerance=.tol[["count"]])
   expect_equivalent(df.residual(res1), df.residual(res2), tolerance=.tol[["count"]])
   expect_equivalent(deviance(res1), deviance(res2), tolerance=.tol[["fit"]])
   expect_equivalent(logLik(res1), logLik(res2), tolerance=.tol[["fit"]])
   expect_equivalent(logLik(res1, REML=TRUE), logLik(res2, REML=TRUE), tolerance=.tol[["fit"]])
   expect_equivalent(AIC(res1)-2, AIC(res2), tolerance=.tol[["fit"]])
   expect_equivalent(BIC(res1)-log(nrow(dat)), BIC(res2), tolerance=.tol[["fit"]])
   expect_equivalent(AICc(res2), AIC(res2, correct=TRUE), tolerance=.tol[["fit"]])
   expect_equivalent(confint(res1)[-1,], confint(res2)$tab[-1,2:3], tolerance=.tol[["ci"]])

   dat[] <- scale(dat)
   res1 <- lm(mpg ~ 0 + hp + wt + am, data=dat)
   R <- cor(dat)
   res2 <- matreg(y="mpg", x=c("hp","wt","am"), R=R, n=nrow(dat))
   expect_equivalent(coef(res1), coef(res2), tolerance=.tol[["coef"]])
   expect_equivalent(vcov(res1), vcov(res2), tolerance=.tol[["var"]])
   expect_equivalent(se(res1), se(res2), tolerance=.tol[["se"]])
   expect_equivalent(sigma(res1), sigma(res2), tolerance=.tol[["var"]])
   expect_equivalent(nobs(res1), nobs(res2), tolerance=.tol[["count"]])
   expect_equivalent(df.residual(res1), df.residual(res2), tolerance=.tol[["count"]])
   expect_equivalent(deviance(res1), deviance(res2), tolerance=.tol[["fit"]])
   expect_equivalent(logLik(res1), logLik(res2), tolerance=.tol[["fit"]])
   expect_equivalent(logLik(res1, REML=TRUE), logLik(res2, REML=TRUE), tolerance=.tol[["fit"]])
   expect_equivalent(AIC(res1), AIC(res2), tolerance=.tol[["fit"]])
   expect_equivalent(BIC(res1), BIC(res2), tolerance=.tol[["fit"]])
   expect_equivalent(AICc(res1), AICc(res2), tolerance=.tol[["fit"]])
   expect_equivalent(AICc(res2), AIC(res2, correct=TRUE), tolerance=.tol[["fit"]])
   expect_equivalent(confint(res1), confint(res2)$tab[,2:3], tolerance=.tol[["ci"]])

})

test_that("matreg() works correctly for 'dat.craft2003'.", {

   dat <- dat.craft2003

   ### construct dataset and var-cov matrix of the correlations
   tmp <- rcalc(ri ~ var1 + var2 | study, ni=ni, data=dat)
   V <- tmp$V
   dat <- tmp$dat

   out <- capture.output(print(tmp))

   sav <- structure(list(study = c("1", "1", "1", "1", "1", "1"), var1 = c("acog", "asom", "conf", "acog", "acog", "asom"), var2 = c("perf", "perf", "perf", "asom", "conf", "conf"), var1.var2 = c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf"), yi = c(-0.55, -0.48, 0.66, 0.47, -0.38, -0.46), ni = c(142L, 142L, 142L, 142L, 142L, 142L)), row.names = c(NA, 6L), class = "data.frame")
   expect_equivalent(dat[1:6,], sav, tolerance=.tol[["coef"]])

   sav <- structure(c(0.00345039893617021, 0.00132651489361702, -0.000554579787234042, -0.00139678475177305, 0.00250189539007092, 0.000932237234042553, 0.00132651489361702, 0.00420059687943262, -0.000952140709219857, -0.00194335914893617, 0.00126485617021277, 0.00251607829787234, -0.000554579787234042, -0.000952140709219857, 0.00225920113475177, 0.00057910914893617, -0.00153379787234043, -0.00106924595744681, -0.00139678475177305, -0.00194335914893617, 0.00057910914893617, 0.00430494191489362, -0.00180268914893617, -0.00120505595744681, 0.00250189539007092, 0.00126485617021277, -0.00153379787234043, -0.00180268914893617, 0.00519185361702128, 0.00188440468085106, 0.000932237234042553, 0.00251607829787234, -0.00106924595744681, -0.00120505595744681, 0.00188440468085106, 0.00440833021276596), .Dim = c(6L, 6L), .Dimnames = list(c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf"), c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf")))
   expect_equivalent(V[1:6,1:6], sav, tolerance=.tol[["var"]])

   ### turn var1.var2 into a factor with the desired order of levels
   dat$var1.var2 <- factor(dat$var1.var2,
      levels=c("acog.perf", "asom.perf", "conf.perf", "acog.asom", "acog.conf", "asom.conf"))

   ### multivariate random-effects model
   expect_warning(res <- rma.mv(yi, V, mods = ~ 0 + var1.var2, random = ~ var1.var2 | study, struct="UN", data=dat, sparse=.sparse))

   ### restructure estimated mean correlations into a 4x4 matrix
   R <- matrix(NA, nrow=4, ncol=4)
   R[lower.tri(R)] <- coef(res)
   rownames(R) <- colnames(R) <- c("perf", "acog", "asom", "conf")

   ### fit regression model with 'perf' as outcome and 'acog', 'asom', and 'conf' as predictors
   fit <- matreg(1, 2:4, R=R, V=vcov(res))

   out <- capture.output(print(fit))

   sav <- structure(list(estimate = c(0.14817903234559, -0.0536342615587582, 0.363679177420187), se = c(0.156551433378687, 0.0768472434859867, 0.0909539697381244), zval = c(0.946519805967891, -0.697933447262015, 3.99849702511387), pval = c(0.343883525131896, 0.485218815885662, 0.0000637459821320369), ci.lb = c(-0.158656138804758, -0.204252091102472, 0.185412672482517), ci.ub = c(0.455014203495939, 0.0969835679849561, 0.541945682357857)), class = "data.frame", row.names = c("acog", "asom", "conf"))
   expect_equivalent(fit$tab, sav, tolerance=.tol[["misc"]])

   expect_equivalent(coef(fit), c(acog = 0.148179, asom = -0.053634, conf = 0.363679), tolerance=.tol[["coef"]])
   expect_equivalent(vcov(fit), structure(c(0.024508, -0.001433, 0.004527, -0.001433, 0.005905, 0.004089, 0.004527, 0.004089, 0.008273), dim = c(3L, 3L), dimnames = list(c("acog", "asom", "conf"), c("acog", "asom", "conf"))), tolerance=.tol[["var"]])
   expect_equivalent(se(fit), c(acog = 0.156552, asom = 0.076847, conf = 0.090954), tolerance=.tol[["se"]])
   expect_equivalent(confint(fit)$tab, structure(c(0.148179, -0.053634, 0.363679, -0.158656, -0.204252, 0.185413, 0.455015, 0.096984, 0.541945), dim = c(3L, 3L), dimnames = list(c("acog", "asom", "conf"), c("estimate", "ci.lb", "ci.ub"))), tolerance=.tol[["ci"]])

   ### use variable names
   fit <- matreg("perf", c("acog","asom","conf"), R=R, V=vcov(res))
   expect_equivalent(fit$tab, sav, tolerance=.tol[["misc"]])

})

rm(list=ls())
