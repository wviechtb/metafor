### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:viechtbauer2007b

context("Checking analysis example: viechtbauer2007b")

### create dataset for example
data(dat.linde2005, package="metafor")
dat <- escalc(measure="RR", ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=dat.linde2005)
dat <- dat[c(7:10,13:25), c(13:16,18:19,11,6,7,9)]
dat$dosage <- (dat$dosage * 7) / 1000

test_that("results are correct for the CIs.", {

   sav <- summary(dat, transf=exp)[c(13,17),]

   ### compare with results on page 106
   tmp <- round(sav$ci.lb, 2)
   expect_equivalent(tmp, c(.74, 1.00)) ### 1.01 in article
   tmp <- round(sav$ci.ub, 2)
   expect_equivalent(tmp, c(1.28, 1.54))

})

test_that("results are correct for the fixed-effects model.", {

   res <- rma(yi, vi, data=dat, method="FE")
   sav <- predict(res, transf=exp)
   tmp <- round(c(sav$pred, sav$ci.lb, sav$ci.ub), 2)

   ### compare with results on page 107
   expect_equivalent(tmp, c(1.38, 1.26, 1.52)) ### 1.39 in article
   expect_equivalent(round(res$QE, 2), 51.55) ### 55.54 in article

})

test_that("results are correct for the random-effects model.", {

   res <- rma(yi, vi, data=dat, method="DL")
   sav <- predict(res, transf=exp)

   ### compare with results on page 109
   tmp <- round(c(sav$pred, sav$ci.lb, sav$ci.ub), 2)
   expect_equivalent(tmp, c(1.57, 1.31, 1.89)) ### 1.90 in article
   tmp <- round(c(sav$cr.lb, sav$cr.ub), 2)
   expect_equivalent(tmp, c(.85, 2.91)) ### .87, 2.83 in article (but this was calculated without taking Var[hat(mu)] into consideration)
   expect_equivalent(round(res$tau2, 3), .090) ### .091 in article

})

test_that("results are correct for the mixed-effects model.", {

   dat$dosage <- dat$dosage * dat$duration
   res <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="DL")

   ### compare with results on page 112
   expect_equivalent(round(res$tau2, 3), .047)
   expect_equivalent(round(res$R2, 2), 47.38) ### 48% in article

   sav <- structure(list(estimate = c(0.476, -0.006, -0.067, -0.002), se = c(0.088, 0.01, 0.035, 0.003),
                         zval = c(5.434, -0.585, -1.909, -0.456), pval = c(0, 0.559, 0.056, 0.649)),
                         .Names = c("estimate", "se", "zval", "pval"),
                         row.names = c("Intercept", "Dosage", "Baseline", "Dosage x Baseline"), class = "data.frame")

   ### compare with results in Table II on page 113
   expect_equivalent(round(coef(summary(res))[,1:4], 3), sav)

   ### compare with results on page 113
   sav <- predict(res, newmods=c(34-34, 12.5-20, (34-34)*(12.5-20)), transf=exp)
   tmp <- round(c(sav$pred, sav$ci.lb, sav$ci.ub), 2)
   expect_equivalent(tmp, c(2.67, 1.46, 4.88)) ### 2.66, 1.46, 4.90 in article
   sav <- predict(res, newmods=c(34-34, 23.6-20, (34-34)*(23.6-20)), transf=exp)
   tmp <- round(c(sav$pred, sav$ci.lb, sav$ci.ub), 2)
   expect_equivalent(tmp, c(1.26, 0.99, 1.61)) ### 1.61 in article

   skip_on_cran()

   size <- 1 / sqrt(dat$vi)
   size <- 0.15 * size / max(size)

   modvals <- cbind(0, cbind(seq(12, 24, by=.1)) - 20, 0)
   preds   <- predict(res, modvals, transf=exp)

   opar <- par(no.readonly=TRUE)

   plot(NA, NA, xlab="Baseline HRSD Score", ylab="Relative Rate", xlim=c(12,24), ylim=c(0.5,4.0), bty="l")
   abline(h=seq(1, 4, by=0.5), col="lightgray")
   abline(v=seq(14, 24, by=2), col="lightgray")
   lines(modvals[,2] + 20, preds$pred, col="darkgray", lwd=2)
   lines(modvals[,2] + 20, preds$ci.lb, col="darkgray", lty="dashed", lwd=2)
   lines(modvals[,2] + 20, preds$ci.ub, col="darkgray", lty="dashed", lwd=2)
   symbols(dat$baseline, exp(dat$yi), circles=size, inches=FALSE, add=TRUE, bg="black")

   par(opar)

   ### check results for all tau^2 estimators

   res.HS   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="HS")
   res.HE   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="HE")
   res.DL   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="DL")
   res.GENQ <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="GENQ", weights = n1i + n2i)
   res.SJ   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="SJ")
   res.DLIT <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="DLIT", control=list(maxiter=500))
   res.SJIT <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="SJIT")
   res.PM   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="PM")
   res.ML   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="ML")
   res.REML <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="REML")
   res.EB   <- rma(yi, vi, mods = ~ I(dosage-34) * I(baseline-20), data=dat, method="EB")

   res <- list(res.HS, res.HE, res.DL, res.GENQ, res.SJ, res.DLIT, res.SJIT, res.PM, res.ML, res.REML, res.EB)

   res <- data.frame(method=sapply(res, function(x) x$method),
                     tau2=sapply(res, function(x) round(x$tau2,3)),
                     se.tau2=sapply(res, function(x) round(x$se.tau2, 3)))

   expect_equivalent(res$tau2,    c(0.025, 0.039, 0.047, 0.060, 0.091, 0.030, 0.063, 0.063, 0.024, 0.056, 0.063))
   expect_equivalent(res$se.tau2, c(0.020, 0.076, 0.038, 0.053, 0.044, 0.044, 0.046, 0.046, 0.022, 0.041, 0.046))

})
