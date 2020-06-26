### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:vanhouwelingen2002

context("Checking analysis example: vanhouwelingen2002")

source("tolerances.r") # read in tolerances

### load data
dat <- dat.colditz1994

### calculate log(OR)s and corresponding sampling variances
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### 'center' year variable
dat$year <- dat$year - 1900

test_that("results for the fixed-effects model are correct.", {

   res <- rma(yi, vi, data=dat, method="FE")
   tmp <- predict(res, transf=exp, digits=3)

   ### compare with results on page 596 (in text)
   expect_equivalent(tmp$pred,  .6465, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, .5951, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, .7024, tolerance=.tol[["ci"]]) ### .703 in paper

})

test_that("results for the random-effects model are correct.", {

   res <- rma(yi, vi, data=dat, method="ML")
   tmp <- predict(res, transf=exp, digits=3)

   ### compare with results on page 597 (in text)
   expect_equivalent(tmp$pred,  .4762, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, .3360, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, .6749, tolerance=.tol[["ci"]])
   expect_equivalent(res$tau2,  .3025, tolerance=.tol[["var"]])

   ### CI for tau^2
   tmp <- confint(res)
   expect_equivalent(tmp$random[1,2], 0.1302, tolerance=.tol[["var"]]) ### 0.135 based on a Satterthwaite approximation (page 597)
   expect_equivalent(tmp$random[1,3], 1.1812, tolerance=.tol[["var"]]) ### 1.181 based on a Satterthwaite approximation (page 597)

   ### CI for mu with Knapp & Hartung method
   res <- rma(yi, vi, data=dat, method="ML", test="knha")
   tmp <- predict(res, transf=exp, digits=3)

   ### (results for this not given in paper)
   expect_equivalent(tmp$ci.lb, .3175, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, .7141, tolerance=.tol[["ci"]])

})

test_that("profile plot for tau^2 can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma(yi, vi, data=dat, method="ML")

   opar <- par(no.readonly=TRUE)
   profile(res, xlim=c(.01,2), steps=200, log="x", cex=0, lwd=2, progbar=FALSE)
   abline(h=logLik(res) - 1.92, lwd=2)
   abline(v=c(.12, .89), lty="dashed")
   par(opar)

})

test_that("forest plot of observed log(OR)s and corresponding BLUPs can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma(yi, vi, data=dat, method="ML")
   sav <- blup(res)

   opar <- par(no.readonly=TRUE)
   par(family="mono", mar=c(5,4,1,2))
   forest(res, refline=res$beta, addpred=TRUE, xlim=c(-7,8), alim=c(-3,3), slab=1:13, psize=.8,
          ilab=paste0("(n = ", formatC(apply(dat[,c(4:7)], 1, sum), width=7, big.mark=","), ")"),
          ilab.xpos=-3.5, ilab.pos=2, rows=13:1+.15, header="Trial (total n)")
   arrows(sav$pi.lb, 13:1 - .15, sav$pi.ub, 13:1 -.15, length=.03, angle=90, code=3, lty="dotted")
   points(sav$pred, 13:1 - .15, pch=15, cex=.8)
   par(opar)

})

test_that("the prediction interval is correct.", {

   res <- rma(yi, vi, data=dat, method="ML")

   ### computation as done in the paper
   tmp <- c(res$beta) + c(-1,+1) * qnorm(.975) * sqrt(res$tau2)

   ### compare with results on page 599 (in text)
   expect_equivalent(tmp, c(-1.8199, 0.3359), tolerance=.tol[["ci"]])

   ### computation done with metafor
   tmp <- predict(res, digits=3)

   ### (results for this not given in paper)
   expect_equivalent(tmp$pi.lb, -1.875, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$pi.ub,  0.391, tolerance=.tol[["ci"]])

})

test_that("L'Abbe plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, method="FE")

   opar <- par(no.readonly=TRUE)
   labbe(res, xlim=c(-7,-1), ylim=c(-7,-1), xlab="ln(odds) not-vaccinated group", ylab="ln(odds) vaccinated group")
   par(opar)

})

############################################################################

### create dataset in long format
dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.colditz1994)
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)
dat.long$tpos <- dat.long$tneg <- dat.long$cpos <- dat.long$cneg <- NULL
levels(dat.long$group) <- c("EXP", "CON")

test_that("results for the bivariate model are correct.", {

   res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | trial, struct="UN", data=dat.long, method="ML")

   ### compare with results on pages 604-605 (in text)
   expect_equivalent(coef(res), c(-4.8337, -4.0960), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, c(1.4314, 2.4073), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, .9467, tolerance=.tol[["cor"]])

   res <- rma.mv(yi, vi, mods = ~ relevel(group, ref="CON"), random = ~ group | trial, struct="UN", data=dat.long, method="ML")

   ### compare with results on pages 604-605 (in text)
   expect_equivalent(coef(res), c(-4.0960, -0.7378), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.4347, 0.1797), tolerance=.tol[["se"]])

   ### estimated odds ratio
   tmp <- predict(res, newmods=1, intercept=FALSE, transf=exp, digits=3)
   expect_equivalent(tmp$pred,  0.4782, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.3362, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.6801, tolerance=.tol[["ci"]])

   ### amount of heterogeneity in log odds ratios
   tmp <- res$tau2[1] + res$tau2[2] - 2*res$rho*sqrt(res$tau2[1]*res$tau2[2])
   expect_equivalent(tmp, 0.3241, tolerance=.tol[["var"]])

})

############################################################################

test_that("results for the meta-regression analyses are correct.", {

   res <- rma(yi, vi, mods = ~ ablat, data=dat, method="ML")

   ### compare with results on pages 608-609 (in text)
   expect_equivalent(coef(res), c(0.3710, -0.0327), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.1061, 0.0034), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.0040, tolerance=.tol[["var"]])
   expect_equivalent(res$R2, 98.6691, tolerance=.tol[["r2"]])

   res <- rma.mv(yi, vi, mods = ~ group + group:I(ablat-33) - 1, random = ~ group | trial, struct="UN", data=dat.long, method="ML")

   ### compare with results on pages 612-613 (in text)
   expect_equivalent(coef(res), c(-4.8257, -4.1174, 0.0391, 0.0725), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.3129, 0.3061, 0.0224, 0.0219), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, c(1.2262, 1.1819), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 1.0000, tolerance=.tol[["cor"]])

   res <- rma.mv(yi, vi, mods = ~ relevel(group, ref="CON")*I(ablat-33), random = ~ group | trial, struct="UN", data=dat.long, method="ML")

   ### compare with results on pages 612-613 (in text)
   expect_equivalent(coef(res), c(-4.1174, -0.7083, 0.0725, -0.0333), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.3061, 0.0481, 0.0219, 0.0028), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, c(1.2262, 1.1819), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 1.0000, tolerance=.tol[["cor"]])

})
