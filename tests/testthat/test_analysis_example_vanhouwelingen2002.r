### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:vanhouwelingen2002

context("Checking analysis example vanhouwelingen2002")

### load data
dat <- get(data(dat.colditz1994, package="metafor"))

### calculate log(OR)s and corresponding sampling variances
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### 'center' year variable
dat$year <- dat$year - 1900

test_that("results for the fixed-effects model are correct.", {

   res <- rma(yi, vi, data=dat, method="FE")
   tmp <- predict(res, transf=exp, digits=3)

   ### compare with results on page 596 (in text)
   expect_equivalent(round(tmp$pred,3), .647)
   expect_equivalent(round(tmp$ci.lb,3), .595)
   expect_equivalent(round(tmp$ci.ub,3), .702) ### .703 in paper

})

test_that("results for the random-effects model are correct.", {

   res <- rma(yi, vi, data=dat, method="ML")
   tmp <- predict(res, transf=exp, digits=3)

   ### compare with results on page 597 (in text)
   expect_equivalent(round(tmp$pred,3), .476)
   expect_equivalent(round(tmp$ci.lb,3), .336)
   expect_equivalent(round(tmp$ci.ub,3), .675)
   expect_equivalent(round(res$tau2,3), .302)

   ### CI for tau^2
   tmp <- confint(res)
   expect_equivalent(round(tmp$random[1,2],3), 0.130) ### 0.135 based on a Satterthwaite approximation (page 597)
   expect_equivalent(round(tmp$random[1,3],3), 1.181) ### 1.181 based on a Satterthwaite approximation (page 597)

   ### CI for mu with Knapp & Hartung method
   res <- rma(yi, vi, data=dat, method="ML", knha=TRUE)
   tmp <- predict(res, transf=exp, digits=3)

   ### (results for this not given in paper)
   expect_equivalent(round(tmp$ci.lb,3), .318)
   expect_equivalent(round(tmp$ci.ub,3), .714)

})

test_that("profile plot for tau^2 can be drawn.", {

   skip_on_cran()

   res <- rma(yi, vi, data=dat, method="ML")

   opar <- par(no.readonly=TRUE)
   profile(res, xlim=c(.01,2), steps=200, log="x", cex=0, lwd=2, progbar=FALSE)
   abline(h=logLik(res) - 1.92, lwd=2)
   abline(v=c(.12, .89), lty="dashed")
   par(opar)

})

test_that("forest plot of observed log(OR)s and corresponding BLUPs can be drawn.", {

   skip_on_cran()

   res <- rma(yi, vi, data=dat, method="ML")
   sav <- blup(res)

   opar <- par(no.readonly=TRUE)
   par(family="mono", mar=c(5,4,1,2))
   forest(res, refline=res$b, addcred=TRUE, xlim=c(-7,8), alim=c(-3,3), slab=1:13, psize=.8,
          ilab=paste0("(n = ", formatC(apply(dat[,c(4:7)], 1, sum), width=7, big.mark=","), ")"),
          ilab.xpos=-3.5, ilab.pos=2, rows=13:1+.15)
   arrows(sav$pi.lb, 13:1 - .15, sav$pi.ub, 13:1 -.15, length=.03, angle=90, code=3, lty="dotted")
   points(sav$pred, 13:1 - .15, pch=15, cex=.8)
   text(-7, 15, "Trial (total n)", pos=4)
   text( 8, 15, "Log Odds Ratio [95% CI]", pos=2)
   par(opar)

})

test_that("the credibility/prediction interval is correct.", {

   res <- rma(yi, vi, data=dat, method="ML")

   ### computation as done in the paper
   tmp <- round(res$b + c(-1,+1) * qnorm(.975) * sqrt(res$tau2), 3)

   ### compare with results on page 599 (in text)
   expect_equivalent(tmp, c(-1.820, 0.336))

   ### computation done with metafor
   tmp <- predict(res, digits=3)

   ### (results for this not given in paper)
   expect_equivalent(round(tmp$cr.lb,3), -1.875)
   expect_equivalent(round(tmp$cr.ub,3),  0.391)

})

test_that("L'Abbe plot can be drawn.", {

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
   expect_equivalent(round(coef(res),3), c(-4.834, -4.096))
   expect_equivalent(round(res$tau2,3), c(1.431, 2.407))
   expect_equivalent(round(res$rho,3), .947)

   res <- rma.mv(yi, vi, mods = ~ relevel(group, ref="CON"), random = ~ group | trial, struct="UN", data=dat.long, method="ML")

   ### compare with results on pages 604-605 (in text)
   expect_equivalent(round(coef(res),3), c(-4.096, -0.738))
   expect_equivalent(round(res$se,3), c(0.435, 0.180))

   ### estimated odds ratio
   tmp <- predict(res, newmods=1, intercept=FALSE, transf=exp, digits=3)
   expect_equivalent(round(tmp$pred,3), 0.478)
   expect_equivalent(round(tmp$ci.lb,3), 0.336)
   expect_equivalent(round(tmp$ci.ub,3), 0.680)

   ### amount of heterogeneity in log odds ratios
   tmp <- res$tau2[1] + res$tau2[2] - 2*res$rho*sqrt(res$tau2[1]*res$tau2[2])
   expect_equivalent(round(tmp,3), 0.324)

})

############################################################################

test_that("results for the meta-regression analyses are correct.", {

   res <- rma(yi, vi, mods = ~ ablat, data=dat, method="ML")

   ### compare with results on pages 608-609 (in text)
   expect_equivalent(round(coef(res),3), c(0.371, -0.033))
   expect_equivalent(round(res$se,3), c(0.106, 0.003))
   expect_equivalent(round(res$tau2,3), 0.004)
   expect_equivalent(round(res$R2,1), 98.7)

   res <- rma.mv(yi, vi, mods = ~ group + group:I(ablat-33) - 1, random = ~ group | trial, struct="UN", data=dat.long, method="ML")

   ### compare with results on pages 612-613 (in text)
   expect_equivalent(round(coef(res),3), c(-4.826, -4.117, 0.039, 0.072))
   expect_equivalent(round(res$se,3), c(0.313, 0.306, 0.022, 0.022))
   expect_equivalent(round(res$tau2,2), c(1.23, 1.18)) ### rounded a bit more heavily, so 32-bit and 64-bit versions give same result
   expect_equivalent(round(res$rho,3), 1)

   res <- rma.mv(yi, vi, mods = ~ relevel(group, ref="CON")*I(ablat-33), random = ~ group | trial, struct="UN", data=dat.long, method="ML")

   ### compare with results on pages 612-613 (in text)
   expect_equivalent(round(coef(res),3), c(-4.117, -0.708, 0.072, -0.033))
   expect_equivalent(round(res$se,3), c(0.306, 0.048, 0.022, 0.003))
   expect_equivalent(round(res$tau2,2), c(1.23, 1.18)) ### rounded a bit more heavily, so 32-bit and 64-bit versions give same result
   expect_equivalent(round(res$rho,3), 1)

})
