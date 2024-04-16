### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see also: https://www.metafor-project.org/doku.php/analyses:vanhouwelingen2002

context("Checking analysis example: vanhouwelingen2002")

source("settings.r")

### load data
dat <- dat.colditz1994

### calculate log(OR)s and corresponding sampling variances
dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat)

### 'center' year variable
dat$year <- dat$year - 1900

test_that("results for the equal-effects model are correct.", {

   res <- rma(yi, vi, data=dat, method="EE")
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

   ### CI for tau^2 (profile likelihood method)
   tmp <- confint(res, type="PL")
   expect_equivalent(tmp$random[1,2], 0.1151, tolerance=.tol[["var"]])
   expect_equivalent(tmp$random[1,3], 0.8937, tolerance=.tol[["var"]])

   ### CI for tau^2 (Q-profile method)
   tmp <- confint(res)
   expect_equivalent(tmp$random[1,2], 0.1302, tolerance=.tol[["var"]]) ### 0.1350 based on a Satterthwaite approximation (page 597)
   expect_equivalent(tmp$random[1,3], 1.1812, tolerance=.tol[["var"]]) ### 1.1810 based on a Satterthwaite approximation (page 597)

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

   png(filename="images/test_analysis_example_vanhouwelingen2002_profile_test.png", res=200, width=1800, height=1600, type="cairo")
   profile(res, xlim=c(0.01,2), steps=200, log="x", cex=0, lwd=2, cline=TRUE, progbar=FALSE)
   abline(v=c(0.1151, 0.8937), lty="dotted")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_vanhouwelingen2002_profile_test.png", "images/test_analysis_example_vanhouwelingen2002_profile.png"))

})

test_that("forest plot of observed log(OR)s and corresponding BLUPs can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   res <- rma(yi, vi, data=dat, method="ML")
   sav <- blup(res)

   png(filename="images/test_analysis_example_vanhouwelingen2002_forest_light_test.png", res=200, width=1800, height=1400, family="mono")
   par(mar=c(5,5,1,2))
   forest(res, refline=res$b, addcred=TRUE, xlim=c(-7,7), alim=c(-3,3), slab=1:13, psize=0.8,
          ilab=paste0("(n = ", formatC(apply(dat[,c(4:7)], 1, sum), width=7, big.mark=","), ")"),
          ilab.xpos=-3.5, ilab.pos=2, rows=13:1+0.15, header="Trial (total n)", lty="dashed")
   arrows(sav$pi.lb, 13:1 - 0.15, sav$pi.ub, 13:1 - 0.15, length=0.035, angle=90, code=3)
   points(sav$pred, 13:1 - 0.15, pch=15, cex=0.8)
   dev.off()

   expect_true(.vistest("images/test_analysis_example_vanhouwelingen2002_forest_light_test.png", "images/test_analysis_example_vanhouwelingen2002_forest_light.png"))

   png(filename="images/test_analysis_example_vanhouwelingen2002_forest_dark_test.png", res=200, width=1800, height=1400, family="mono")
   setmfopt(theme="dark")
   par(mar=c(5,5,1,2))
   forest(res, refline=res$b, addcred=TRUE, xlim=c(-7,7), alim=c(-3,3), slab=1:13, psize=0.8,
          ilab=paste0("(n = ", formatC(apply(dat[,c(4:7)], 1, sum), width=7, big.mark=","), ")"),
          ilab.xpos=-3.5, ilab.pos=2, rows=13:1+0.15, header="Trial (total n)", lty="dashed")
   arrows(sav$pi.lb, 13:1 - 0.15, sav$pi.ub, 13:1 - 0.15, length=0.035, angle=90, code=3)
   points(sav$pred, 13:1 - 0.15, pch=15, cex=0.8)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_vanhouwelingen2002_forest_dark_test.png", "images/test_analysis_example_vanhouwelingen2002_forest_dark.png"))

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

   res <- rma(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, method="EE")

   png(filename="images/test_analysis_example_vanhouwelingen2002_labbe_light_test.png", res=200, width=1800, height=1400, type="cairo")
   par(mar=c(5,5,1,2))
   labbe(res, xlim=c(-7,-1), ylim=c(-7,-1),
         xlab="ln(odds) not-vaccinated group", ylab="ln(odds) vaccinated group")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_vanhouwelingen2002_labbe_light_test.png", "images/test_analysis_example_vanhouwelingen2002_labbe_light.png"))

   png(filename="images/test_analysis_example_vanhouwelingen2002_labbe_dark_test.png", res=200, width=1800, height=1400, type="cairo")
   setmfopt(theme="dark")
   par(mar=c(5,5,1,2))
   labbe(res, xlim=c(-7,-1), ylim=c(-7,-1),
         xlab="ln(odds) not-vaccinated group", ylab="ln(odds) vaccinated group")
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_vanhouwelingen2002_labbe_dark_test.png", "images/test_analysis_example_vanhouwelingen2002_labbe_dark.png"))

})

############################################################################

### create dataset in long format
dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.colditz1994)
dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)
dat.long$tpos <- dat.long$tneg <- dat.long$cpos <- dat.long$cneg <- NULL
levels(dat.long$group) <- c("CON", "EXP")

test_that("results for the bivariate model are correct.", {

   res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | trial, struct="UN", data=dat.long, method="ML", sparse=.sparse)

   ### compare with results on pages 604-605 (in text)
   expect_equivalent(coef(res), c(-4.0960, -4.8337), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, c(2.4073, 1.4314), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, .9467, tolerance=.tol[["cor"]])

   res <- rma.mv(yi, vi, mods = ~ group, random = ~ group | trial, struct="UN", data=dat.long, method="ML", sparse=.sparse)

   ### compare with results on pages 604-605 (in text)
   expect_equivalent(coef(res), c(-4.0960, -0.7378), tolerance=.tol[["coef"]])
   expect_equivalent(se(res), c(0.4347, 0.1797), tolerance=.tol[["se"]])

   ### estimated odds ratio
   tmp <- predict(res, newmods=1, intercept=FALSE, transf=exp, digits=3)
   expect_equivalent(tmp$pred,  0.4782, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.3362, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 0.6801, tolerance=.tol[["ci"]])

   ### amount of heterogeneity in log odds ratios
   tmp <- res$tau2[1] + res$tau2[2] - 2*res$rho*sqrt(res$tau2[1]*res$tau2[2])
   expect_equivalent(tmp, 0.3241, tolerance=.tol[["var"]])

   ### regression of log(odds)_EXP on log(odds)_CON
   res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | trial, struct="UN", data=dat.long, method="ML", sparse=.sparse)
   reg <- matreg(y=2, x=1, R=res$G, cov=TRUE, means=coef(res), n=res$g.levels.comb.k)
   expect_equivalent(reg$tab$beta, c(-1.8437, 0.7300), tolerance=.tol[["coef"]])
   expect_equivalent(reg$tab$se,   c( 0.3265, 0.0749), tolerance=.tol[["se"]])

   ### same idea but now use var-cov matrix of tau^2_1, tau_12, tau^2_2 for this
   res <- rma.mv(yi, vi, mods = ~ group - 1, random = ~ group | trial, struct="UN", data=dat.long, method="ML", cvvc="varcov", control=list(nearpd=TRUE), sparse=.sparse)
   reg <- matreg(y=2, x=1, R=res$G, cov=TRUE, means=coef(res), V=res$vvc)
   expect_equivalent(reg$tab$beta, c(-1.8437, 0.7300), tolerance=.tol[["coef"]])
   expect_equivalent(reg$tab$se,   c( 0.3548, 0.0866), tolerance=.tol[["se"]])

})

############################################################################

test_that("results for the meta-regression analyses are correct.", {

   res <- rma(yi, vi, mods = ~ ablat, data=dat, method="ML")

   ### compare with results on pages 608-609 (in text)
   expect_equivalent(coef(res), c(0.3710, -0.0327), tolerance=.tol[["coef"]])
   expect_equivalent(se(res), c(0.1061, 0.0034), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, 0.0040, tolerance=.tol[["var"]])
   expect_equivalent(res$R2, 98.6691, tolerance=.tol[["r2"]])

   res <- rma.mv(yi, vi, mods = ~ group + group:I(ablat-33) - 1, random = ~ group | trial, struct="UN", data=dat.long, method="ML", sparse=.sparse)

   ### compare with results on pages 612-613 (in text)
   expect_equivalent(coef(res), c(-4.1174, -4.8257, 0.0725, 0.0391), tolerance=.tol[["coef"]])
   expect_equivalent(se(res), c(0.3061, 0.3129, 0.0219, 0.0224), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, c(1.1819, 1.2262), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 1.0000, tolerance=.tol[["cor"]])

   res <- rma.mv(yi, vi, mods = ~ group*I(ablat-33), random = ~ group | trial, struct="UN", data=dat.long, method="ML", sparse=.sparse)

   ### compare with results on pages 612-613 (in text)
   expect_equivalent(coef(res), c(-4.1174, -0.7083, 0.0725, -0.0333), tolerance=.tol[["coef"]])
   expect_equivalent(se(res), c(0.3061, 0.0481, 0.0219, 0.0028), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, c(1.1819, 1.2262), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 1.0000, tolerance=.tol[["cor"]])

})

rm(list=ls())
