### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

### see: https://www.metafor-project.org/doku.php/analyses:vanhouwelingen1993

context("Checking analysis example: vanhouwelingen1993")

source("settings.r")

### load data
dat <- dat.collins1985a

test_that("the log likelihood plot can be created.", {

   skip_on_cran()

   png(filename="images/test_analysis_example_vanhouwelingen1993_llplot_light_test.png", res=200, width=1800, height=1200, type="cairo")
   par(mar=c(5,5,1,2))
   expect_warning(llplot(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat,
                  xlim=c(-4,4), lwd=1, col="black", refline=NA, drop00=FALSE))
   dev.off()

   expect_true(.vistest("images/test_analysis_example_vanhouwelingen1993_llplot_light_test.png", "images/test_analysis_example_vanhouwelingen1993_llplot_light.png"))

   png(filename="images/test_analysis_example_vanhouwelingen1993_llplot_dark_test.png", res=200, width=1800, height=1200, type="cairo")
   setmfopt(theme="dark")
   par(mar=c(5,5,1,2))
   expect_warning(llplot(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat,
                  xlim=c(-4,4), lwd=1, col="white", refline=NA, drop00=FALSE))
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_analysis_example_vanhouwelingen1993_llplot_dark_test.png", "images/test_analysis_example_vanhouwelingen1993_llplot_dark.png"))

})

test_that("results of the equal-effects conditional logistic model are correct.", {

   skip_on_cran()

   expect_warning(res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="EE"))

   ### compare with results on page 2275 (in text)
   expect_equivalent(coef(res), 0.1216, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.0995, tolerance=.tol[["se"]])
   expect_equivalent(res$ci.lb, -0.0734, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 0.3165, tolerance=.tol[["ci"]]) ### 0.31 in paper (rounded a bit more heavily, so 32-bit and 64-bit versions give same result)
   expect_equivalent(c(logLik(res)), -53.6789, tolerance=.tol[["fit"]])

   ### run with control(dnchgcalc="dnoncenhypergeom")
   expect_warning(res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="EE", control=list(dnchgcalc="dnoncenhypergeom")))

   ### some very minor discrepancies
   expect_equivalent(coef(res), 0.1216, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.0996, tolerance=.tol[["se"]])
   expect_equivalent(res$ci.lb, -0.0735, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 0.3167, tolerance=.tol[["ci"]])
   expect_equivalent(c(logLik(res)), -53.6789, tolerance=.tol[["fit"]])

})

test_that("results of the random-effects conditional logistic model are correct.", {

   skip_on_cran()

   expect_warning(res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="ML"))

   ### compare with results on page 2277 (in text)
   expect_equivalent(coef(res), 0.1744, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.1364, tolerance=.tol[["se"]])
   expect_equivalent(res$ci.lb, -0.0929, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 0.4417, tolerance=.tol[["ci"]])
   expect_equivalent(c(logLik(res)), -52.99009, tolerance=.tol[["fit"]])
   expect_equivalent(res$tau2, 0.1192, tolerance=.tol[["var"]])

   ### run with control(dnchgcalc="dnoncenhypergeom")
   expect_warning(res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="ML", control=list(dnchgcalc="dnoncenhypergeom")))

   ### no discrepancies
   expect_equivalent(coef(res), 0.1744, tolerance=.tol[["coef"]])
   expect_equivalent(se(res), 0.1364, tolerance=.tol[["se"]])
   expect_equivalent(res$ci.lb, -0.0930, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 0.4418, tolerance=.tol[["ci"]])
   expect_equivalent(c(logLik(res)), -52.99009, tolerance=.tol[["fit"]])
   expect_equivalent(res$tau2, 0.1192, tolerance=.tol[["var"]])

})

rm(list=ls())
