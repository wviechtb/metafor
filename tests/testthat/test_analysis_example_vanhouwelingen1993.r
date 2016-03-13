### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:vanhouwelingen1993

context("Checking analysis example: vanhouwelingen1993")

### load data
dat <- get(data(dat.collins1985a, package="metafor"))

test_that("the log likelihood plot can be created.", {

   skip_on_cran()

   opar <- par(no.readonly=TRUE)
   llplot(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat,
          xlim=c(-4,4), lwd=1, col="black", refline=NA, drop00=FALSE)
   par(opar)

})

test_that("results of the fixed-effects conditional logistic model are correct.", {

   skip_on_cran()

   res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="FE")

   ### compare with results on page 2275 (in text)
   expect_equivalent(round(coef(res),3), 0.122)
   expect_equivalent(round(res$se,3), 0.100)
   expect_equivalent(round(res$ci.lb,3), -0.074)
   expect_equivalent(round(res$ci.ub,3), 0.317) ### 0.31 in paper
   expect_equivalent(round(c(logLik(res)), 3), -53.679)

})

test_that("results of the random-effects conditional logistic model are correct.", {

   skip_on_cran()

   res <- rma.glmm(measure="OR", ai=b.xci, n1i=nci, ci=b.xti, n2i=nti, data=dat, model="CM.EL", method="ML")

   ### compare with results on page 2277 (in text)
   expect_equivalent(round(coef(res),3), 0.175)
   expect_equivalent(round(res$se,3), 0.134)
   expect_equivalent(round(res$ci.lb,3), -0.088)
   expect_equivalent(round(res$ci.ub,3), 0.437)
   expect_equivalent(round(c(logLik(res)), 3), -52.989)
   expect_equivalent(round(res$tau2,3), 0.119)

})
