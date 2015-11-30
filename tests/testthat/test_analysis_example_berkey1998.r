### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:berkey1998

context("Checking analysis example berkey1998")

### load data
dat <- get(data(dat.berkey1998, package="metafor"))

### construct variance-covariance matrix of the observed outcomes
V <- bldiag(lapply(split(dat[,c("v1i", "v2i")], dat$trial), as.matrix))

test_that("results are correct for the multiple outcomes random-effects model.", {

   ### multiple outcomes random-effects model (with ML estimation)
   res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")

   ### (results for this model not given in paper)
   expect_equivalent(round(coef(res),3), c(-0.338, 0.345))
   expect_equivalent(round(res$se,3), c(0.080, 0.049))
   expect_equivalent(round(res$tau2,3), c(0.026, 0.007))
   expect_equivalent(round(res$rho,3), 0.699)

})

test_that("results are correct for the multiple outcomes mixed-effects (meta-regression) model.", {

   ### multiple outcomes mixed-effects (meta-regression) model (with ML estimation)
   res <- rma.mv(yi, V, mods = ~ outcome + outcome:I(year - 1983) - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")

   ### compare with results on page 2545 (Table II)
   expect_equivalent(round(coef(res),3), c(-0.335, 0.348, -0.011, 0.001))
   expect_equivalent(round(res$se,3), c(0.079, 0.052, 0.024, 0.015))
   expect_equivalent(round(res$tau2,3), c(0.025, 0.008))
   expect_equivalent(round(res$rho,3), 0.659)

   ### compute the covariance
   tmp <- round(res$rho*sqrt(res$tau2[1]*res$tau2[2]),3)
   expect_equivalent(tmp, 0.009)

   ### test the difference in slopes
   res <- rma.mv(yi, V, mods = ~ outcome*I(year - 1983) - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")

   ### (results for this model not given in paper)
   expect_equivalent(round(coef(res),3), c(-0.335, 0.348, -0.011, 0.012))
   expect_equivalent(round(res$se,3), c(0.079, 0.052, 0.024, 0.020))
   expect_equivalent(round(res$pval,3), c(0.000, 0.000, 0.656, 0.553))

})

test_that("results are correct when testing var-cov structures against each other with LRTs.", {

   ### test whether the amount of heterogeneity is the same in the two outcomes
   res1 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")
   res0 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="CS", data=dat, method="ML")
   tmp <- anova(res0, res1)

   ### (results for this not given in paper)
   expect_equivalent(round(tmp$pval,4), 0.2597)

   ### test the correlation among the true effects
   res1 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")
   res0 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML", rho=0)
   tmp <- anova(res0, res1)

   ### (results for this not given in paper)
   expect_equivalent(round(tmp$pval,4), 0.2452)

})
