### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:berkey1998

source("tolerances.r") # read in tolerances

context("Checking analysis example: berkey1998")

### load data
dat <- dat.berkey1998

### construct variance-covariance matrix of the observed outcomes
V <- bldiag(lapply(split(dat[,c("v1i", "v2i")], dat$trial), as.matrix))

test_that("results are correct for the multiple outcomes random-effects model.", {

   ### multiple outcomes random-effects model (with ML estimation)
   res <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")
   out <- capture.output(print(res)) ### so that print.rma.mv() is run (at least once)

   ### (results for this model not given in paper)
   expect_equivalent(coef(res), c(-0.3379, 0.3448), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.0798, 0.0495), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, c(0.0261, 0.0070), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 0.6992, tolerance=.tol[["cor"]])

})

test_that("results are correct for the multiple outcomes mixed-effects (meta-regression) model.", {

   ### multiple outcomes mixed-effects (meta-regression) model (with ML estimation)
   res <- rma.mv(yi, V, mods = ~ outcome + outcome:I(year - 1983) - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")

   ### compare with results on page 2545 (Table II)
   expect_equivalent(coef(res), c(-0.3351, 0.3479, -0.0108, 0.0010), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.0787, 0.0520, 0.0243, 0.0154), tolerance=.tol[["se"]])
   expect_equivalent(res$tau2, c(0.0250, 0.0080), tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 0.6587, tolerance=.tol[["cor"]])

   ### compute the covariance
   tmp <- res$rho*sqrt(res$tau2[1]*res$tau2[2])
   expect_equivalent(tmp, 0.0093, tolerance=.tol[["cov"]])

   ### test the difference in slopes
   res <- rma.mv(yi, V, mods = ~ outcome*I(year - 1983) - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")

   ### (results for this model not given in paper)
   expect_equivalent(coef(res), c(-0.3351, 0.3479, -0.0108, 0.0118), tolerance=.tol[["coef"]])
   expect_equivalent(res$se, c(0.0787, 0.0520, 0.0243, 0.0199), tolerance=.tol[["se"]])
   expect_equivalent(res$pval, c(0.0000, 0.0000, 0.6563, 0.5534), tolerance=.tol[["pval"]])

})

test_that("results are correct when testing var-cov structures against each other with LRTs.", {

   ### test whether the amount of heterogeneity is the same in the two outcomes
   res1 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")
   res0 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="CS", data=dat, method="ML")
   tmp <- anova(res0, res1)
   out <- capture.output(print(tmp)) ### so that print.anova.rma() is run (at least once)

   ### (results for this not given in paper)
   expect_equivalent(tmp$pval, 0.2597, tolerance=.tol[["pval"]])

   ### test the correlation among the true effects
   res1 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML")
   res0 <- rma.mv(yi, V, mods = ~ outcome - 1, random = ~ outcome | trial, struct="UN", data=dat, method="ML", rho=0)
   tmp <- anova(res0, res1)

   ### (results for this not given in paper)
   expect_equivalent(tmp$pval, 0.2452, tolerance=.tol[["pval"]])

})
