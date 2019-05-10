### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:raudenbush1985

context("Checking analysis example: raudenbush1985")

source("tolerances.r") # read in tolerances

### load data
dat <- dat.raudenbush1985

test_that("results are correct for the random-effects model.", {

   ### random-effects model
   res <- rma(yi, vi, data=dat, digits=3)

   ### compare with results on pages 83, 85, and 86 (in text)
   expect_equivalent(res$tau2,   0.0188, tolerance=.tol[["var"]])
   expect_equivalent(coef(res),  0.0837, tolerance=.tol[["coef"]])
   expect_equivalent(res$QE,    35.8295, tolerance=.tol[["test"]]) ### 35.85 in paper
   expect_equivalent(res$zval,   1.6208, tolerance=.tol[["test"]])

   ### empirical Bayes estimates
   tmp <- blup(res)

   out <- capture.output(print(tmp)) ### so that print.list.rma() is run (at least once)

   ### compare with results in Figure 2
   expect_equivalent(tmp$pred,  c(0.0543, 0.1006, -0.0064, 0.2144, 0.1051, -0.0082, 0.0174, -0.0293, 0.1604, 0.2485, 0.1618, 0.1102, 0.0646, 0.1105, -0.0288, 0.0258, 0.1905, 0.0744, 0.0248), tolerance=.tol[["pred"]])
   expect_equivalent(tmp$pi.lb, c(-0.1324, -0.1033, -0.2228, -0.0533, -0.1622, -0.1737, -0.1481, -0.2689, -0.0543, 0, -0.097, -0.1303, -0.192, -0.1463, -0.2405, -0.1906, -0.0076, -0.0808, -0.1954), tolerance=.tol[["ci"]])
   expect_equivalent(tmp$pi.ub, c(0.2411, 0.3045, 0.21, 0.4821, 0.3724, 0.1572, 0.1828, 0.2102, 0.3751, 0.497, 0.4206, 0.3507, 0.3212, 0.3672, 0.1829, 0.2422, 0.3886, 0.2295, 0.245), tolerance=.tol[["ci"]])

   ### empirical Bayes estimates (just the random effects)
   tmp <- ranef(res)

   expect_equivalent(tmp$pred, c(-0.0294, 0.0169, -0.0901, 0.1307, 0.0214, -0.0919, -0.0664, -0.1131, 0.0767, 0.1648, 0.0781, 0.0265, -0.0191, 0.0268, -0.1125, -0.0579, 0.1068, -0.0093, -0.0589), tolerance=.tol[["pred"]])
   expect_equivalent(tmp$pi.lb, c(-0.2187, -0.1852, -0.3019, -0.122, -0.231, -0.2659, -0.2403, -0.343, -0.1337, -0.0723, -0.1674, -0.2043, -0.2627, -0.217, -0.3207, -0.2697, -0.091, -0.1761, -0.2736), tolerance=.tol[["ci"]])
   expect_equivalent(tmp$pi.ub, c(0.1599, 0.219, 0.1216, 0.3834, 0.2738, 0.082, 0.1076, 0.1169, 0.2871, 0.4019, 0.3235, 0.2572, 0.2246, 0.2706, 0.0956, 0.1539, 0.3046, 0.1574, 0.1558), tolerance=.tol[["ci"]])

   skip_on_cran()

   ### profile tau^2
   opar <- par(no.readonly=TRUE)
   profile(res, xlim=c(0,.20), progbar=FALSE)
   par(opar)

   ### profile tau^2 (without 'xlim' specified)
   opar <- par(no.readonly=TRUE)
   profile(res, progbar=FALSE)
   par(opar)

   ### profile tau^2 (with parallel processing)
   opar <- par(no.readonly=TRUE)
   profile(res, xlim=c(0,.20), progbar=FALSE, parallel="snow")
   par(opar)

})

test_that("results are correct for the mixed-effects model.", {

   ### recode weeks variable
   dat$weeks.c <- ifelse(dat$weeks > 3, 3, dat$weeks)

   ### mixed-effects model
   res <- rma(yi, vi, mods = ~ weeks.c, data=dat, digits=3)

   ### compare with results on pages 90 and 92 (in text)
   expect_equivalent(res$tau2, 0.0000, tolerance=.tol[["var"]])
   expect_equivalent(coef(res), c(0.4072, -0.1572), tolerance=.tol[["coef"]])
   expect_equivalent(res$QE, 16.5708, tolerance=.tol[["test"]]) ### 16.58 in paper
   expect_equivalent(res$zval, c(4.6782, -4.3884), tolerance=.tol[["test"]])

   ### empirical Bayes estimates
   tmp <- blup(res)

   ### (results for this not given in chapter)
   expect_equivalent(tmp$pred,  c(0.0927, -0.0645, -0.0646, 0.4072, 0.4072, -0.0645, -0.0645, -0.0646, 0.4072, 0.2499, 0.4072, 0.4072, 0.2499, 0.0927, -0.0646, -0.0645, 0.2499, 0.0927, -0.0645), tolerance=.tol[["pred"]])
   expect_equivalent(tmp$pi.lb, c(0.0198, -0.1552, -0.1552, 0.2366, 0.2366, -0.1552, -0.1552, -0.1552, 0.2366, 0.1391, 0.2366, 0.2366, 0.1391, 0.0198, -0.1552, -0.1552, 0.1391, 0.0198, -0.1552), tolerance=.tol[["ci"]])
   expect_equivalent(tmp$pi.ub, c(0.1656, 0.0261, 0.0261, 0.5778, 0.5778, 0.0261, 0.0261, 0.0261, 0.5778, 0.3608, 0.5778, 0.5778, 0.3608, 0.1656, 0.0261, 0.0261, 0.3608, 0.1656, 0.0261), tolerance=.tol[["ci"]])

   ### empirical Bayes estimates (just the random effects)
   tmp <- ranef(res)

   expect_equivalent(tmp$pred,  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), tolerance=.tol[["pred"]])
   expect_equivalent(tmp$pi.lb, c(-0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016, -0.0016), tolerance=.tol[["ci"]])
   expect_equivalent(tmp$pi.ub, c(0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016, 0.0016), tolerance=.tol[["ci"]])

   ### predicted/fitted values
   tmp <- predict(res)

   ### (results for this not given in chapter)
   expect_equivalent(tmp$pred,  c(0.0927, -0.0645, -0.0645, 0.4072, 0.4072, -0.0645, -0.0645, -0.0645, 0.4072, 0.2499, 0.4072, 0.4072, 0.2499, 0.0927, -0.0645, -0.0645, 0.2499, 0.0927, -0.0645), tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, c(0.0198, -0.1552, -0.1552, 0.2366, 0.2366, -0.1552, -0.1552, -0.1552, 0.2366, 0.1391, 0.2366, 0.2366, 0.1391, 0.0198, -0.1552, -0.1552, 0.1391, 0.0198, -0.1552), tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, c(0.1656, 0.0261, 0.0261, 0.5778, 0.5778, 0.0261, 0.0261, 0.0261, 0.5778, 0.3607, 0.5778, 0.5778, 0.3607, 0.1656, 0.0261, 0.0261, 0.3607, 0.1656, 0.0261), tolerance=.tol[["ci"]])

   skip_on_cran()

   ### profile tau^2
   opar <- par(no.readonly=TRUE)
   profile(res, xlim=c(0,.06), progbar=FALSE)
   par(opar)

})
