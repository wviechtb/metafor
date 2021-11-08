### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: selmodel() function")

source("tolerances.r") # read in tolerances

test_that("results are correct for a step function model.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   dat <- dat.hackshaw1998
   res <- rma(yi, vi, data=dat)

   sav <- selmodel(res, type="stepfun", steps=c(0.05, 0.10, 0.50, 1.00))
   out <- capture.output(print(sav))

   expect_equivalent(coef(sav)$delta, c(1, 2.422079, 0.977543, 0.396713), tolerance=.tol[["coef"]])
   expect_equivalent(sav$se.delta, c(NA, 1.66085, 0.820387, 0.469235), tolerance=.tol[["se"]])
   expect_equivalent(sav$LRT, 7.066137, tolerance=.tol[["test"]])
   expect_identical(sav$LRTdf, 3L)
   expect_equivalent(sav$tau2, 0.03071325, tolerance=.tol[["var"]])

   opar <- par(no.readonly=TRUE)
   plot(sav)
   par(opar)

   tmp <- confint(sav)
   expect_equivalent(tmp[[1]]$random[1,], c(0.030713, 0.000224, 0.135284), tolerance=.tol[["var"]])
   expect_equivalent(tmp[[2]]$random[1,], c(2.422079, 0.665133, 9.915798), tolerance=.tol[["coef"]])
   expect_equivalent(tmp[[3]]$random[1,], c(0.977543, 0.209558, 5.386044), tolerance=.tol[["coef"]])
   expect_equivalent(tmp[[4]]$random[1,], c(0.396713, 0.040198, 4.119681), tolerance=.tol[["coef"]])

})

test_that("results are correct for the beta function model.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   # data from Baskerville, N. B., Liddy, C., & Hogg, W. (2012). Systematic review and meta-analysis of practice facilitation within primary care settings. Annals of Family Medicine, 10(1), 63-74.

   yi  <- c(1.01, 0.82, 0.59, 0.44, 0.84, 0.73, 1.12, 0.04, 0.24, 0.32, 1.04, 1.31, 0.59, 0.66, 0.62, 0.47, 1.08, 0.98, 0.26, 0.39, 0.60, 0.94, 0.11)
   sei <- c(0.52, 0.46, 0.23, 0.18, 0.29, 0.29, 0.36, 0.37, 0.15, 0.40, 0.32, 0.57, 0.29, 0.19, 0.31, 0.27, 0.32, 0.32, 0.18, 0.18, 0.31, 0.53, 0.27)
   xi  <- c(1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1)

   res <- rma(yi, sei^2, method="ML")

   sav <- selmodel(res, type="beta", delta=c(1,1))
   expect_equivalent(logLik(res), logLik(sav), tolerance=.tol[["fit"]])

   sav <- selmodel(res, type="beta")
   out <- capture.output(print(sav))

   expect_equivalent(coef(sav)$delta, c(0.4731131, 4.4613093), tolerance=.tol[["coef"]])
   expect_equivalent(sav$se.delta, c(0.2352523, 2.1845971), tolerance=.tol[["se"]])
   expect_equivalent(sav$LRT, 7.846906, tolerance=.tol[["test"]])
   expect_identical(sav$LRTdf, 2L)
   expect_equivalent(sav$tau2, 0.00000243, tolerance=.tol[["var"]])

   opar <- par(no.readonly=TRUE)
   plot(sav)
   par(opar)

   res <- rma(yi, sei^2, mods = ~ xi, method="ML")

   sav <- selmodel(res, type="beta")
   out <- capture.output(print(sav))

   expect_equivalent(coef(sav)$delta, c(0.4200973, 5.0959707), tolerance=.tol[["coef"]])
   expect_equivalent(sav$se.delta, c(0.239128, 2.410997), tolerance=.tol[["se"]])
   expect_equivalent(sav$LRT, 9.044252, tolerance=.tol[["test"]])
   expect_identical(sav$LRTdf, 2L)
   expect_equivalent(sav$tau2, 0.00000193, tolerance=.tol[["var"]])
   expect_equivalent(coef(sav)$beta, c(0.1343001, -0.1363559), tolerance=.tol[["coef"]])
   expect_equivalent(sav$se, c(0.1707418, 0.1244394), tolerance=.tol[["se"]])

})

test_that("results are correct for the various exponential function models.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   # data from Preston, C., Ashby, D., & Smyth, R. (2004). Adjusting for publication bias: Modelling the selection process. Journal of Evaluation in Clinical Practice, 10(2), 313-322.   ai  <- c(4,0,34,7,6,1,0,11,2,0,0,33)

   ai  <- c(4,0,34,7,6,1,0,11,2,0,0,33)
   n1i <- c(19,18,341,71,45,94,22,88,82,33,15,221)
   ci  <- c(5,0,50,16,5,8,0,12,7,0,1,43)
   n2i <- c(19,18,334,69,44,96,22,82,84,30,20,218)
   dat <- escalc(measure="OR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, drop00=TRUE)

   expect_warning(res <- rma(yi, vi, data=dat, method="EE"))

   alternative <- "less"

   sav1 <- selmodel(res, type="halfnorm", alternative=alternative)
   sav2 <- selmodel(res, type="negexp",   alternative=alternative)
   sav3 <- selmodel(res, type="logistic", alternative=alternative)
   sav4 <- selmodel(res, type="power",    alternative=alternative)

   expect_equivalent(c(sav1$delta, sav2$delta, sav3$delta, sav4$delta), c(3.162948, 2.656714, 3.339338, 1.458923), tolerance=.tol[["coef"]])
   expect_equivalent(c(sav1$se.delta, sav2$se.delta, sav3$se.delta, sav4$se.delta), c(2.988922, 2.347468, 2.388776, 1.393725), tolerance=.tol[["se"]])

   opar <- par(no.readonly=TRUE)
   tmp <- profile(sav1, progbar=FALSE)
   par(opar)

   expect_equivalent(tmp$ll, c(NA, -6.569986, -6.35659, -6.210436, -6.121035, -6.07939, -6.077928, -6.110356, -6.171488, -6.257068, -6.363607, -6.488238, -6.628599, -6.782733, -6.949015, -7.126075, -7.312763, -7.508097, -7.711241, -7.921472), tolerance=.tol[["fit"]])

   sav1 <- selmodel(res, type="halfnorm", prec="sei", alternative=alternative, scaleprec=FALSE)
   sav2 <- selmodel(res, type="negexp",   prec="sei", alternative=alternative, scaleprec=FALSE)
   sav3 <- selmodel(res, type="logistic", prec="sei", alternative=alternative, scaleprec=FALSE)
   sav4 <- selmodel(res, type="power",    prec="sei", alternative=alternative, scaleprec=FALSE)

   expect_equivalent(c(sav1$delta, sav2$delta, sav3$delta, sav4$delta), c(3.506329, 2.279336, 3.017851, 1.444174), tolerance=.tol[["coef"]])
   expect_equivalent(c(sav1$se.delta, sav2$se.delta, sav3$se.delta, sav4$se.delta), c(3.387300, 2.133013, 2.315789, 1.381633), tolerance=.tol[["se"]])

   sav1 <- selmodel(res, type="halfnorm", prec="sei", alternative=alternative, steps=.05)
   sav2 <- selmodel(res, type="negexp",   prec="sei", alternative=alternative, steps=.05)
   sav3 <- selmodel(res, type="logistic", prec="sei", alternative=alternative, steps=.05)
   sav4 <- selmodel(res, type="power",    prec="sei", alternative=alternative, steps=.05, control=list(hessianCtrl=list(r=8)))

   expect_equivalent(c(sav1$delta, sav2$delta, sav3$delta, sav4$delta), c(5.832106, 3.819847, 5.041039, 2.399645), tolerance=.tol[["coef"]])
   expect_equivalent(c(sav1$se.delta, sav2$se.delta, sav3$se.delta, sav4$se.delta), c(5.644466, 3.627467, 2.306998, 2.134629), tolerance=.tol[["se"]])

   sav <- selmodel(res, type="negexppow", alternative=alternative)
   expect_equivalent(sav$delta, c(2.673818, 1.153199), tolerance=.tol[["coef"]])
   expect_equivalent(sav$se.delta, c(2.363403, 2.143849), tolerance=.tol[["se"]])

})

test_that("results are correct for a pirori chosen step function models.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   tab <- data.frame(
      steps = c(0.005, 0.01, 0.05, 0.10, 0.25, 0.35, 0.50, 0.65, 0.75, 0.90, 0.95, 0.99, 0.995, 1),
      delta.mod.1 = c(1, 0.99, 0.95, 0.80, 0.75, 0.65, 0.60, 0.55, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50),
      delta.sev.1 = c(1, 0.99, 0.90, 0.75, 0.60, 0.50, 0.40, 0.35, 0.30, 0.25, 0.10, 0.10, 0.10, 0.10),
      delta.mod.2 = c(1, 0.99, 0.95, 0.90, 0.80, 0.75, 0.60, 0.60, 0.75, 0.80, 0.90, 0.95, 0.99, 1.00),
      delta.sev.2 = c(1, 0.99, 0.90, 0.75, 0.60, 0.50, 0.25, 0.25, 0.50, 0.60, 0.75, 0.90, 0.99, 1.00))

   dat <- dat.cohen1981
   dat <- escalc(measure="ZCOR", ri=ri, ni=ni, data=dat[c(1,4,5)])

   res <- rma(yi, vi, data=dat, method="ML")

   sav <- lapply(tab[-1], function(x) selmodel(res, type="stepfun", steps=tab$steps, delta=x, defmap=TRUE))
   expect_equivalent(sapply(sav, function(x) x$beta), c(0.351894, 0.321518, 0.362019, 0.33218), tolerance=.tol[["coef"]])
   expect_equivalent(sapply(sav, function(x) x$tau2), c(0.0045, 0.009544, 0.002774, 0.005652),  tolerance=.tol[["var"]])

})
