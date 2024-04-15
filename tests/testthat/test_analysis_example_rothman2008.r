### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/analyses:rothman2008

context("Checking analysis example: rothman2008")

source("settings.r")

############################################################################

### create dataset (Table 15-1)
dat <- data.frame(
age = c("Age <55", "Age 55+"),
ai = c(8,22),
bi = c(98,76),
ci = c(5,16),
di = c(115,69),
stringsAsFactors=FALSE)

test_that("the to.table() function works.", {

   tmp <- to.table(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", slab=age, rows=c("Tolbutamide", "Placebo"), cols=c("Dead", "Surviving"))

   expected <- structure(c(8, 5, 98, 115, 22, 16, 76, 69), .Dim = c(2L, 2L, 2L), .Dimnames = list(c("Tolbutamide", "Placebo"), c("Dead", "Surviving"), c("Age <55", "Age 55+")))

   ### compare with data in Table 15-1
   expect_equivalent(tmp, expected)

})

test_that("the to.long() function works.", {

   tmp <- to.long(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", slab=age)

   expected <- structure(list(age = c("Age <55", "Age <55", "Age 55+", "Age 55+"),
      ai = c(8, 8, 22, 22), bi = c(98, 98, 76, 76), ci = c(5, 5, 16, 16), di = c(115, 115, 69, 69),
      study = structure(c(2L, 2L, 1L, 1L), .Label = c("Age 55+", "Age <55"), class = "factor"),
      group = structure(c(2L, 1L, 2L, 1L), .Label = c("2", "1"), class = "factor"),
      out1 = c(8, 5, 22, 16), out2 = c(98, 115, 76, 69)), class = "data.frame", row.names = c(NA, 4L))

   expect_equivalent(tmp, expected)

})

test_that("the stratum-specific and crude risk differences are computed correctly.", {

   ### stratum-specific risk differences
   tmp <- summary(escalc(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RD", digits=3, append=FALSE))
   tmp <- as.matrix(tmp[1:4])

   expected <- structure(c(0.0338, 0.0363, 0.001, 0.0036, 0.0315, 0.0598, 1.0738, 0.6064), .Dim = c(2L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-1
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

   ### crude risk difference
   tmp <- summary(escalc(ai=sum(ai), bi=sum(bi), ci=sum(ci), di=sum(di), data=dat, measure="RD", digits=3, append=FALSE))
   tmp <- as.matrix(tmp[1:4])

   expected <- structure(c(0.0446, 0.0011, 0.0326, 1.3683), .Dim = c(1L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-1
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

})

test_that("the stratum-specific and crude risk ratios are computed correctly.", {

   ### stratum-specific risk ratios
   tmp <- summary(escalc(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RR", digits=2), transf=exp, append=FALSE)
   tmp <- as.matrix(tmp)

   expected <- structure(c(1.8113, 1.1926, 0.6112, 0.6713, 5.3679, 2.1188), .Dim = 2:3, .Dimnames = list(NULL, c("yi", "ci.lb", "ci.ub")))

   ### compare with data in Table 15-1
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

   ### crude risk ratio
   tmp <- summary(escalc(ai=sum(ai), bi=sum(bi), ci=sum(ci), di=sum(di), data=dat, measure="RR", digits=2, append=FALSE), transf=exp)
   tmp <- as.matrix(tmp)

   expected <- structure(c(1.4356, 0.851, 2.4216), .Dim = c(1L, 3L), .Dimnames = list(NULL, c("yi", "ci.lb", "ci.ub")))

   ### compare with data in Table 15-1
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

})

test_that("results are correct for Mantel-Haenszel method.", {

   ### Mantel-Haenszel method with risk differences
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RD", digits=3, level=90)
   out <- capture.output(print(res)) ### so that print.rma.mh() is used

   expect_equivalent(coef(res),  0.0349, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.0176, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub,  0.0874, tolerance=.tol[["ci"]]) ### 0.088 in chapter
   expect_equivalent(res$QE,  0.0017, tolerance=.tol[["test"]]) ### 0.001 in chapter
   expect_equivalent(res$QEp, 0.9669, tolerance=.tol[["pval"]])

   ### Mantel-Haenszel method with risk ratios
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="RR", digits=2, level=90)
   out <- capture.output(print(res)) ### so that print.rma.mh() is used

   expect_equivalent(coef(res),  0.2818, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.1442, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub,  0.7078, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE,  0.4472, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0.5037, tolerance=.tol[["pval"]])

   tmp <- c(confint(res, transf=exp)$fixed)
   expect_equivalent(tmp, c(1.3256, 0.8658, 2.0296), tolerance=.tol[["ci"]])

   ### Mantel-Haenszel method with odds ratios
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", correct=FALSE, digits=2, level=90)
   out <- capture.output(print(res)) ### so that print.rma.mh() is used

   expect_equivalent(coef(res),  0.3387, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.1731, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub,  0.8505, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE,  0.3474, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0.5556, tolerance=.tol[["pval"]])
   expect_equivalent(res$CO,  1.1976, tolerance=.tol[["test"]])
   expect_equivalent(res$COp, 0.2738, tolerance=.tol[["pval"]])
   expect_equivalent(res$MH,  1.1914, tolerance=.tol[["test"]])
   expect_equivalent(res$MHp, 0.2750, tolerance=.tol[["pval"]])
   expect_equivalent(res$TA,  0.3489, tolerance=.tol[["test"]])
   expect_equivalent(res$TAp, 0.5547, tolerance=.tol[["pval"]])

   tmp <- c(confint(res, transf=exp)$fixed)
   expect_equivalent(tmp, c(1.4031, 0.8411, 2.3409), tolerance=.tol[["ci"]])

   skip_on_cran()

   ### conditional MLE of the odds ratio
   res <- rma.glmm(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", model="CM.EL", method="EE")

   expect_equivalent(coef(res),  0.3381, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, -0.2699, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub,  0.9461, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE.Wld,  0.3480, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.Wld, 0.5552, tolerance=.tol[["pval"]])
   expect_equivalent(res$QE.LRT,  0.3502, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.LRT, 0.5540, tolerance=.tol[["pval"]])

   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred,  1.4022, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 0.7634, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 2.5756, tolerance=.tol[["ci"]])

})

############################################################################

### create dataset (Table 15-2)
dat <- data.frame(
age = c("35-44", "45-54", "55-64", "65-74", "75-84"),
x1i = c(32, 104, 206, 186, 102),
t1i = c(52407, 43248, 28612, 12663, 5317) / 10000,
x2i = c(2, 12, 28, 28, 31),
t2i = c(18790, 10673, 5710, 2585, 1462) / 10000,
stringsAsFactors=FALSE)

test_that("the to.table() function works.", {

   tmp <- to.table(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", slab=age, rows=c("Smokers", "Nonsmokers"), cols=c("Deaths", "Years"))

   expected <- structure(c(32, 2, 5.2407, 1.879, 104, 12, 4.3248, 1.0673, 206, 28, 2.8612, 0.571, 186, 28, 1.2663, 0.2585, 102, 31, 0.5317, 0.1462), .Dim = c(2L, 2L, 5L), .Dimnames = list(c("Smokers", "Nonsmokers"), c("Deaths", "Years"), c("35-44", "45-54", "55-64", "65-74", "75-84")))

   ### compare with data in Table 15-2
   expect_equivalent(tmp, expected)

})

test_that("the to.long() function works.", {

   tmp <- to.long(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", slab=age)

   expected <- structure(list(age = c("35-44", "35-44", "45-54", "45-54", "55-64", "55-64", "65-74", "65-74", "75-84", "75-84"),
      x1i = c(32, 32, 104, 104, 206, 206, 186, 186, 102, 102),
      t1i = c(5.2407, 5.2407, 4.3248, 4.3248, 2.8612, 2.8612, 1.2663, 1.2663, 0.5317, 0.5317),
      x2i = c(2, 2, 12, 12, 28, 28, 28, 28, 31, 31),
      t2i = c(1.879, 1.879, 1.0673, 1.0673, 0.571, 0.571, 0.2585, 0.2585, 0.1462, 0.1462),
      study = structure(c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L), .Label = c("35-44", "45-54", "55-64", "65-74", "75-84"), class = "factor"),
      group = structure(c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L), .Label = c("2", "1"), class = "factor"),
      events = c(32, 2, 104, 12, 206, 28, 186, 28, 102, 31),
      ptime = c(5.2407, 1.879, 4.3248, 1.0673, 2.8612, 0.571, 1.2663, 0.2585, 0.5317, 0.1462)),
      class = "data.frame", row.names = c(NA, 10L))

   expect_equivalent(tmp, expected)

})

test_that("the stratum-specific and crude rate differences are computed correctly.", {

   ### stratum-specific rate differences
   tmp <- summary(escalc(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRD", digits=1, append=FALSE))
   tmp <- as.matrix(tmp[1:4])

   expected <- structure(c(5.0417, 12.804, 22.961, 38.5674, -20.2008, 1.7316, 16.0947, 111.0423, 535.0172, 1811.1307, 1.3159, 4.0118, 10.5377, 23.1304, 42.5574, 3.8313, 3.1916, 2.1789, 1.6674, -0.4747), .Dim = c(5L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-2
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

   ### crude rate difference
   tmp <- summary(escalc(x1i=sum(x1i), x2i=sum(x2i), t1i=sum(t1i), t2i=sum(t2i), data=dat, measure="IRD", digits=1, append=FALSE))
   tmp <- as.matrix(tmp[1:4])

   expected <- structure(c(18.537, 9.6796, 3.1112, 5.9581), .Dim = c(1L, 4L), .Dimnames = list(NULL, c("yi", "vi", "sei", "zi")))

   ### compare with data in Table 15-2
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

})

test_that("the stratum-specific and crude rate ratios are computed correctly.", {

   ### stratum-specific rate ratios
   tmp <- summary(escalc(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=1, append=FALSE), transf=exp)
   tmp <- as.matrix(tmp)

   expected <- structure(c(5.7366, 2.1388, 1.4682, 1.3561, 0.9047, 1.3748, 1.1767, 0.9894, 0.9115, 0.6053, 23.9371, 3.8876, 2.1789, 2.0176, 1.3524), .Dim = c(5L, 3L), .Dimnames = list(NULL, c("yi", "ci.lb", "ci.ub")))

   ### compare with data in Table 15-2
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

   ### crude rate ratio
   tmp <- summary(escalc(x1i=sum(x1i), x2i=sum(x2i), t1i=sum(t1i), t2i=sum(t2i), data=dat, measure="IRR", digits=1, append=FALSE), transf=exp)
   tmp <- as.matrix(tmp)

   expected <- structure(c(1.7198, 1.394, 2.1219), .Dim = c(1L, 3L), .Dimnames = list(NULL, c("yi", "ci.lb", "ci.ub")))

   ### compare with data in Table 15-2
   expect_equivalent(tmp, expected, tolerance=.tol[["misc"]])

})

test_that("results are correct for Mantel-Haenszel method.", {

   ### Mantel-Haenszel method with rate differences
   res <- rma.mh(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRD", digits=2, level=90)

   expect_equivalent(coef(res), 11.4392, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb,  6.3498, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 16.5286, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE, 26.8758, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0.0000, tolerance=.tol[["pval"]])

   ### Mantel-Haenszel method with rate ratios
   res <- rma.mh(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=2, level=90)

   expect_equivalent(coef(res), 0.3539, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.1776, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 0.5303, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE, 10.4117, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0.0340, tolerance=.tol[["pval"]])
   expect_equivalent(res$MH, 10.7021, tolerance=.tol[["test"]])
   expect_equivalent(res$MHp, 0.0011, tolerance=.tol[["pval"]])

   tmp <- c(confint(res, transf=exp)$fixed)
   expect_equivalent(tmp, c(1.4247, 1.1944, 1.6994), tolerance=.tol[["ci"]])

   ### Mantel-Haenszel test without continuity correction
   res <- rma.mh(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", level=90, correct=FALSE)

   expect_equivalent(res$MH, 11.0162, tolerance=.tol[["test"]])
   expect_equivalent(res$MHp, 0.0009, tolerance=.tol[["pval"]])

   skip_on_cran()

   ### unconditional MLE of the rate ratio
   res <- rma.glmm(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=2, level=90, model="UM.FS", method="EE")

   expect_equivalent(coef(res), 0.3545, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.1779, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 0.5312, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE.Wld, 10.1991, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.Wld, 0.0372, tolerance=.tol[["pval"]])
   expect_equivalent(res$QE.LRT, 12.1324, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.LRT, 0.0164, tolerance=.tol[["pval"]])

   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred,  1.4255, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 1.1947, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 1.7009, tolerance=.tol[["ci"]])

   ### conditional MLE of the rate ratio
   res <- rma.glmm(x1i=x1i, x2i=x2i, t1i=t1i, t2i=t2i, data=dat, measure="IRR", digits=2, level=90, model="CM.EL", method="EE")

   expect_equivalent(coef(res), 0.3545, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.1779, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 0.5312, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE.Wld, 10.1991, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.Wld, 0.0372, tolerance=.tol[["pval"]])
   expect_equivalent(res$QE.LRT, 12.1324, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.LRT, 0.0164, tolerance=.tol[["pval"]])

   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred,  1.4255, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 1.1947, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 1.7009, tolerance=.tol[["ci"]])

})

############################################################################

### create dataset (Table 15-5)
dat <- data.frame(
age = c("<35", "35+"),
ai = c(3,1),
bi = c(9,3),
ci = c(104,5),
di = c(1059,86),
stringsAsFactors=FALSE)

test_that("the to.table() function works.", {

   tmp <- to.table(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", slab=age, rows=c("Down Syndrome", "Control"), cols=c("Spermicide Use", "No Spermicide"))

   expected <- structure(c(3, 104, 9, 1059, 1, 5, 3, 86), .Dim = c(2L, 2L, 2L), .Dimnames = list(c("Down Syndrome", "Control"), c("Spermicide Use", "No Spermicide"), c("<35", "35+")))

   ### compare with data in Table 15-5
   expect_equivalent(tmp, expected)

})

test_that("the to.long() function works.", {

   tmp <- to.long(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", slab=age)

   expected <- structure(list(age = c("<35", "<35", "35+", "35+"),
      ai = c(3, 3, 1, 1), bi = c(9, 9, 3, 3),
      ci = c(104, 104, 5, 5), di = c(1059, 1059, 86, 86),
      study = structure(c(2L, 2L, 1L, 1L), .Label = c("35+", "<35"), class = "factor"),
      group = structure(c(1L, 2L, 1L, 2L), .Label = c("1", "2"), class = "factor"),
      out1 = c(3, 104, 1, 5), out2 = c(9, 1059, 3, 86)),
      class = "data.frame", row.names = c(NA, 4L))

   expect_equivalent(tmp, expected)

})

test_that("results are correct for Mantel-Haenszel method.", {

   ### Mantel-Haenszel method with odds ratios
   res <- rma.mh(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", digits=2, level=90, correct=FALSE)

   expect_equivalent(coef(res), 1.3300, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.3579, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 2.3021, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE,  0.1378, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0.7105, tolerance=.tol[["pval"]])
   expect_equivalent(res$CO,  5.8248, tolerance=.tol[["test"]])
   expect_equivalent(res$COp, 0.0158, tolerance=.tol[["pval"]])
   expect_equivalent(res$MH,  5.8092, tolerance=.tol[["test"]])
   expect_equivalent(res$MHp, 0.0159, tolerance=.tol[["pval"]])
   expect_equivalent(res$TA,  0.1391, tolerance=.tol[["test"]])
   expect_equivalent(res$TAp, 0.7092, tolerance=.tol[["pval"]])

   tmp <- c(confint(res, transf=exp)$fixed)
   expect_equivalent(tmp, c(3.7812, 1.4304, 9.9954), tolerance=.tol[["ci"]])

   skip_on_cran()

   ### unconditional MLE of the odds ratio
   res <- rma.glmm(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", digits=2, level=90, model="UM.FS", method="EE")

   expect_equivalent(coef(res), 1.3318, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.3582, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 2.3053, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE.Wld,  0.1374, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.Wld, 0.7109, tolerance=.tol[["pval"]])
   expect_equivalent(res$QE.LRT,  0.1324, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.LRT, 0.7160, tolerance=.tol[["pval"]])

   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred,   3.7878, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb,  1.4308, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 10.0276, tolerance=.tol[["ci"]])

   ### conditional MLE of the odds ratio
   #res <- rma.glmm(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", digits=2, level=90, model="CM.EL", method="EE", control=list(optimizer="bobyqa"))
   res <- rma.glmm(ai=ai, bi=bi, ci=ci, di=di, data=dat, measure="OR", digits=2, level=90, model="CM.EL", method="EE")

   expect_equivalent(coef(res), 1.3257, tolerance=.tol[["coef"]])
   expect_equivalent(res$ci.lb, 0.3559, tolerance=.tol[["ci"]])
   expect_equivalent(res$ci.ub, 2.2954, tolerance=.tol[["ci"]])
   expect_equivalent(res$QE.Wld,  0.1327, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.Wld, 0.7156, tolerance=.tol[["pval"]])
   expect_equivalent(res$QE.LRT,  0.1188, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp.LRT, 0.7304, tolerance=.tol[["pval"]])

   tmp <- predict(res, transf=exp)
   expect_equivalent(tmp$pred,  3.7647, tolerance=.tol[["pred"]])
   expect_equivalent(tmp$ci.lb, 1.4274, tolerance=.tol[["ci"]])
   expect_equivalent(tmp$ci.ub, 9.9287, tolerance=.tol[["ci"]])

})

############################################################################

rm(list=ls())
