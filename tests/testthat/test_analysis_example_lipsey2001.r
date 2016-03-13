### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:lipsey2001

context("Checking analysis example: lipsey2001")

### create dataset
dat <- data.frame(
id = c(100, 308, 1596, 2479, 9021, 9028, 161, 172, 537, 7049),
yi = c(-0.33, 0.32, 0.39, 0.31, 0.17, 0.64, -0.33, 0.15, -0.02, 0.00),
vi = c(0.084, 0.035, 0.017, 0.034, 0.072, 0.117, 0.102, 0.093, 0.012, 0.067),
random = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
intensity = c(7, 3, 7, 5, 7, 7, 4, 4, 5, 6))

test_that("results are correct for the fixed-effects model.", {

   res <- rma(yi, vi, data=dat, method="FE")

   ### compare with results on page 133 (Exhibit 7.3)
   expect_equivalent(round(c(as.matrix(coef(summary(res)))), 4), c(0.1549, 0.0609, 2.5450, 0.0109, 0.0356, 0.2742))
   expect_equivalent(round(res$QE,4), 14.7640)
   expect_equivalent(round(res$QEp,4), 0.0976)

})

test_that("results are correct for the random-effects model.", {

   res <- rma(yi, vi, data=dat, method="DL")

   ### compare with results on page 133 (Exhibit 7.3)
   expect_equivalent(round(c(as.matrix(coef(summary(res)))), 4), c(0.1534, 0.0858, 1.7893, 0.0736, -0.0146, 0.3215))
   expect_equivalent(round(res$tau2,6), 0.025955)

})

test_that("results are correct for the ANOVA-type analysis.", {

   res <- rma(yi, vi, mods = ~ random, data=dat, method="FE")

   res0 <- rma(yi, vi, data=dat, method="FE", subset=random==0)
   res1 <- rma(yi, vi, data=dat, method="FE", subset=random==1)

   tmp <- predict(res, newmods=c(0,1))
   tmp <- do.call(cbind, unclass(tmp)[1:4])
   tmp <- round(tmp,4)

   ### compare with results on page 138 (Exhibit 7.4)
   expect_equivalent(tmp[1,], c( 0.2984, 0.0813,  0.1390, 0.4578))
   expect_equivalent(tmp[2,], c(-0.0277, 0.0917, -0.2075, 0.1521))
   expect_equivalent(round(res$QM,4), 7.0739) ### 7.0738 in chapter
   expect_equivalent(round(res$QMp,4), 0.0078)
   expect_equivalent(round(res$QE,4), 7.6901) ### 7.6902 in chapter
   expect_equivalent(round(res$QEp,4), 0.4643)
   expect_equivalent(round(res0$QE,4), 6.4382) ### 6.4383 in chapter
   expect_equivalent(round(res0$QEp,4), 0.2659)
   expect_equivalent(round(res1$QE,4), 1.2519)
   expect_equivalent(round(res1$QEp,4), 0.7406)

})

test_that("results are correct for the meta-regression analysis (fixed-effects with moderators model).", {

   res <- rma(yi, vi, mods = ~ random + intensity, data=dat, method="FE")

   expected <- structure(list(estimate = c(0.3223, -0.3298, -0.0041), se = c(0.2998, 0.1304, 0.0493),
                              zval = c(1.0752, -2.5286, -0.0829), pval = c(0.2823, 0.0115, 0.9339),
                              ci.lb = c(-0.2652, -0.5854, -0.1007), ci.ub = c(0.9099, -0.0742, 0.0925)),
                              .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub"),
                              row.names = c("intrcpt", "random", "intensity"), class = "data.frame")

   ### compare with results on page 141 (Exhibit 7.6)

   expect_equivalent(round(coef(summary(res)),4), expected)
   expect_equivalent(round(res$QM,4), 7.0807)
   expect_equivalent(round(res$QMp,4), 0.0290)
   expect_equivalent(round(res$QE,4), 7.6832) ### 7.6833 in chapter
   expect_equivalent(round(res$QEp,4), 0.3614) ### 0.3613 in chapter

})

test_that("results are correct for the meta-regression analysis (mixed-effects model).", {

   res <- rma(yi, vi, mods = ~ random + intensity, data=dat, method="DL")

   expected <- structure(list(estimate = c(0.3311, -0.3269, -0.0068), se = c(0.3198, 0.1439, 0.0528),
                              zval = c(1.0351, -2.2712, -0.1292), pval = c(0.3006, 0.0231, 0.8972),
                              ci.lb = c(-0.2958, -0.609, -0.1103), ci.ub = c(0.9579, -0.0448, 0.0967)),
                              .Names = c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub"),
                              row.names = c("intrcpt", "random", "intensity"), class = "data.frame")

   ### compare with results on page 141 (Exhibit 7.7)

   expect_equivalent(round(coef(summary(res)),4), expected)
   expect_equivalent(round(res$QM,4), 5.5711) ### 5.5709 in chapter
   expect_equivalent(round(res$QMp,4), 0.0617)
   expect_equivalent(round(res$tau2,5), 0.00488) ### 5.5709 in chapter

})

test_that("results are correct for the comutation of R^2 via the anova() function.", {

   res.ME <- rma(yi, vi, mods = ~ random + intensity, data=dat, method="DL")
   res.RE <- rma(yi, vi, data=dat, method="DL")
   tmp <- anova(res.RE, res.ME)

   expect_equivalent(tmp$R2, 81.2)

})
