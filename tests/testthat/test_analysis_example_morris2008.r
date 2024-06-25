### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

### see: https://www.metafor-project.org/doku.php/analyses:morris2008

context("Checking analysis example: morris2008")

source("settings.r")

### create datasets

datT <- data.frame(
m_pre   = c(30.6, 23.5, 0.5, 53.4, 35.6),
m_post  = c(38.5, 26.8, 0.7, 75.9, 36.0),
sd_pre  = c(15.0, 3.1, 0.1, 14.5, 4.7),
sd_post = c(11.6, 4.1, 0.1, 4.4, 4.6),
ni      = c(20, 50, 9, 10, 14),
ri      = c(.47, .64, .77, .89, .44))

datC <- data.frame(
m_pre   = c(23.1, 24.9, 0.6, 55.7, 34.8),
m_post  = c(19.7, 25.3, 0.6, 60.7, 33.4),
sd_pre  = c(13.8, 4.1, 0.2, 17.3, 3.1),
sd_post = c(14.8, 3.3, 0.2, 17.9, 6.9),
ni      = c(20, 42, 9, 11, 14),
ri      = c(.47, .64, .77, .89, .44))

test_that("calculations of escalc() are correct for measure='SMCR'.", {

   ### compute standardized mean changes using raw-score standardization
   datT <- escalc(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ni=ni, ri=ri, data=datT)
   datC <- escalc(measure="SMCR", m1i=m_post, m2i=m_pre, sd1i=sd_pre, ni=ni, ri=ri, data=datC)

   ### (results for this not given in paper)
   expect_equivalent(datT$yi, c( 0.5056, 1.0481, 1.8054, 1.4181,  0.0801), tolerance=.tol[["est"]])
   expect_equivalent(datT$vi, c( 0.0594, 0.0254, 0.2322, 0.1225,  0.0802), tolerance=.tol[["var"]])
   expect_equivalent(datC$yi, c(-0.2365, 0.0958, 0.0000, 0.2667, -0.4250), tolerance=.tol[["est"]])
   expect_equivalent(datC$vi, c( 0.0544, 0.0173, 0.0511, 0.0232,  0.0864), tolerance=.tol[["var"]])

   ### compute difference between treatment and control groups
   dat <- data.frame(yi = datT$yi - datC$yi, vi = datT$vi + datC$vi)

   ### compare with results on page 382 (Table 5)
   expect_equivalent(dat$yi, c(0.7421, 0.9524, 1.8054, 1.1514, 0.5050), tolerance=.tol[["est"]])

   ### (results for this not given in paper)
   expect_equivalent(dat$vi, c(0.1138, 0.0426, 0.2833, 0.1458, 0.1667), tolerance=.tol[["var"]])

   ### use pooled pretest SDs

   sd_pool <- sqrt((with(datT, (ni-1)*sd_pre^2) + with(datC, (ni-1)*sd_pre^2)) / (datT$ni + datC$ni - 2))

   dat <- data.frame(yi = metafor:::.cmicalc(datT$ni + datC$ni - 2) * (with(datT, m_post - m_pre) - with(datC, m_post - m_pre)) / sd_pool)
   dat$vi <- 2*(1-datT$ri) * (1/datT$ni + 1/datC$ni) + dat$yi^2 / (2*(datT$ni + datC$ni))

   ### compare with results on page 382 (Table 5)
   expect_equivalent(dat$yi, c(0.7684, 0.8010, 1.2045, 1.0476, 0.4389), tolerance=.tol[["est"]])

   ### (results for this not given in paper)
   expect_equivalent(dat$vi, c(0.1134, 0.0350, 0.1425, 0.0681, 0.1634), tolerance=.tol[["var"]])

})

test_that("calculations of escalc() are correct for measure='SMCC'.", {

   ### compute standardized mean changes using change-score standardization
   datT <- escalc(measure="SMCC", m1i=m_post, m2i=m_pre, sd1i=sd_post, sd2i=sd_pre, ni=ni, ri=ri, data=datT)
   datC <- escalc(measure="SMCC", m1i=m_post, m2i=m_pre, sd1i=sd_post, sd2i=sd_pre, ni=ni, ri=ri, data=datC)

   ### (results for this not given in paper)
   expect_equivalent(datT$yi, c( 0.5417, 1.0198, 2.6619, 1.9088,  0.0765), tolerance=.tol[["est"]])
   expect_equivalent(datT$vi, c( 0.0573, 0.0304, 0.5048, 0.2822,  0.0716), tolerance=.tol[["var"]])
   expect_equivalent(datC$yi, c(-0.2213, 0.1219, 0.0000, 0.5575, -0.2126), tolerance=.tol[["est"]])
   expect_equivalent(datC$vi, c( 0.0512, 0.0240, 0.1111, 0.1050,  0.0730), tolerance=.tol[["var"]])

   ### compute difference between treatment and control groups
   dat <- data.frame(yi = datT$yi - datC$yi, vi = datT$vi + datC$vi)

   ### (results for this not given in paper)
   expect_equivalent(dat$yi, c(0.7630, 0.8979, 2.6619, 1.3513, 0.2891), tolerance=.tol[["est"]])
   expect_equivalent(dat$vi, c(0.1086, 0.0544, 0.6159, 0.3872, 0.1447), tolerance=.tol[["var"]])

})

rm(list=ls())
