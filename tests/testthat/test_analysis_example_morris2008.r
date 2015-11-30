### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:morris2008

context("Checking analysis example morris2008")

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
   expect_equivalent(round(datT$yi,4), c( 0.5056, 1.0481, 1.8054, 1.4181,  0.0801))
   expect_equivalent(round(datT$vi,4), c( 0.0594, 0.0254, 0.2322, 0.1225,  0.0802))
   expect_equivalent(round(datC$yi,4), c(-0.2365, 0.0958, 0.0000, 0.2667, -0.4250))
   expect_equivalent(round(datC$vi,4), c( 0.0544, 0.0173, 0.0511, 0.0232,  0.0864))

   ### compute difference between treatment and control groups
   dat <- data.frame(yi = datT$yi - datC$yi, vi = datT$vi + datC$vi)

   ### compare with results on page 368 (Table 1)
   expect_equivalent(round(dat$yi,2), c(0.74, 0.95, 1.81, 1.15, 0.51))

   ### (results for this not given in paper)
   expect_equivalent(round(dat$vi,4), c(0.1138, 0.0426, 0.2833, 0.1458, 0.1667))

})

test_that("calculations of escalc() are correct for measure='SMCC'.", {

   ### compute standardized mean changes using change-score standardization
   datT <- escalc(measure="SMCC", m1i=m_post, m2i=m_pre, sd1i=sd_post, sd2i=sd_pre, ni=ni, ri=ri, data=datT)
   datC <- escalc(measure="SMCC", m1i=m_post, m2i=m_pre, sd1i=sd_post, sd2i=sd_pre, ni=ni, ri=ri, data=datC)

   ### (results for this not given in paper)
   expect_equivalent(round(datT$yi,4), c( 0.5417, 1.0198, 2.6619, 1.9088,  0.0765))
   expect_equivalent(round(datT$vi,4), c( 0.0573, 0.0304, 0.5048, 0.2822,  0.0716))
   expect_equivalent(round(datC$yi,4), c(-0.2213, 0.1219, 0.0000, 0.5575, -0.2126))
   expect_equivalent(round(datC$vi,4), c( 0.0512, 0.0240, 0.1111, 0.1050,  0.0730))

   ### compute difference between treatment and control groups
   dat <- data.frame(yi = datT$yi - datC$yi, vi = datT$vi + datC$vi)

   ### (results for this not given in paper)
   expect_equivalent(round(dat$yi,2), c(0.76, 0.90, 2.66, 1.35, 0.29))
   expect_equivalent(round(dat$vi,4), c(0.1086, 0.0544, 0.6159, 0.3872, 0.1447))

})
