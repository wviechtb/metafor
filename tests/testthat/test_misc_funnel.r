### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: funnel() functions")

source("settings.r")

test_that("funnel() works correctly.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   ### simulate a large meta-analytic dataset (correlations with rho = 0.0)
   ### with no heterogeneity or publication bias; then try out different
   ### versions of the funnel plot

   gencor <- function(rhoi, ni) {
      x1 <- rnorm(ni, mean=0, sd=1)
      x2 <- rnorm(ni, mean=0, sd=1)
      x3 <- rhoi*x1 + sqrt(1-rhoi^2)*x2
      cor(x1, x3)
   }

   set.seed(78123)

   k  <- 200                               ### number of studies to simulate
   ni <- round(rchisq(k, df=2) * 20 + 20)  ### simulate sample sizes (skewed distribution)
   ri <- mapply(gencor, rep(0.0,k), ni)    ### simulate correlations

   dat <- escalc(measure="ZCOR", ri=ri, ni=ni) ### compute r-to-z transformed correlations
   res <- rma(yi, vi, data=dat, method="EE")

   png(filename="test_misc_funnel_1.png", res=200, width=1800, height=2000, type="cairo")

   par(mfrow=c(5,2), mar=c(5,4,1,1), cex=0.5)

   funnel(res, yaxis="sei")
   funnel(res, yaxis="vi")
   funnel(res, yaxis="seinv")
   funnel(res, yaxis="vinv")
   funnel(res, yaxis="ni")
   funnel(res, yaxis="ninv")
   funnel(res, yaxis="sqrtni")
   funnel(res, yaxis="sqrtninv")
   funnel(res, yaxis="lni")
   funnel(res, yaxis="wi")

   dev.off()

   expect_true(.vistest("test_misc_funnel_1.png", "images/test_misc_funnel_1.png"))

   png(filename="test_misc_funnel_2.png", res=200, width=1800, height=2000, type="cairo")

   par(mfrow=c(5,2), mar=c(5,4,1,1), cex=0.5)

   funnel(dat$yi, dat$vi, yaxis="sei")
   funnel(dat$yi, dat$vi, yaxis="vi")
   funnel(dat$yi, dat$vi, yaxis="seinv")
   funnel(dat$yi, dat$vi, yaxis="vinv")
   funnel(dat$yi, dat$vi, yaxis="ni")
   funnel(dat$yi, dat$vi, yaxis="ninv")
   funnel(dat$yi, dat$vi, yaxis="sqrtni")
   funnel(dat$yi, dat$vi, yaxis="sqrtninv")
   funnel(dat$yi, dat$vi, yaxis="lni")
   funnel(dat$yi, dat$vi, yaxis="wi")

   dev.off()

   expect_true(.vistest("test_misc_funnel_2.png", "images/test_misc_funnel_2.png"))

})

rm(list=ls())
