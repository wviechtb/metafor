### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: https://www.metafor-project.org/doku.php/plots:baujat_plot

context("Checking plots example: Baujat plot")

test_that("plot can be drawn.", {

   expect_equivalent(TRUE, TRUE) # avoid 'Empty test' message

   skip_on_cran()

   opar <- par(no.readonly=TRUE)

   ### adjust margins so the space is better used
   par(mar=c(5,4,2,2))

   ### load data from Pignon et al. (2000)
   dat <- dat.pignon2000

   ### compute estimated log hazard ratios and sampling variances
   dat$yi <- with(dat, OmE/V)
   dat$vi <- with(dat, 1/V)

   ### meta-analysis based on all 65 trials
   res <- rma(yi, vi, data=dat, method="EE", slab=id)

   ### create Baujat plot
   baujat(res, xlim=c(0,20), ylim=c(0,.20))

   par(opar)

})
