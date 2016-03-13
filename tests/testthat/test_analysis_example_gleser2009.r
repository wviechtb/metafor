### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

### see also: http://www.metafor-project.org/doku.php/analyses:gleser2009

context("Checking analysis example: gleser2009")

############################################################################

### create dataset

dat <- data.frame(study=c(1,1,2,3,3,3), trt=c(1,2,1,1,2,3),
                  ai=c( 40, 40, 10,150,150,150), n1i=c(1000,1000,200,2000,2000,2000),
                  ci=c(100,150, 15, 40, 80, 50), n2i=c(4000,4000,400,1000,1000,1000))
dat$pti <- with(dat, ci / n2i)
dat$pci <- with(dat, ai / n1i)

test_that("results are correct for the multiple-treatment studies example with risk differences.", {

   dat <- escalc(measure="RD", ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=dat)

   ### compare with results on page 360 (Table 19.2)
   expect_equivalent(round(dat$yi,4), c(0.0150, 0.0025, 0.0125, 0.0350, -0.0050, 0.0250))

   calc.v <- function(x) {
      v <- matrix(x$pci[1]*(1-x$pci[1])/x$n1i[1], nrow=nrow(x), ncol=nrow(x))
      diag(v) <- x$vi
      v
   }

   V <- bldiag(lapply(split(dat, dat$study), calc.v))

   res <- rma.mv(yi, V, mods = ~ factor(trt) - 1, data=dat, method="FE")

   ### compare with results on page 361 (eq. 19.6)
   expect_equivalent(round(coef(res),4), c(0.0200, 0.0043, 0.0211))

   ### compare with results on page 361 (eq. 19.7)
   tmp <- round(vcov(res) * 10^6, 3)
   expected <- structure(c(24.612, 19.954, 13.323, 19.954, 28.538, 13.255, 13.323, 13.255, 69.806), .Dim = c(3L, 3L),
                         .Dimnames = list(c("factor(trt)1", "factor(trt)2", "factor(trt)3"), c("factor(trt)1", "factor(trt)2", "factor(trt)3")))
   expect_equivalent(tmp, expected)

   ### compare with results on page 362 (eq. 19.8)
   expect_equivalent(round(res$QE,3), 7.191)

})

test_that("results are correct for the multiple-treatment studies example with log odds ratios.", {

   dat <- escalc(measure="OR", ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=dat)

   ### compare with results on page 362
   expect_equivalent(round(dat$yi,4), c(0.4855, 0.0671, 0.3008, 0.6657, -0.0700, 0.4321))

   calc.v <- function(x) {
      v <- matrix(1/(x$n1i[1]*x$pci[1]*(1-x$pci[1])), nrow=nrow(x), ncol=nrow(x))
      diag(v) <- x$vi
      v
   }

   V <- bldiag(lapply(split(dat, dat$study), calc.v))

   res <- rma.mv(yi, V, mods = ~ factor(trt) - 1, data=dat, method="FE")

   ### compare with results on page 363
   expect_equivalent(round(coef(res),4), c(0.5099, 0.0044, 0.4301))

   ### compare with results on page 363
   tmp <- round(vcov(res), 5)
   expected <- structure(c(0.01412, 0.00712, 0.00425, 0.00712, 0.01178, 0.00455, 0.00425, 0.00455, 0.02703), .Dim = c(3L, 3L),
                         .Dimnames = list(c("factor(trt)1", "factor(trt)2", "factor(trt)3"), c("factor(trt)1", "factor(trt)2", "factor(trt)3")))
   expect_equivalent(tmp, expected)

   ### compare with results on page 363
   expect_equivalent(round(res$QE,3), 2.056) ### 2.057 in chapter

})

test_that("results are correct for the multiple-treatment studies example with log risk ratios.", {

   dat <- escalc(measure="RR", ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=dat)

   ### compare with results on page 364
   expect_equivalent(round(dat$yi,4), c(0.4700, 0.0645, 0.2877, 0.6286, -0.0645, 0.4055))

   calc.v <- function(x) {
      v <- matrix((1-x$pci[1])/(x$n1i[1]*x$pci[1]), nrow=nrow(x), ncol=nrow(x))
      diag(v) <- x$vi
      v
   }

   V <- bldiag(lapply(split(dat, dat$study), calc.v))

   res <- rma.mv(yi, V, mods = ~ factor(trt) - 1, data=dat, method="FE")

   ### compare with results on page 363
   expect_equivalent(round(coef(res),4), c(0.4875, 0.0006, 0.4047))

   ### (results for this not given in chapter)
   tmp <- round(vcov(res), 5)
   expected <- structure(c(0.01287, 0.00623, 0.00371, 0.00623, 0.01037, 0.00399, 0.00371, 0.00399, 0.02416), .Dim = c(3L, 3L),
                         .Dimnames = list(c("factor(trt)1", "factor(trt)2", "factor(trt)3"), c("factor(trt)1", "factor(trt)2", "factor(trt)3")))
   expect_equivalent(tmp, expected)

   ### (results for this not given in chapter)
   expect_equivalent(round(res$QE,3), 1.895)

})

test_that("results are correct for the multiple-treatment studies example with difference of arcsine transformed risks.", {

   dat <- escalc(measure="AS", ai=ai, ci=ci, n1i=n1i, n2i=n2i, data=dat)

   ### compare with results on page 364
   expect_equivalent(round(dat$yi*2,4), c(0.0852, 0.0130, 0.0613, 0.1521, -0.0187, 0.1038)) ### need *2 factor due to difference in definition of measure

   calc.v <- function(x) {
      v <- matrix(1/(4*x$n1i[1]), nrow=nrow(x), ncol=nrow(x))
      diag(v) <- x$vi
      v
   }

   V <- bldiag(lapply(split(dat, dat$study), calc.v))

   res <- rma.mv(yi, V, mods = ~ factor(trt) - 1, data=dat, method="FE")

   ### compare with results on page 365
   expect_equivalent(round(coef(res)*2,4), c(0.1010, 0.0102, 0.0982))

   ### compare with results on page 365
   tmp <- round(vcov(res)*2^2, 5)
   expected <- structure(c(0.00058, 4e-04, 0.00024, 4e-04, 0.00061, 0.00025, 0.00024, 0.00025, 0.00137), .Dim = c(3L, 3L),
                         .Dimnames = list(c("factor(trt)1", "factor(trt)2", "factor(trt)3"), c("factor(trt)1", "factor(trt)2", "factor(trt)3")))
   expect_equivalent(tmp, expected)

   ### compare with results on page 365
   expect_equivalent(round(res$QE,3), 4.263) ### 4.264 in chapter

})

############################################################################

### create dataset

dat <- data.frame(study=c(1,1,2,3,4,4), trt=c(1,2,1,1,1,2),
                  m1i=c(7.87, 4.35, 9.32, 8.08, 7.44, 5.34),
                  m2i=c(-1.36, -1.36, 0.98, 1.17, 0.45, 0.45),
                  sdpi=c(4.2593,4.2593,2.8831,3.1764,2.9344,2.9344),
                  n1i=c(25,22,38,50,30,30), n2i=c(25,25,40,50,30,30))

test_that("results are correct for the multiple-treatment studies example with standardized mean differences.", {

   dat$Ni <- unlist(lapply(split(dat, dat$study), function(x) rep(sum(x$n1i) + x$n2i[1], each=nrow(x))))

   dat$yi <- with(dat, (m1i-m2i)/sdpi)
   dat$vi <- with(dat, 1/n1i + 1/n2i + yi^2/(2*Ni))

   ### compare with results on page 364
   expect_equivalent(round(dat$yi,4), c(2.1670, 1.3406, 2.8927, 2.1754, 2.3821, 1.6664))

   calc.v <- function(x) {
      v <- matrix(1/x$n2i[1] + outer(x$yi, x$yi, "*")/(2*x$Ni[1]), nrow=nrow(x), ncol=nrow(x))
      diag(v) <- x$vi
      v
   }

   V <- bldiag(lapply(split(dat, dat$study), calc.v))

   res <- rma.mv(yi, V, mods = ~ factor(trt) - 1, data=dat, method="FE")

   ### compare with results on page 367
   expect_equivalent(round(coef(res),3), c(2.374, 1.570))

   ### compare with results on page 367
   tmp <- round(vcov(res), 5)
   expected <- structure(c(0.02257, 0.01244, 0.01244, 0.03554), .Dim = c(2L, 2L), .Dimnames = list(c("factor(trt)1", "factor(trt)2"), c("factor(trt)1", "factor(trt)2")))
   expect_equivalent(tmp, expected)

   ### compare with results on page 367
   expect_equivalent(round(res$QE,3), 3.945)

})

############################################################################

### create dataset

dat <- data.frame(school=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7),
                  outcome=rep(c("math", "reading"), times=7),
                  m1i=c(2.3,2.5,2.4,1.3,2.5,2.4,3.3,1.7,1.1,2.0,2.8,2.1,1.7,0.6),
                  m2i=c(10.3,6.6,9.7,3.1,8.7,3.7,7.5,8.5,2.2,2.1,3.8,1.4,1.8,3.9),
                  sdpi=c(8.2,7.3,8.3,8.9,8.5,8.3,7.7,9.8,9.1,10.4,9.6,7.9,9.2,10.2),
                  ri=rep(c(.55,.43,.57,.66,.51,.59,.49), each=2),
                  n1i=rep(c(22,21,26,18,38,42,39), each=2),
                  n2i=rep(c(24,21,23,18,36,42,38), each=2))

test_that("results are correct for the multiple-endpoint studies example with standardized mean differences.", {

   dat$yi <- round(with(dat, (m2i-m1i)/sdpi), 3)
   dat$vi <- round(with(dat, 1/n1i + 1/n2i + yi^2/(2*(n1i+n2i))), 4)
   dat$covi <- round(with(dat, (1/n1i + 1/n2i) * ri + (rep(sapply(split(dat$yi, dat$school), prod), each=2) / (2*(n1i+n2i))) * ri^2), 4)

   V <- bldiag(lapply(split(dat, dat$school), function(x) matrix(c(x$vi[1], x$covi[1], x$covi[2], x$vi[2]), nrow=2)))

   ### fit model

   res <- rma.mv(yi, V, mods = ~ outcome - 1, data=dat, method="FE")

   ### (results for this not given in chapter)
   expect_equivalent(round(coef(res),3), c(0.362, 0.205))

   ### (results for this not given in chapter)
   tmp <- round(vcov(res), 5)
   expected <- structure(c(0.01008, 0.00537, 0.00537, 0.00989), .Dim = c(2L, 2L), .Dimnames = list(c("outcomemath", "outcomereading"), c("outcomemath", "outcomereading")))
   expect_equivalent(tmp, expected)

   ### compare with results on page 371
   expect_equivalent(round(res$QE,2), 19.63) ### 19.62 in chapter

})

############################################################################
