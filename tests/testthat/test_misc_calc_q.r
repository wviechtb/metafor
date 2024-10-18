### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking misc: computation of Q-test")

source("settings.r")

test_that("computation is correct for 'dat.bcg'.", {

   dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

   res <- rma(yi, vi, data=dat)
   expect_equivalent(res$QE, 152.23301, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0, tolerance=.tol[["pval"]])

   res <- rma(yi, vi, mods = ~ ablat, data=dat)
   expect_equivalent(res$QE, 30.73309, tolerance=.tol[["test"]])
   expect_equivalent(res$QEp, 0.001214, tolerance=.tol[["pval"]])

})

perm <- function(v) {
   n <- length(v)
   if (n == 1) {
      v
   } else {
      X <- NULL
      for (i in 1:n)
         X <- rbind(X, cbind(v[i], perm(v[-i])))
      X
   }
}

test_that("the computation is correct for measurements of the Planck constant.", {

   dat <- read.table(header=TRUE, text="
   exp     h             uh
   NRC-17  6.62607013300 6.00e-08
   NMIJ-17 6.62607005883 1.65e-07
   NIST-17 6.62606993400 8.90e-08")

   perms <- perm(1:nrow(dat))

   QE  <- rep(NA_real_, nrow(dat))
   QEp <- rep(NA_real_, nrow(dat))

   for (i in 1:nrow(perms)) {
       tmp <- dat[perms[i,],]
       res <- rma(yi=h, sei=uh, data=tmp, method="DL")
       QE[i]  <- res$QE
       QEp[i] <- res$QEp
   }

   expect_equivalent(QE, rep(3.442127, length(QE)), tolerance=.tol[["test"]])
   expect_equivalent(QEp, rep(0.1788758, length(QEp)), tolerance=.tol[["pval"]])

})

test_that("the computation is correct for measurements of the Newtonian gravitational constant.", {

   dat <- read.table(header=TRUE, text="
   label         G          uG
   NIST-82       6.67248    0.00043
   TR&D-96       6.6729     0.00050
   LANL-97       6.67398    0.00070
   UWash-00      6.674255   0.000092
   BIPM-01       6.67559    0.00027
   UWup-02       6.67422    0.00098
   MSL-03        6.67387    0.00027
   HUST-05       6.67222    0.00087
   UZur-06       6.67425    0.00012
   HUST-09       6.67349    0.00018
   BIPM-14       6.67554    0.00016
   LENS-14       6.67191    0.00099
   UCI-14        6.67435    0.00013
   HUSTT-18      6.674184   0.000078
   HUSTA-18      6.674484   0.000077
   JILA-18       6.67260    0.00025")

   QE  <- rep(NA_real_, 100)
   QEp <- rep(NA_real_, 100)

   set.seed(1234)

   for (i in 1:100) {
       tmp <- dat[sample(nrow(dat)),]
       res <- rma(yi=G, sei=uG, data=tmp, method="DL")
       QE[i]  <- res$QE
       QEp[i] <- res$QEp
   }

   expect_equivalent(QE, rep(197.8399, length(QE)), tolerance=.tol[["test"]])
   expect_equivalent(QEp, rep(0, length(QEp)), tolerance=.tol[["pval"]])

})

test_that("the computation is correct for measurements Planck constant.", {

   dat <- read.table(header=TRUE, text="
      label           h       uh
     NPL-79 6.626073000 6.70e-06
    NIST-80 6.626065800 8.80e-06
     NMI-89 6.626068400 3.60e-06
     NPL-90 6.626068200 1.30e-06
     PTB-91 6.626067000 4.20e-06
     NIM-95 6.626071000 1.10e-05
    NIST-98 6.626068910 5.80e-07
     IAC-11 6.626069890 2.00e-07
   METAS-11 6.626069100 2.00e-06
     NPL-12 6.626071200 1.30e-06
     IAC-15 6.626070150 1.30e-07
     LNE-15 6.626068800 1.70e-06
    NIST-15 6.626069360 3.80e-07
     NRC-17 6.626070133 6.00e-08
     LNE-17 6.626070410 3.80e-07
    NMIJ-17 6.626070059 1.65e-07
     NIM-17 6.626069200 1.60e-06
     IAC-17 6.626070404 7.92e-08")

   QE  <- rep(NA_real_, 100)
   QEp <- rep(NA_real_, 100)

   set.seed(1234)

   for (i in 1:100) {
       tmp <- dat[sample(nrow(dat)),]
       res <- rma(yi=h, sei=uh, data=tmp, method="DL")
       QE[i]  <- res$QE
       QEp[i] <- res$QEp
   }

   expect_equivalent(QE, rep(26.63226, length(QE)), tolerance=.tol[["test"]])
   expect_equivalent(QEp, rep(0.06368617, length(QEp)), tolerance=.tol[["pval"]])

})

rm(list=ls())
