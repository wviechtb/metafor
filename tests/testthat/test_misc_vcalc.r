### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: vcalc() function")

source("settings.r")

test_that("vcov() works correctly for 'dat.assink2016' example.", {

   dat <- dat.assink2016

   ### assume that the effect sizes within studies are correlated with rho=0.6
   V <- vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)
   sav <- blsplit(V, dat$study, round, 4)[1:2]

   expected <- list(`1` = structure(c(0.074, 0.0326, 0.0358, 0.0252, 0.0297, 0.0486, 0.0326, 0.0398, 0.0263, 0.0185, 0.0218, 0.0356, 0.0358, 0.0263, 0.0481, 0.0203, 0.0239, 0.0392, 0.0252, 0.0185, 0.0203, 0.0239, 0.0169, 0.0276, 0.0297, 0.0218, 0.0239, 0.0169, 0.0331, 0.0325, 0.0486, 0.0356, 0.0392, 0.0276, 0.0325, 0.0886), .Dim = c(6L,6L)), `2` = structure(c(0.0115, 0.0056, 0.0052, 0.0056, 0.0076, 0.0042, 0.0052, 0.0042, 0.0065), .Dim = c(3L, 3L)))
   expect_equivalent(sav, expected, tolerance=.tol[["var"]])

   ### use a correlation of 0.7 for effect sizes corresponding to the same type of
   ### delinquent behavior and a correlation of 0.5 for effect sizes corresponding
   ### to different types of delinquent behavior
   V <- vcalc(vi, cluster=study, type=deltype, obs=esid, data=dat, rho=c(0.7, 0.5))
   sav <- blsplit(V, dat$study, round, 3)[16]

   expected <- list(`16` = structure(c(0.091, 0.045, 0.027, 0.044, 0.03, 0.039, 0.076, 0.028, 0.034, 0.03, 0.039, 0.043, 0.039, 0.067, 0.028, 0.032, 0.045, 0.087, 0.027, 0.061, 0.03, 0.039, 0.053, 0.027, 0.047, 0.041, 0.053, 0.059, 0.053, 0.046, 0.038, 0.043, 0.027, 0.027, 0.033, 0.027, 0.025, 0.033, 0.033, 0.023, 0.021, 0.018, 0.023, 0.026, 0.023, 0.029, 0.017, 0.019, 0.044, 0.061, 0.027, 0.086, 0.029, 0.038, 0.053, 0.027, 0.047, 0.041, 0.053, 0.058, 0.053, 0.046, 0.038, 0.043, 0.03, 0.03, 0.025, 0.029, 0.04, 0.037, 0.036, 0.026, 0.023, 0.02, 0.026, 0.028, 0.026, 0.031, 0.018, 0.021, 0.039, 0.039, 0.033, 0.038, 0.037, 0.068, 0.047, 0.033, 0.03, 0.026, 0.034, 0.037, 0.034, 0.041, 0.024, 0.027, 0.076, 0.053, 0.033, 0.053, 0.036, 0.047, 0.129, 0.033, 0.041, 0.035, 0.046, 0.051, 0.046, 0.079, 0.033, 0.037, 0.028, 0.027, 0.023, 0.027, 0.026, 0.033, 0.033, 0.033, 0.021, 0.018, 0.023, 0.026, 0.024, 0.029, 0.017, 0.019, 0.034, 0.047, 0.021, 0.047, 0.023, 0.03, 0.041, 0.021, 0.052, 0.031, 0.041, 0.045, 0.041, 0.036, 0.029, 0.033, 0.03, 0.041, 0.018, 0.041, 0.02, 0.026, 0.035, 0.018, 0.031, 0.039, 0.036, 0.039, 0.036, 0.031, 0.025, 0.029, 0.039, 0.053, 0.023, 0.053, 0.026, 0.034, 0.046, 0.023, 0.041, 0.036, 0.066, 0.051, 0.047, 0.04, 0.033, 0.038, 0.043, 0.059, 0.026, 0.058, 0.028, 0.037, 0.051, 0.026, 0.045, 0.039, 0.051, 0.081, 0.051, 0.045, 0.037, 0.042, 0.039, 0.053, 0.023, 0.053, 0.026, 0.034, 0.046, 0.024, 0.041, 0.036, 0.047, 0.051, 0.067, 0.041, 0.033, 0.038, 0.067, 0.046, 0.029, 0.046, 0.031, 0.041, 0.079, 0.029, 0.036, 0.031, 0.04, 0.045, 0.041, 0.099, 0.029, 0.033, 0.028, 0.038, 0.017, 0.038, 0.018, 0.024, 0.033, 0.017, 0.029, 0.025, 0.033, 0.037, 0.033, 0.029, 0.034, 0.027, 0.032, 0.043, 0.019, 0.043, 0.021, 0.027, 0.037, 0.019, 0.033, 0.029, 0.038, 0.042, 0.038, 0.033, 0.027, 0.044), .Dim = c(16L, 16L)))
   expect_equivalent(sav, expected, tolerance=.tol[["var"]])

})

test_that("vcov() works correctly for 'dat.ishak2007' example.", {

   dat <- dat.ishak2007

   ### create long format dataset
   dat <- reshape(dat, direction="long", idvar="study", v.names=c("yi","vi"),
                  varying=list(c(2,4,6,8), c(3,5,7,9)))
   dat <- dat[order(study, time),]

   ### remove missing measurement occasions from dat
   dat <- dat[!is.na(yi),]
   rownames(dat) <- NULL

   ### construct the full (block diagonal) V matrix with an AR(1) structure
   ### assuming an autocorrelation of 0.97 as estimated by Ishak et al. (2007)
   V <- vcalc(vi, cluster=study, time1=time, phi=0.97, data=dat)
   sav <- blsplit(V, dat$study)[1:5]

   expected <- list(`Alegret (2001)` = structure(14.3, .Dim = c(1L, 1L)), `Barichella (2003)` = structure(c(7.3, 6.0693520102314, 6.0693520102314, 5.7), .Dim = c(2L, 2L)), `Berney (2002)` = structure(7.3, .Dim = c(1L, 1L)), `Burchiel (1999)` = structure(c(8, 7.76, 5.95077410090486, 7.76, 8, 6.13481866072665, 5.95077410090486, 6.13481866072665, 5), .Dim = c(3L, 3L)), `Chen (2003)` = structure(125, .Dim = c(1L, 1L)))
   expect_equivalent(sav, expected, tolerance=.tol[["var"]])

})

test_that("vcov() works correctly for 'dat.kalaian1996' example.", {

   dat <- dat.kalaian1996

   ### construct the variance-covariance matrix assuming rho = 0.66 for effect sizes
   ### corresponding to the 'verbal' and 'math' outcome types
   V <- vcalc(vi, cluster=study, type=outcome, data=dat, rho=0.66)
   sav <- round(V[1:12,1:12], 4)

   expected <- structure(c(0.0817, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0507, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1045, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0442, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0557, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0561, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1151, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0147, 0.0097, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0097, 0.0147, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0218, 0.0143, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0143, 0.0216), .Dim = c(12L, 12L))
   expect_equivalent(sav, expected, tolerance=.tol[["var"]])

})

test_that("vcov() works correctly for 'dat.berkey1998' example.", {

   dat <- dat.berkey1998

   ### variables v1i and v2i correspond to the 2x2 var-cov matrices of the studies;
   ### so use these variables to construct the V matrix (note: since v1i and v2i are
   ### var-cov matrices and not correlation matrices, set vi=1 for all rows)
   V <- vcalc(vi=1, cluster=author, rvars=c(v1i, v2i), data=dat)
   sav <- blsplit(V, dat$author, function(x) round(cov2cor(x), 2))

   expected <- list(`Pihlstrom et al.` = structure(c(1, 0.39, 0.39, 1), .Dim = c(2L, 2L)), `Lindhe et al.` = structure(c(1, 0.42, 0.42, 1), .Dim = c(2L, 2L)), `Knowles et al.` = structure(c(1, 0.41, 0.41, 1), .Dim = c(2L, 2L)), `Ramfjord et al.` = structure(c(1, 0.43, 0.43, 1), .Dim = c(2L, 2L)), `Becker et al.` = structure(c(1, 0.34, 0.34, 1), .Dim = c(2L, 2L)))
   expect_equivalent(sav, expected, tolerance=.tol[["var"]])

})

test_that("vcov() works correctly for 'dat.knapp2017' example.", {

   dat <- dat.knapp2017

   ### create variable that indicates the task and difficulty combination as increasing integers
   dat$task.diff <- unlist(lapply(split(dat, dat$study), function(x) {
      task.int <- as.integer(factor(x$task))
      diff.int <- as.integer(factor(x$difficulty))
      diff.int[is.na(diff.int)] <- 1
      paste0(task.int, ".", diff.int)}))

   ### construct correlation matrix for two tasks with four different difficulties where the
   ### correlation is 0.4 for different difficulties of the same task, 0.7 for the same
   ### difficulty of different tasks, and 0.28 for different difficulties of different tasks
   R <- matrix(0.4, nrow=8, ncol=8)
   R[5:8,1:4] <- R[1:4,5:8] <- 0.28
   diag(R[1:4,5:8]) <- 0.7
   diag(R[5:8,1:4]) <- 0.7
   diag(R) <- 1
   rownames(R) <- colnames(R) <- paste0(rep(1:2, each=4), ".", 1:4)

   ### construct an approximate V matrix accounting for the use of shared groups and
   ### for correlations among tasks/difficulties as specified in the R matrix above
   V <- vcalc(vi, cluster=study, grp1=group1, grp2=group2, w1=n_sz, w2=n_hc,
              obs=task.diff, rho=R, data=dat)
   Vs <- blsplit(V, dat$study)
   sav <- Vs[c(3,6,12,24,29)]

   expected <- list(`3` = structure(c(0.062, 0.0313021866879515, 0.0305960523769429, 0.0306223534669685, 0.0313021866879515, 0.073, 0.0301021398261882, 0.0301280163373072, 0.0305960523769429, 0.0301021398261882, 0.102, 0.029448369695669, 0.0306223534669685, 0.0301280163373072, 0.029448369695669, 0.084), .Dim = c(4L, 4L)), `6` = structure(c(0.17, 0.07485452558129, 0.0675988165576883, 0.0711280535372648, 0.120045408075445, 0.0489799959167005, 0.0511105468567888, 0.0495212277715325, 0.07485452558129, 0.206, 0.0744129021070943, 0.0782978926919493, 0.0528584827629398, 0.134793174901402, 0.0562625843700767, 0.0545130589858981, 0.0675988165576884, 0.0744129021070943, 0.168, 0.0707084153407499, 0.0477348677593224, 0.0486910258671965, 0.127022517688794, 0.0492290645858725, 0.0711280535372648, 0.0782978926919493, 0.0707084153407499, 0.186, 0.0502270365440765, 0.0512331142914424, 0.0534616722521846, 0.129498108094288, 0.120045408075445, 0.0528584827629398, 0.0477348677593224, 0.0502270365440765, 0.173, 0.0705861176152932, 0.0736565000526091, 0.0713660983941255, 0.0489799959167005, 0.134793174901402, 0.0486910258671965, 0.0512331142914424, 0.0705861176152932, 0.18, 0.0751318840439929, 0.0727956042628949, 0.0511105468567888, 0.0562625843700767, 0.127022517688794, 0.0534616722521846, 0.0736565000526091, 0.0751318840439929, 0.196, 0.075962095811003, 0.0495212277715325, 0.0545130589858981, 0.0492290645858725, 0.129498108094288, 0.0713660983941255, 0.0727956042628949, 0.075962095811003, 0.184), .Dim = c(8L, 8L)), `12` = structure(c(0.02, 0.00819756061276768, 0.008, 0.00839047078536121, 0.00819756061276768, 0.021, 0.00819756061276768, 0.00859767410408187, 0.008, 0.00819756061276768, 0.02, 0.00839047078536121, 0.00839047078536121, 0.00859767410408187, 0.00839047078536121, 0.022), .Dim = c(4L, 4L)), `24` = structure(c(0.022, 0, 0, 0.03), .Dim = c(2L, 2L)), `29` = structure(c(0.039, 0, 0, 0, 0.039, 0, 0, 0, 0.121), .Dim = c(3L, 3L)))
   expect_equivalent(sav, expected, tolerance=.tol[["var"]])

})

rm(list=ls())
