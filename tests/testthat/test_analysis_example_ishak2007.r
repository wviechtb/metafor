### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking analysis example: ishak2007")

source("settings.r")

### load dataset
dat <- dat.ishak2007

### create long format dataset
dat.long <- reshape(dat, direction="long", idvar="study", v.names=c("yi","vi"),
                    varying=list(c(2,4,6,8), c(3,5,7,9)))
dat.long <- dat.long[order(dat.long$study, dat.long$time),]
rownames(dat.long) <- 1:nrow(dat.long)

### remove missing measurement occasions from dat.long
is.miss  <- is.na(dat.long$yi)
dat.long <- dat.long[!is.miss,]

### construct the full (block diagonal) V matrix with an AR(1) structure
rho.within <- .97 ### value as estimated by Ishak et al. (2007)
V <- lapply(split(with(dat, cbind(v1i, v2i, v3i, v4i)), dat$study), diag)
V <- lapply(V, function(v) sqrt(v) %*% toeplitz(ARMAacf(ar=rho.within, lag.max=3)) %*% sqrt(v))
V <- bldiag(V)
V <- V[!is.miss,!is.miss] ### remove missing measurement occasions from V

test_that("results are correct for diag(V) and struct='DIAG'.", {

   res <- rma.mv(yi, diag(V), mods = ~ factor(time) - 1, random = ~ factor(time) | study,
                 struct = "DIAG", data = dat.long, sparse=.sparse)

   ### Table 1, column "Time-specific (Independence)"
   expect_equivalent(coef(res), c(-24.8686, -27.4728, -28.5239, -24.1415), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, c(23.0537, 27.8113, 27.6767, 29.9405), tolerance=.tol[["var"]])

})

test_that("results are correct for diag(V) and random study effects.", {

   res <- rma.mv(yi, diag(V), mods = ~ factor(time) - 1, random = ~ 1 | study, data = dat.long, sparse=.sparse)

   ### Table 1, column "Random study effects"
   expect_equivalent(coef(res), c(-26.2127, -27.1916, -28.5464, -25.6339), tolerance=.tol[["coef"]])
   expect_equivalent(res$sigma2, 26.6829, tolerance=.tol[["var"]])

})

test_that("results are correct for diag(V) and struct='ID'.", {

   res <- rma.mv(yi, diag(V), mods = ~ factor(time) - 1, random = ~ factor(time) | study,
                 struct = "ID", data = dat.long, sparse=.sparse)

   ### not in paper
   expect_equivalent(coef(res), c(-24.8792, -27.4670, -28.5185, -24.1502), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, 26.6847, tolerance=.tol[["var"]])

})

test_that("results are correct for diag(V) and struct='HAR'.", {

   res <- rma.mv(yi, diag(V), mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "HAR", data = dat.long, sparse=.sparse)

   ### Table 1, column "Correlated random time effects"
   expect_equivalent(coef(res), c(-25.9578, -27.3100, -28.5543, -25.7923), tolerance=.tol[["coef"]]) # -27.5 in Table vs -27.3
   expect_equivalent(res$tau2, c(20.3185, 35.9720, 26.4233, 30.1298), tolerance=.tol[["var"]]) # 20.4 in Table vs 20.3
   expect_equivalent(res$rho, 1.0000, tolerance=.tol[["cor"]])

})

test_that("results are correct for struct='HAR'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "HAR", data = dat.long, sparse=.sparse)

   ### Table 1, column "Multivariate model"
   expect_equivalent(coef(res), c(-25.9047, -27.4608, -28.6559, -26.4934), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, c(22.7258, 33.7295, 26.1426, 31.1803), tolerance=.tol[["var"]]) # 22.6 in Table vs 22.7; 31.1 in Table vs 31.2
   expect_equivalent(res$rho, 0.8832, tolerance=.tol[["cor"]])

})

test_that("results are correct for struct='AR'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "AR", data = dat.long, sparse=.sparse)

   ### not in paper
   expect_equivalent(coef(res), c(-25.9418, -27.3937, -28.7054, -26.3970), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, 26.6874, tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 0.8656, tolerance=.tol[["cor"]])

})

test_that("results are correct for struct='HCS'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ factor(time) | study,
                 struct = "HCS", data = dat.long, sparse=.sparse)

   ### not in paper
   expect_equivalent(coef(res), c(-25.8814, -27.3293, -28.6510, -26.6631), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, c(20.8629, 32.7429, 27.6593, 32.1908), tolerance=.tol[["var"]])

})

test_that("results are correct for struct='CAR'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "CAR", data = dat.long, sparse=.sparse)

   ### not in paper
   expect_equivalent(coef(res), c(-25.9418, -27.3937, -28.7054, -26.3970), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, 26.6875, tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 0.8656, tolerance=.tol[["cor"]])

})

test_that("results are correct for struct='CAR' with unequally spaced time points.", {

   dat.long$time[dat.long$time == 4] <- 24/3
   dat.long$time[dat.long$time == 3] <- 12/3
   dat.long$time[dat.long$time == 2] <-  6/3
   dat.long$time[dat.long$time == 1] <-  3/3

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "CAR", data = dat.long, sparse=.sparse)

   ### not in paper
   expect_equivalent(coef(res), c(-26.0293, -27.3838, -28.7339, -26.0515), tolerance=.tol[["coef"]])
   expect_equivalent(res$tau2, 26.9825, tolerance=.tol[["var"]])
   expect_equivalent(res$rho, 0.9171, tolerance=.tol[["cor"]])

})

rm(list=ls())
