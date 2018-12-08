### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking analysis example: ishak2007")

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
                 struct = "DIAG", data = dat.long)

   ### Table 1, column "Time-specific (Independence)"
   expect_equivalent(round(coef(res),1), c(-24.9, -27.5, -28.5, -24.1))
   expect_equivalent(round(res$tau2,1), c(23.1, 27.8, 27.7, 29.9))

})

test_that("results are correct for diag(V) and random study effects.", {

   res <- rma.mv(yi, diag(V), mods = ~ factor(time) - 1, random = ~ 1 | study, data = dat.long)

   ### Table 1, column "Random study effects"
   expect_equivalent(round(coef(res),1), c(-26.2, -27.2, -28.5, -25.6))
   expect_equivalent(round(res$sigma2,1), 26.7)

})

test_that("results are correct for diag(V) and struct='ID'.", {

   res <- rma.mv(yi, diag(V), mods = ~ factor(time) - 1, random = ~ factor(time) | study,
                 struct = "ID", data = dat.long)

   ### not in paper
   expect_equivalent(round(coef(res),1), c(-24.9, -27.5, -28.5, -24.2))
   expect_equivalent(round(res$tau2,1), 26.7)

})

test_that("results are correct for diag(V) and struct='HAR'.", {

   res <- rma.mv(yi, diag(V), mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "HAR", data = dat.long)

   ### Table 1, column "Correlated random time effects"
   expect_equivalent(round(coef(res),1), c(-26.0, -27.3, -28.6, -25.8)) # -27.5 in Table vs -27.3
   expect_equivalent(round(res$tau2,1), c(20.3, 36.0, 26.4, 30.1)) # 20.4 in Table vs 20.3
   expect_equivalent(round(res$rho,2), 1.00)

})

test_that("results are correct for struct='HAR'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "HAR", data = dat.long)

   ### Table 1, column "Multivariate model"
   expect_equivalent(round(coef(res),1), c(-25.9, -27.5, -28.7, -26.5))
   expect_equivalent(round(res$tau2,1), c(22.7, 33.7, 26.1, 31.2)) # 22.6 in Table vs 22.7; 31.1 in Table vs 31.2
   expect_equivalent(round(res$rho,2), 0.88)

})

test_that("results are correct for struct='AR'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "AR", data = dat.long)

   ### not in paper
   expect_equivalent(round(coef(res),1), c(-25.9, -27.4, -28.7, -26.4))
   expect_equivalent(round(res$tau2,1), 26.7)
   expect_equivalent(round(res$rho,2), 0.87)

})

test_that("results are correct for struct='HCS'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ factor(time) | study,
                 struct = "HCS", data = dat.long)

   ### not in paper
   expect_equivalent(round(coef(res),1), c(-25.9, -27.3, -28.7, -26.7))
   expect_equivalent(round(res$tau2,1), c(20.9, 32.7, 27.7, 32.2))

})

test_that("results are correct for struct='CAR'.", {

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "CAR", data = dat.long)

   ### not in paper
   expect_equivalent(round(coef(res),1), c(-25.9, -27.4, -28.7, -26.4))
   expect_equivalent(round(res$tau2,1), 26.7)
   expect_equivalent(round(res$rho,2), 0.87)

})

test_that("results are correct for struct='CAR' with unequally spaced time points.", {

   dat.long$time[dat.long$time == 4] <- 24/3
   dat.long$time[dat.long$time == 3] <- 12/3
   dat.long$time[dat.long$time == 2] <-  6/3
   dat.long$time[dat.long$time == 1] <-  3/3

   res <- rma.mv(yi, V, mods = ~ factor(time) - 1, random = ~ time | study,
                 struct = "CAR", data = dat.long)

   ### not in paper
   expect_equivalent(round(coef(res),1), c(-26.0, -27.4, -28.7, -26.1))
   expect_equivalent(round(res$tau2,1), 27.0)
   expect_equivalent(round(res$rho,2), 0.92)

})
