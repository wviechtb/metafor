skip_if_not_installed("clubSandwich")
skip_if_not_installed("microbenchmark")

# concoct fake dataset
set.seed(20240815)
studies <- 50
vi <- 4 / (20 + rpois(studies, lambda = 30))
n_es <- sample(1:50, size = studies, replace = TRUE)
dat <- data.frame(
  study = rep(1:studies, n_es),
  esid = 1:sum(n_es),
  vi = rep(vi, n_es)
)

# compare timings
timings <- microbenchmark::microbenchmark(
  "metafor-old" = vcalc(vi, cluster = study, obs = esid, data = dat, rho = 0.6, sparse = TRUE, new = FALSE),
  "metafor-new" = vcalc(vi, cluster = study, obs = esid, data = dat, rho = 0.6, sparse = TRUE, new = TRUE),
  "clubSandwich" = clubSandwich::impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, return_list = TRUE),
  times = 10
)

summary(timings)

test_that("vcalc(new = TRUE) returns results equivalent to new = FALSE.", {
  V_old <- vcalc(vi, cluster = study, obs = esid, data = dat, rho = 0.6, sparse = TRUE, new = FALSE)
  V_new <- vcalc(vi, cluster = study, obs = esid, data = dat, rho = 0.6, sparse = TRUE, new = TRUE)
  
  expect_equal(V_old, V_new)
})
