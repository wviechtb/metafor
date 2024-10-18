### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true"); .tol[1:length(.tol)] <- .0001

context("Checking misc: proper handling of missing values")

source("settings.r")

test_that("rma.glmm() handles NAs correctly.", {

   skip_on_cran()

   dat <- data.frame(ni = rep(20, 10),
                     xi =   c(NA, 4, 0, 0, 2, 2, 3, 8, 9, 2),
                     mod1 = c(0, NA, 0, 0, 0, 0, 0, 1, 1, 1),
                     mod2 = c(0,  0, 0, 1, 0, 0, 0, 0, 0, 0))

   ### 1) NA in table data for study 1
   ### 2) NA for mod1 in study 2
   ### 3) if add=0, then yi/vi pair will be NA/NA for study 3
   ### 4) if add=0, then yi/vi pair will be NA/NA for study 4, which causes the X.yi matrix to be rank deficient after row 4 is removed
   ### note: even for the model fitting itself, study 4 is a problem, because the log(odds) for study 4 is -Inf, so the coefficient for
   ###       mod2 is in essence also -Inf; on x86_64-w64-mingw32/x64 (64-bit) with lme4 version 1.1-7, this just barely converges, but
   ###       may fail in other cases; so checks with both moderators included are skipped on CRAN

   expect_warning(res <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ mod1, data=dat))

   ### k, length of xi/mi, and number of rows in X must be equal to 8 (studies 1 and 2 removed due to NAs in table data)
   expect_equivalent(res$k, 8)
   expect_equivalent(length(res$outdat$xi), 8)
   expect_equivalent(length(res$outdat$mi), 8)
   expect_equivalent(nrow(res$X), 8)

   ### k.yi and length of yi/vi must be equal to 8 (studies 1 and 2 removed due to NAs in table data)
   expect_equivalent(res$k.yi, 8)
   expect_equivalent(length(res$yi), 8)
   expect_equivalent(length(res$vi), 8)

   ### full data saved in .f elements
   expect_equivalent(res$k.f, 10)
   expect_equivalent(length(res$outdat.f$xi), 10)
   expect_equivalent(length(res$outdat.f$mi), 10)
   expect_equivalent(nrow(res$X.f), 10)
   expect_equivalent(length(res$yi.f), 10)
   expect_equivalent(length(res$vi.f), 10)

   ### now use add=0, so that studies 3 and 4 have NA/NA for yi/vi

   expect_warning(res <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ mod1, data=dat, add=0))

   ### k, length of xi/mi, and number of rows in X must be equal to 8 (studies 1 and 2 removed due to NAs in table data, but studies 3 and 4 included in the model fitting)
   expect_equivalent(res$k, 8)
   expect_equivalent(length(res$outdat$xi), 8)
   expect_equivalent(length(res$outdat$mi), 8)
   expect_equivalent(nrow(res$X), 8)

   ### k.yi and length of yi/vi must be equal to 6 (studies 1 and 2 removed due to NAs in table data and studies 3 and 4 have NA/NA for yi/vi)
   expect_equivalent(res$k.yi, 6)
   expect_equivalent(length(res$yi), 6)
   expect_equivalent(length(res$vi), 6)

   ### full data saved in .f elements
   expect_equivalent(res$k.f, 10)
   expect_equivalent(length(res$outdat.f$xi), 10)
   expect_equivalent(length(res$outdat.f$mi), 10)
   expect_equivalent(nrow(res$X.f), 10)
   expect_equivalent(length(res$yi.f), 10)
   expect_equivalent(length(res$vi.f), 10)

   ### include both mod1 and mod2 in the model and use add=0, so that studies 3 and 4 have NA/NA for yi/vi
   ### as a result, the model matrix for X.yi is rank deficient, so that in essence mod2 needs to be removed for the I^2/H^2 computation
   ### also note that the coefficient for mod2 is technically -Inf (since xi=0 for the only study where mod2=1); glmer() therefore issues
   ### several warnings

   expect_warning(res <- rma.glmm(measure="PLO", xi=xi, ni=ni, mods = ~ mod1 + mod2, data=dat, add=0))

   ### k, length of xi/mi, and number of rows in X must be equal to 8 (studies 1 and 2 removed due to NAs in table data, but studies 3 and 4 included in the model fitting)
   expect_equivalent(res$k, 8)
   expect_equivalent(length(res$outdat$xi), 8)
   expect_equivalent(length(res$outdat$mi), 8)
   expect_equivalent(nrow(res$X), 8)

   ### k.yi and length of yi/vi must be equal to 6 (studies 1 and 2 removed due to NAs in table data and studies 3 and 4 have NA/NA for yi/vi)
   expect_equivalent(res$k.yi, 6)
   expect_equivalent(length(res$yi), 6)
   expect_equivalent(length(res$vi), 6)

   ### full data saved in .f elements
   expect_equivalent(res$k.f, 10)
   expect_equivalent(length(res$outdat.f$xi), 10)
   expect_equivalent(length(res$outdat.f$mi), 10)
   expect_equivalent(nrow(res$X.f), 10)
   expect_equivalent(length(res$yi.f), 10)
   expect_equivalent(length(res$vi.f), 10)

})

rm(list=ls())
