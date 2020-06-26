### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: replmiss() function")

source("tolerances.r") # read in tolerances

test_that("replmiss() works correctly.", {

   var1 <- c(1:4,NA,6,NA,8:10)
   var2 <- as.numeric(1:10)

   expect_identical(replmiss(var1, 0), c(1, 2, 3, 4, 0, 6, 0, 8, 9, 10))
   expect_identical(replmiss(var1, var2), as.numeric(1:10))
   expect_error(replmiss(var1, 1:9))

})
