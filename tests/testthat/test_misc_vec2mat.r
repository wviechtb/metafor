### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: vec2mat() function")

source("settings.r")

test_that("vec2mat() works correctly.", {

   sav <- vec2mat(1:6, corr=FALSE)
   expect_identical(sav, structure(c(NA, 1, 2, 3, 1, NA, 4, 5, 2, 4, NA, 6, 3, 5, 6, NA), .Dim = c(4L, 4L)))

   sav <- vec2mat(round(seq(0.2, 0.7, by=0.1), 1), corr=TRUE)
   expect_identical(sav, structure(c(1, 0.2, 0.3, 0.4, 0.2, 1, 0.5, 0.6, 0.3, 0.5, 1, 0.7, 0.4, 0.6, 0.7, 1), .Dim = c(4L, 4L)))

   sav <- vec2mat(1:10, diag=TRUE)
   expect_identical(sav, structure(c(1, 2, 3, 4, 2, 5, 6, 7, 3, 6, 8, 9, 4, 7, 9, 10), .Dim = c(4L, 4L)))

   sav <- vec2mat(1:6, corr=FALSE, dimnames=c("A","B","C","D"))
   expect_identical(sav, structure(c(NA, 1, 2, 3, 1, NA, 4, 5, 2, 4, NA, 6, 3, 5, 6, NA), .Dim = c(4L, 4L), .Dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))))

})

rm(list=ls())
