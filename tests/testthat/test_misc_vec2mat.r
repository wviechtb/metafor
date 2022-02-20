### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: vec2mat() function")

source("settings.r")

test_that("vec2mat() works correctly.", {

   sav <- vec2mat(1:6, corr=FALSE)
   expect_identical(sav, structure(c(NA, 1L, 2L, 3L, 1L, NA, 4L, 5L, 2L, 4L, NA, 6L, 3L, 5L, 6L, NA), .Dim = c(4L, 4L)))

   sav <- vec2mat(round(seq(0.2, 0.7, by=0.1), 1), corr=TRUE)
   expect_identical(sav, structure(c(1, 0.2, 0.3, 0.4, 0.2, 1, 0.5, 0.6, 0.3, 0.5, 1, 0.7, 0.4, 0.6, 0.7, 1), .Dim = c(4L, 4L)))

   sav <- vec2mat(1:10, diag=TRUE)
   expect_identical(sav, structure(c(1L, 2L, 3L, 4L, 2L, 5L, 6L, 7L, 3L, 6L, 8L, 9L, 4L, 7L, 9L, 10L), .Dim = c(4L, 4L)))

   sav <- vec2mat(1:6, corr=FALSE, dimnames=c("A","B","C","D"))
   expect_identical(sav, structure(c(NA, 1L, 2L, 3L, 1L, NA, 4L, 5L, 2L, 4L, NA, 6L, 3L, 5L, 6L, NA), .Dim = c(4L, 4L), .Dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))))

})

rm(list=ls())
