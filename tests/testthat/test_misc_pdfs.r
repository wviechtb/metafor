### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: pdfs of various measures")

source("tolerances.r") # read in tolerances

test_that(".dsmd() works correctly.", {

   d <- metafor:::.dsmd(0.5, n1=15, n2=15, theta=0.8, correct=TRUE)
   expect_equivalent(d, 0.8208, tolerance=.tol[["den"]])

   d <- metafor:::.dsmd(0.5, n1=15, n2=15, theta=0.8, correct=FALSE)
   expect_equivalent(d, 0.7757, tolerance=.tol[["den"]])

})

test_that(".dcor() works correctly.", {

   d <- metafor:::.dcor(0.5, n=15, rho=0.8)
   expect_equivalent(d, 0.2255, tolerance=.tol[["den"]])

})

test_that(".dzcor() works correctly.", {

   d <- metafor:::.dzcor(0.5, n=15, rho=0.8)
   expect_equivalent(d, 0.1183, tolerance=.tol[["den"]])

   d <- metafor:::.dzcor(0.5, n=15, zrho=transf.rtoz(0.8))
   expect_equivalent(d, 0.1183, tolerance=.tol[["den"]])

})
