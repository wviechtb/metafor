### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: pdfs of various measures")

test_that(".dsmd() works correctly.", {

   d <- metafor:::.dsmd(0.5, n1=15, n2=15, theta=0.8, correct=TRUE)
   expect_equivalent(round(d,3), 0.821)

   d <- metafor:::.dsmd(0.5, n1=15, n2=15, theta=0.8, correct=FALSE)
   expect_equivalent(round(d,3), 0.776)

})

test_that(".dcor() works correctly.", {

   d <- metafor:::.dcor(0.5, n=15, rho=0.8)
   expect_equivalent(round(d,3), 0.225)

})

test_that(".dzcor() works correctly.", {

   d <- metafor:::.dzcor(0.5, n=15, rho=0.8)
   expect_equivalent(round(d,3), 0.118)

   d <- metafor:::.dzcor(0.5, n=15, zrho=transf.rtoz(0.8))
   expect_equivalent(round(d,3), 0.118)

})
