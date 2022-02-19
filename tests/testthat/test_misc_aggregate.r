### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: aggregate() function")

source("tolerances.r") # read in tolerances

test_that("aggregate() works correctly for 'dat.konstantopoulos2011'.", {

   dat <- dat.konstantopoulos2011
   agg <- aggregate(dat, cluster=district, struct="ID")

   expect_equivalent(c(agg$yi), c(-0.125687, 0.06654, 0.350303, 0.499691, 0.051008, -0.041842, 0.885529, -0.02875, 0.250475, 0.015033, 0.161917), tolerance=.tol[["est"]])
   expect_equivalent(c(agg$vi), c(0.032427, 0.003981, 0.006664, 0.001443, 0.001549, 0.000962, 0.003882, 0.000125, 0.001799, 0.006078, 0.018678), tolerance=.tol[["var"]])


})

test_that("aggregate() works correctly for 'dat.assink2016'.", {

   dat <- dat.assink2016
   dat <- escalc(yi=yi, vi=vi, data=dat)

   agg <- aggregate(dat, cluster=study, rho=0.6)

   expect_equivalent(c(agg$yi), c(0.162877, 0.406036, 1.079003, -0.0447, 1.549, -0.054978, 1.007244, 0.3695, 0.137862, 0.116737, 0.525765, 0.280461, 0.301829, 0.035593, 0.090821, 0.018099, -0.055203), tolerance=.tol[["est"]])
   expect_equivalent(c(agg$vi), c(0.019697, 0.005572, 0.083174, 0.0331, 0.1384, 0.02139, 0.054485, 0.0199, 0.027057, 0.010729, 0.011432, 0.002814, 0.011, 0.001435, 0.126887, 0.016863, 0.007215), tolerance=.tol[["var"]])

})

test_that("aggregate() works correctly for 'dat.ishak2007'.", {

   dat <- dat.ishak2007
   dat <- reshape(dat.ishak2007, direction="long", idvar="study", v.names=c("yi","vi"),
                       varying=list(c(2,4,6,8), c(3,5,7,9)))
   dat <- dat[order(study, time),]
   dat <- dat[!is.na(yi),]
   rownames(dat) <- NULL

   agg <- aggregate(dat, cluster=study, struct="CAR", time=time, phi=0.9)

   expect_equivalent(c(agg$yi), c(-33.4, -28.137183, -21.1, -17.22908, -32.9, -26.342019, -31.37934, -25, -36, -21.275427, -8.6, -28.830656, -28.00566, -35.277625, -28.02381, -24.818713, -36.3, -29.4, -33.552998, -20.6, -33.9, -35.4, -34.9, -32.7, -26.471326, -32.753685, -18.412199, -29.2, -31.7, -32.46738, -31.7, -35.274832, -30.189494, -17.6, -22.9, -36, -22.5, -20.67624, -9.3, -25.52315, -16.7, -29.440786, -31.221009, -20.73355, -37.982183, -22.1), tolerance=.tol[["est"]])
   expect_equivalent(c(agg$vi), c(14.3, 5.611511, 7.3, 4.562371, 125, 4.132918, 86.117899, 17, 5, 6.308605, 41, 20.229622, 7.743863, 5.632795, 3.438095, 12.975915, 27.3, 10.7, 1.895013, 25.3, 20.1, 21.2, 18, 16.3, 29.751824, 9.417499, 5.156788, 5.8, 12.4, 24.954806, 19.1, 17.528303, 8.508767, 28.4, 20, 27.7, 20.3, 1.379225, 85.2, 15.281948, 9.8, 179.802277, 3.317364, 15.082821, 20.888464, 40.8), tolerance=.tol[["var"]])

})
