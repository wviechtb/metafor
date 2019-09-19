### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking analysis example: law2016")

source("tolerances.r") # read in tolerances

test_that("results are correct for example 1.", {

   skip_on_cran()

   ### example 1

   EG1 <- read.table(header=TRUE, as.is=TRUE, text="
    study           y ref trt contr design
        1 -0.16561092   C   D    CD     CD
        2 -0.13597406   C   D    CD     CD
        3 -0.08012604   C   E    CE     CE
        4 -0.14746890   C   F    CF     CF
        5  0.09316853   E   F    EF     EF
        6 -0.15859403   E   F    EF     EF
        7 -0.22314355   E   F    EF     EF
        8 -0.06744128   F   G    FG     FG
        9 -0.11888254   C   H    CH     CH
       10 -0.06899287   C   H    CH     CH
       11  0.26917860   B   C    BC     BC
       12 -0.33160986   A   B    AB     AB
       13 -0.26236426   A   B    AB     AB
       14 -0.39319502   F   G    FG     FG
       15 -0.11557703   A   B    AB     AB
       16  0.00000000   E   F    EF     EF
       17 -0.40987456   A   E    AE     AE
   ")

   S1 <- structure(c(0.0294183340466069, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0.147112449467866, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0780588660166125, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.140361934247383, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0479709251030665,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0506583523716436,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.235695187165775,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.04499494438827,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.17968120987923,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.735714285714286,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.184889643463497,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0294022652280727,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.232478632478632,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.857874134296899,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0219285638496459,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.168131868131868,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0826973577700322
   ), .Dim = c(17, 17))

   ### create contrast matrix
   X <- contrmat(EG1, grp1="trt", grp2="ref", append=FALSE, last=NA)[,-1] # remove 'A' to make it the reference level

   ### fit model assuming consistency (tau^2_omega=0)
   modC <- rma.mv(y, S1, mods=X, intercept=FALSE, random = ~ contr | study, rho=1/2, data=EG1)
   ci <- confint(modC)

   expect_equivalent(modC$tau2, 0.0000, tolerance=.tol[["var"]])
   expect_equivalent(coef(modC), c(-0.2243, -0.1667, -0.3274, -0.3152, -0.3520, -0.6489, -0.2758), tolerance=.tol[["coef"]])
   expect_equivalent(ci$random[1,2:3], c(0.0000, 0.0708), tolerance=.tol[["var"]])

   ### fit Jackson's model
   modI <- rma.mv(y, S1, mods=X, intercept=FALSE, random = list(~ contr | study, ~ contr | design), rho=1/2, phi=1/2, data=EG1)
   ci <- confint(modI)

   expect_equivalent(modI$tau2, 0.0000, tolerance=.tol[["var"]])
   expect_equivalent(modI$gamma2, 0.0000, tolerance=.tol[["var"]])
   expect_equivalent(coef(modI), c(-0.2243, -0.1667, -0.3274, -0.3152, -0.3520, -0.6489, -0.2758), tolerance=.tol[["coef"]])
   expect_equivalent(ci[[1]]$random[1,2:3], c(0.0000, 0.0708), tolerance=.tol[["var"]])
   expect_equivalent(ci[[2]]$random[1,2:3], c(0.0000, 0.6153), tolerance=.tol[["var"]])

   sav <- predict(modI, newmods=c(1,0,0,0,0,0,0), transf=exp)
   sav <- c(sav[[1]], sav[[3]], sav[[4]], sav[[5]], sav[[6]])
   expect_equivalent(sav, c(0.7991, 0.6477, 0.9859, 0.6477, 0.9859), tolerance=.tol[["pred"]])

})

test_that("results are correct for example 2.", {

   skip_on_cran()

   ### example 2

   EG2 <- read.table(header=TRUE, as.is=TRUE, text="
    study           y ref trt contr design
        1 -3.61988658   A   B    AB     AB
        2  0.00000000   B   C    BC     BC
        3  0.19342045   B   C    BC     BC
        4  2.79320801   B   C    BC     BC
        5  0.24512246   B   C    BC     BC
        6  0.03748309   B   C    BC     BC
        7  0.86020127   B   D    BD     BD
        8  0.14310084   B   D    BD     BD
        9  0.07598591   C   D    CD     CD
       10 -0.99039870   C   D    CD     CD
       11 -1.74085310   A   B    AB    ABD
       11  0.34830670   A   D    AD    ABD
       12  0.40546510   B   C    BC    BCD
       12  1.91692260   B   D    BD    BCD
       13 -0.32850410   B   C    BC    BCD
       13  1.07329450   B   D    BD    BCD
   ")

   S2 <- structure(c(0.9672619, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0.24987648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.61904762,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.27958937, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.23845689, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.04321419, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.47692308, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.18416468, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.61978022, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.12650164, 0.07397504, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.07397504, 0.1583906, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.389881, 0.2857143,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2857143, 0.5151261,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4361111, 0.2111111,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2111111, 0.5380342
   ), .Dim = c(16, 16))

   ### create contrast matrix
   X <- contrmat(EG2, grp1="trt", grp2="ref", append=FALSE, last=NA)[,-1] # remove 'A' to make it the reference level

   ### fit model assuming consistency (tau^2_omega=0)
   modC <- rma.mv(y, S2, mods=X, intercept=FALSE, random = ~ contr | study, rho=1/2, data=EG2)
   ci <- confint(modC)

   expect_equivalent(modC$tau2, 0.5482, tolerance=.tol[["var"]])
   expect_equivalent(coef(modC), c(-1.8847, -1.3366, -0.7402), tolerance=.tol[["coef"]])
   expect_equivalent(ci$random[1,2:3], c(0.0788, 2.0156), tolerance=.tol[["var"]])

   ### fit Jackson's model
   modI <- rma.mv(y, S2, mods=X, intercept=FALSE, random = list(~ contr | study, ~ contr | design), rho=1/2, phi=1/2, data=EG2)
   ci <- confint(modI)

   expect_equivalent(modI$tau2, 0.1036, tolerance=.tol[["var"]])
   expect_equivalent(modI$gamma2, 0.5391, tolerance=.tol[["var"]])
   expect_equivalent(coef(modI), c(-1.9735, -1.3957, -0.6572), tolerance=.tol[["coef"]])
   expect_equivalent(ci[[1]]$random[1,2:3], c(0.0000, 1.6661), tolerance=.tol[["var"]])
   expect_equivalent(ci[[2]]$random[1,2:3], c(0.0000, 3.9602), tolerance=.tol[["var"]])

   sav <- predict(modI, newmods=c(1,0,0), transf=exp)
   sav <- c(sav[[1]], sav[[3]], sav[[4]], sav[[5]], sav[[6]])
   expect_equivalent(sav, c(0.1390, 0.0369, 0.5230, 0.0178, 1.0856), tolerance=.tol[["pred"]])

   sav <- ranef(modI)

   expect_equivalent(sav[[1]]$intrcpt, c(-0.106, -0.094, -0.078, 0.335, -0.058, -0.128, 0.026, -0.121, 0.013, -0.148, 0.029, 0.13, 0.027, 0.084, -0.101, -0.064), tolerance=.tol[["pred"]])
   expect_equivalent(sav[[1]]$se,      c(0.314, 0.293, 0.283, 0.301, 0.285, 0.282, 0.286, 0.297, 0.297, 0.304, 0.313, 0.315, 0.301, 0.304, 0.302, 0.305), tolerance=.tol[["se"]])
   expect_equivalent(sav[[2]]$intrcpt, c(-0.551, 0.152, 0.675, -0.119, -0.383, 0.104, -0.493, -0.699), tolerance=.tol[["pred"]])
   expect_equivalent(sav[[2]]$se,      c(0.64, 0.619, 0.642, 0.518, 0.543, 0.53, 0.486, 0.54), tolerance=.tol[["se"]])

   out <- capture.output(print(sav))

   sav <- predict(modI)
   expect_equivalent(sav$cr.lb, c(-4.029, -1.2853, -1.2853, -1.2853, -1.2853, -1.2853, -0.4911, -0.4911, -1.137, -1.137, -4.029, -2.7699, -1.2853, -0.4911, -1.2853, -0.4911), tolerance=.tol[["pred"]])

})
