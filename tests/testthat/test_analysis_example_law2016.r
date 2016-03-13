### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking analysis example: law2016")

### function for creating the contrast matrix X

contrmat <- function(trt1, trt2, ref) {
   all.lvls <- sort(unique(c(levels(factor(trt1)), levels(factor(trt2)))))
   trt1 <- factor(trt1, levels=all.lvls)
   trt2 <- factor(trt2, levels=all.lvls)
   X <- model.matrix(~ trt2 - 1) - model.matrix(~ trt1 - 1)
   colnames(X) <- all.lvls
   if (missing(ref))
      ref <- all.lvls[1]
   X[, colnames(X) != ref]
}

test_that("results are correct for example 1.", {

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
   X <- contrmat(EG1$ref, EG1$trt)

   ### fit model assuming consistency (tau^2_omega=0)
   modC <- rma.mv(y, S1, mods=X, intercept=FALSE, random = ~ contr | study, rho=1/2, data=EG1)
   ci <- confint(modC)

   expect_equivalent(round(modC$tau2, digits=3), 0)
   expect_equivalent(round(coef(modC), digits=3), c(-0.224, -0.167, -0.327, -0.315, -0.352, -0.649, -0.276))
   expect_equivalent(round(ci$random[1,2:3], digits=3), c(0.000, 0.071))

   ### fit Jackson's model
   modI <- rma.mv(y, S1, mods=X, intercept=FALSE, random = list(~ contr | study, ~ contr | design), rho=1/2, phi=1/2, data=EG1)
   ci <- confint(modI)

   expect_equivalent(round(modI$tau2, digits=3), 0)
   expect_equivalent(round(modI$gamma2, digits=3), 0)
   expect_equivalent(round(coef(modI), digits=3), c(-0.224, -0.167, -0.327, -0.315, -0.352, -0.649, -0.276))
   expect_equivalent(round(ci[[1]]$random[1,2:3], digits=3), c(0.000, 0.071))
   expect_equivalent(round(ci[[2]]$random[1,2:3], digits=3), c(0.000, 0.615))

})

test_that("results are correct for example 2.", {

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
   X <- contrmat(EG2$ref, EG2$trt)

   ### fit model assuming consistency (tau^2_omega=0)
   modC <- rma.mv(y, S2, mods=X, intercept=FALSE, random = ~ contr | study, rho=1/2, data=EG2)
   ci <- confint(modC)

   expect_equivalent(round(modC$tau2, digits=3), 0.548)
   expect_equivalent(round(coef(modC), digits=3), c(-1.885, -1.337, -0.740))
   expect_equivalent(round(ci$random[1,2:3], digits=3), c(0.079, 2.016))

   ### fit Jackson's model
   modI <- rma.mv(y, S2, mods=X, intercept=FALSE, random = list(~ contr | study, ~ contr | design), rho=1/2, phi=1/2, data=EG2)
   ci <- confint(modI)

   expect_equivalent(round(modI$tau2, digits=3), 0.104)
   expect_equivalent(round(modI$gamma2, digits=3), 0.539)
   expect_equivalent(round(coef(modI), digits=3), c(-1.973, -1.396, -0.657))
   expect_equivalent(round(ci[[1]]$random[1,2:3], digits=3), c(0.000, 1.666))
   expect_equivalent(round(ci[[2]]$random[1,2:3], digits=3), c(0.000, 3.960))

})
