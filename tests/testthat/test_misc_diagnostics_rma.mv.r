### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true")

context("Checking misc: model diagnostic functions for rma.mv()")

dat1 <- dat.konstantopoulos2011

dat1 <- dat1[dat1$district %in% c(11, 12, 18, 71, 108, 644),]
rownames(dat1) <- 1:nrow(dat1)
dat1$yi[dat1$district %in% 12] <- NA ### all values for district 12 are missing
dat1$yi[dat1$district %in% 18 & dat1$school == 2]  <- NA ### second value for district 18 is missing
dat1$yi[dat1$district %in% 108] <- dat1$yi[dat1$district %in% 108] + 1 ### increase district level variance
dat1$district11 <- ifelse(dat1$district == 11, 1, 0) ### dummy for district 11
dat1$study53 <- ifelse(dat1$study == 53, 1, 0) ### dummies for studies in district 644
dat1$study54 <- ifelse(dat1$study == 54, 1, 0) ### dummies for studies in district 644
dat1$study55 <- ifelse(dat1$study == 55, 1, 0) ### dummies for studies in district 644
dat1$study56 <- ifelse(dat1$study == 56, 1, 0) ### dummies for studies in district 644

set.seed(123213)
dat2 <- dat1[sample(nrow(dat1)),] ### reshuffled dataset

res1 <- suppressWarnings(rma.mv(yi, vi, mods = ~ district11 + study53 + study54 + study55 + study56, random = ~ 1 | district/school, data=dat1, slab=study))
res2 <- suppressWarnings(rma.mv(yi, vi, mods = ~ district11 + study53 + study54 + study55 + study56, random = ~ 1 | district/school, data=dat2, slab=study))

test_that("model diagnostic functions work with 'na.omit'.", {

   skip_on_cran()

   options(na.action="na.omit")

   sav1 <- rstandard(res1)
   sav2 <- rstandard(res2)
   sav2 <- sav2[match(sav1$slab, sav2$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$resid), rep(FALSE,18))

   sav1 <- rstandard(res1, cluster=dat1$district)
   sav2 <- rstandard(res2, cluster=dat2$district)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), rep(FALSE,18))
   expect_equivalent(is.na(sav1$cluster$X2), c(TRUE, rep(FALSE,4)))

   sav1 <- rstudent(res1)
   sav2 <- rstudent(res2)
   sav2 <- sav2[match(sav1$slab, sav2$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$resid), c(rep(FALSE,14), rep(TRUE,4)))

   sav1 <- rstudent(res1, cluster=dat1$district)
   sav2 <- rstudent(res2, cluster=dat2$district)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,4), rep(FALSE,10), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- rstudent(res1, cluster=dat1$district, parallel="snow")
   sav2 <- rstudent(res2, cluster=dat2$district, parallel="snow")
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,4), rep(FALSE,10), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- rstudent(res1, cluster=dat1$district, reestimate=FALSE)
   sav2 <- rstudent(res2, cluster=dat2$district, reestimate=FALSE)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,4), rep(FALSE,10), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- rstudent(res1, cluster=dat1$district, parallel="snow", reestimate=FALSE)
   sav2 <- rstudent(res2, cluster=dat2$district, parallel="snow", reestimate=FALSE)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,4), rep(FALSE,10), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1)
   sav2 <- cooks.distance(res2)
   sav2 <- sav2[match(names(sav1), names(sav2))]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(rep(FALSE,14), rep(TRUE,4)))

   sav1 <- cooks.distance(res1, cluster=dat1$district)
   sav2 <- cooks.distance(res2, cluster=dat2$district)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1, cluster=dat1$district, parallel="snow")
   sav2 <- cooks.distance(res2, cluster=dat2$district, parallel="snow")

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1, cluster=dat1$district, reestimate=FALSE)
   sav2 <- cooks.distance(res2, cluster=dat2$district, reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1, cluster=dat1$district, parallel="snow", reestimate=FALSE)
   sav2 <- cooks.distance(res2, cluster=dat2$district, parallel="snow", reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1)
   sav2 <- dfbetas(res2)
   sav2 <- sav2[match(rownames(sav1), rownames(sav2)),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(rep(FALSE,14), rep(TRUE,4)))

   sav1 <- dfbetas(res1, cluster=dat1$district)
   sav2 <- dfbetas(res2, cluster=dat2$district)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1, cluster=dat1$district, parallel="snow")
   sav2 <- dfbetas(res2, cluster=dat2$district, parallel="snow")

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1, cluster=dat1$district, reestimate=FALSE)
   sav2 <- dfbetas(res2, cluster=dat2$district, reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1, cluster=dat1$district, parallel="snow", reestimate=FALSE)
   sav2 <- dfbetas(res2, cluster=dat2$district, parallel="snow", reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(TRUE, rep(FALSE,3), TRUE))

   sav1 <- ranef(res1)
   sav2 <- ranef(res2)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$district$intrcpt), rep(FALSE,5))
   expect_equivalent(is.na(sav1$`district/school`$intrcpt), rep(FALSE,18))

})

test_that("model diagnostic functions work with 'na.pass'.", {

   skip_on_cran()

   options(na.action="na.pass")

   sav1 <- rstandard(res1)
   sav2 <- rstandard(res2)
   sav2 <- sav2[match(sav1$slab, sav2$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$resid), c(rep(FALSE,4), rep(TRUE,4), FALSE, TRUE, rep(FALSE,13)))

   sav1 <- rstandard(res1, cluster=dat1$district)
   sav2 <- rstandard(res2, cluster=dat2$district)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(FALSE,4), rep(TRUE,4), FALSE, TRUE, rep(FALSE,13)))
   expect_equivalent(is.na(sav1$cluster$X2), c(rep(TRUE,2), rep(FALSE,4)))

   sav1 <- rstudent(res1)
   sav2 <- rstudent(res2)
   sav2 <- sav2[match(sav1$slab, sav2$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$resid), c(rep(FALSE,4), rep(TRUE,4), FALSE, TRUE, rep(FALSE,9), rep(TRUE,4)))

   sav1 <- rstudent(res1, cluster=dat1$district)
   sav2 <- rstudent(res2, cluster=dat2$district)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,8), FALSE, TRUE, rep(FALSE,9), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- rstudent(res1, cluster=dat1$district, parallel="snow")
   sav2 <- rstudent(res2, cluster=dat2$district, parallel="snow")
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,8), FALSE, TRUE, rep(FALSE,9), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- rstudent(res1, cluster=dat1$district, reestimate=FALSE)
   sav2 <- rstudent(res2, cluster=dat2$district, reestimate=FALSE)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,8), FALSE, TRUE, rep(FALSE,9), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- rstudent(res1, cluster=dat1$district, parallel="snow", reestimate=FALSE)
   sav2 <- rstudent(res2, cluster=dat2$district, parallel="snow", reestimate=FALSE)
   sav2$obs <- sav2$obs[match(sav1$obs$slab, sav2$obs$slab),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$obs$resid), c(rep(TRUE,8), FALSE, TRUE, rep(FALSE,9), rep(TRUE,4)))
   expect_equivalent(is.na(sav1$cluster$X2), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1)
   sav2 <- cooks.distance(res2)
   sav2 <- sav2[match(names(sav1), names(sav2))]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(rep(FALSE,4), rep(TRUE,4), FALSE, TRUE, rep(FALSE,9), rep(TRUE,4)))

   sav1 <- cooks.distance(res1, cluster=dat1$district)
   sav2 <- cooks.distance(res2, cluster=dat2$district)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1, cluster=dat1$district, parallel="snow")
   sav2 <- cooks.distance(res2, cluster=dat2$district, parallel="snow")

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1, cluster=dat1$district, reestimate=FALSE)
   sav2 <- cooks.distance(res2, cluster=dat2$district, reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- cooks.distance(res1, cluster=dat1$district, parallel="snow", reestimate=FALSE)
   sav2 <- cooks.distance(res2, cluster=dat2$district, parallel="snow", reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1)
   sav2 <- dfbetas(res2)
   sav2 <- sav2[match(rownames(sav1), rownames(sav2)),]

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(rep(FALSE,4), rep(TRUE,4), FALSE, TRUE, rep(FALSE,9), rep(TRUE,4)))

   sav1 <- dfbetas(res1, cluster=dat1$district)
   sav2 <- dfbetas(res2, cluster=dat2$district)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1, cluster=dat1$district, parallel="snow")
   sav2 <- dfbetas(res2, cluster=dat2$district, parallel="snow")

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1, cluster=dat1$district, reestimate=FALSE)
   sav2 <- dfbetas(res2, cluster=dat2$district, reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- dfbetas(res1, cluster=dat1$district, parallel="snow", reestimate=FALSE)
   sav2 <- dfbetas(res2, cluster=dat2$district, parallel="snow", reestimate=FALSE)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$intrcpt), c(rep(TRUE,2), rep(FALSE,3), TRUE))

   sav1 <- ranef(res1)
   sav2 <- ranef(res2)

   expect_equivalent(sav1, sav2)
   expect_equivalent(is.na(sav1$district$intrcpt), c(FALSE, TRUE, rep(FALSE,4)))
   expect_equivalent(is.na(sav1$`district/school`$intrcpt), c(rep(FALSE,4), rep(TRUE,4), FALSE, TRUE, rep(FALSE,13)))

})
