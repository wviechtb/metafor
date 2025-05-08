############################################################################

.tol <- c(est   = .01, # effect size estimates
          coef  = .01, # model coefficients
          pred  = .01, # predicted values, BLUPs, also residuals
          se    = .01, # standard errors
          test  = .01, # test statistics, standardized residuals
          pval  = .01, # p-values
          ci    = .01, # confidence/prediction interval bounds, CI for effects
          var   = .01, # variance components (and CIs thereof), also if sqrt(), var-cov matrices, sampling variances
          cor   = .01, # correlations, ICCs
          cov   = .01, # covariances
          sevar = .01, # SEs of variance components
          fit   = .01, # fit statistics
          r2    = .01, # R^2 type values, model importances
          het   = .01, # heterogeneity statistics (and CIs thereof)
          inf   = .01, # influence statistics, hat values
          den   = .01, # density
          misc  = .01, # miscellaneous, mix of values
          count = 0)   # count

.tol[1:length(.tol)] <- .01
.tol[1:length(.tol)] <- .01

############################################################################

.sparse <- FALSE
#.sparse <- TRUE

############################################################################

.vistest <- function(file1, file2) {
   if (isFALSE(as.logical(Sys.getenv("RUN_VIS_TESTS", "false")))) {
      return(TRUE)
   } else {
      hash1 <- suppressWarnings(system2("md5sum", file1, stdout=TRUE, stderr=TRUE))
      hash2 <- suppressWarnings(system2("md5sum", file2, stdout=TRUE, stderr=TRUE))
      if (isTRUE(attributes(hash1)$status == 1) || isTRUE(attributes(hash2)$status == 1))
         return(FALSE)
      hash1 <- strsplit(hash1, "  ")[[1]][1]
      hash2 <- strsplit(hash2, "  ")[[1]][1]
      return(identical(hash1,hash2))
      #file1 <- readLines(file1, warn=FALSE)
      #file2 <- readLines(file2, warn=FALSE)
      #file1 <- file1[!grepl("CreationDate", file1, fixed=TRUE, useBytes=TRUE)]
      #file2 <- file2[!grepl("CreationDate", file2, fixed=TRUE, useBytes=TRUE)]
      #file1 <- file1[!grepl("ModDate",      file1, fixed=TRUE, useBytes=TRUE)]
      #file2 <- file2[!grepl("ModDate",      file2, fixed=TRUE, useBytes=TRUE)]
      #file1 <- file1[!grepl("Producer",     file1, fixed=TRUE, useBytes=TRUE)]
      #file2 <- file2[!grepl("Producer",     file2, fixed=TRUE, useBytes=TRUE)]
      #return(identical(file1,file2))
   }
}

############################################################################

setmfopt(theme="default")

############################################################################
