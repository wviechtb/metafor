robust.rma.mv <- function(x, cluster, adjust=TRUE, clubSandwich=FALSE, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.mv")

   if (missing(cluster))
      stop(mstyle$stop("Must specify 'cluster' variable."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   level <- .level(x$level)

   ddd <- list(...)

   .chkdots(ddd, c("vcov", "coef_test", "conf_test", "wald_test", "verbose"))

   #########################################################################

   ### process cluster variable
   ### note: cluster variable is assumed to be of the same length as the original data passed to the model fitting function
   ###       so we have to apply the same subsetting (if necessary) and removing of missings as done during model fitting

   mf <- match.call()
   cluster <- .getx("cluster", mf=mf, data=x$data)

   if (length(cluster) != x$k.all)
      stop(mstyle$stop(paste0("Length of variable specified via 'cluster' (", length(cluster), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   if (!is.null(x$subset))
      cluster <- cluster[x$subset]

   cluster <- cluster[x$not.na]

   if (anyNA(cluster))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster) == 0L)
      stop(mstyle$stop(paste0("Cannot find 'cluster' variable (or it has zero length).")))

   ### number of clusters

   n <- length(unique(cluster))

   ### compute degrees of freedom
   ### note: Stata with vce(robust) also uses n-p as the dfs, but with vce(cluster <var>) always uses n-1 (which seems inconsistent)

   dfs <- n - x$p

   ### check if dfs are positive (note: this also handles the case where there is a single cluster)

   if (dfs <= 0)
      stop(mstyle$stop(paste0("Number of clusters (", n, ") must be larger than the number of fixed effects (", x$p, ").")))

   ### use clubSandwich if requested to do so

   if (clubSandwich) {

      if (!suppressMessages(requireNamespace("clubSandwich", quietly=TRUE)))
         stop(mstyle$stop(paste0("Please install the 'clubSandwich' package to make use of its methods.")))

      ### check for vcov, coef_test, conf_test, and wald_test arguments in ... and set values accordingly

      if (is.null(ddd$vcov)) {
         ddd$vcov <- "CR2"
      } else {
         ddd$vcov <- match.arg(ddd$vcov, c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3"))
      }
      if (is.null(ddd$coef_test)) {
         ddd$coef_test <- "Satterthwaite"
      } else {
         ddd$coef_test <- match.arg(ddd$coef_test, c("z", "naive-t", "Satterthwaite", "saddlepoint"))
      }
      if (is.null(ddd$conf_test)) {
         ddd$conf_test <- ddd$coef_test
         if (ddd$conf_test == "saddlepoint") {
            ddd$conf_test <- "Satterthwaite"
            warning(mstyle$warning("Cannot use 'saddlepoint' for conf_test() - using 'Satterthwaite' instead."), call.=FALSE)
         }
      } else {
         ddd$conf_test <- match.arg(ddd$conf_test, c("z", "naive-t", "Satterthwaite"))
      }
      if (is.null(ddd$wald_test)) {
         ddd$wald_test <- "HTZ"
      } else {
         ddd$wald_test <- match.arg(ddd$conf_test, c("chi-sq", "Naive-F", "HTA", "HTB", "HTZ", "EDF", "EDT"))
      }

      ### calculate cluster-robust var-cov matrix of the estimated fixed effects

      vb <- try(clubSandwich::vcovCR(x, cluster=cluster, type=ddd$vcov), silent=!isTRUE(ddd$verbose))

      if (inherits(vb, "try-error"))
         stop(mstyle$stop(paste0("Could not obtain the cluster-robust variance-covariance matrix (use verbose=TRUE for more details).")))

      ### obtain cluster-robust inferences

      cs.coef <- try(clubSandwich::coef_test(x, cluster=cluster, vcov=vb, test=ddd$coef_test), silent=!isTRUE(ddd$verbose))

      if (inherits(cs.coef, "try-error"))
         stop(mstyle$stop(paste0("Could not obtain the cluster-robust tests (use verbose=TRUE for more details).")))

      cs.conf <- try(clubSandwich::conf_int(x,  cluster=cluster, vcov=vb, test=ddd$conf_test, level=1-level), silent=!isTRUE(ddd$verbose))

      if (inherits(cs.conf, "try-error"))
         stop(mstyle$stop(paste0("Could not obtain the cluster-robust confidence intervals (use verbose=TRUE for more details).")))

      if (x$int.only) {

         cs.wald <- NA

      } else {

         cs.wald <- try(clubSandwich::Wald_test(x, cluster=cluster, vcov=vb, test=ddd$wald_test, constraints=clubSandwich::constrain_zero(x$btt)), silent=!isTRUE(ddd$verbose))

         if (inherits(cs.wald, "try-error"))
            stop(mstyle$stop(paste0("Could not obtain the cluster-robust Wald test (use verbose=TRUE for more details).")))

      }

      #return(list(coef_test=cs.coef, conf_int=cs.conf, Wald_test=cs.wald))

      vbest <- ddd$vcov

      beta  <- x$beta
      se    <- cs.coef$SE
      zval  <- cs.coef$tstat
      pval  <- switch(ddd$coef_test, "z" = cs.coef$p_z, "naive-t" = cs.coef$p_t, "Satterthwaite" = cs.coef$p_Satt, "saddlepoint" = cs.coef$p_saddle)
      dfs   <- switch(ddd$coef_test, "z" = Inf, "naive-t" = n-1, "Satterthwaite" = cs.coef$df, "saddlepoint" = NA)
      ci.lb <- cs.conf$CI_L # note: if ddd$coef_test != ddd$conf_test, dfs for CI may be different
      ci.ub <- cs.conf$CI_U

      if (x$int.only) {
         QM   <- max(0, zval^2)
         QMdf <- c(1, dfs)
         QMp  <- pval
      } else {
         QM   <- max(0, cs.wald$Fstat)
         QMdf <- c(cs.wald$df_num, cs.wald$df_denom)
         QMp  <- cs.wald$p_val
      }

      x$sandwiches <- list(coef_test=cs.coef, conf_int=cs.conf, Wald_test=cs.wald)
      x$coef_test <- ddd$coef_test
      x$conf_test <- ddd$conf_test
      x$wald_test <- ddd$wald_test

   } else {

      ### note: since we use split() below and then put things back together into a block-diagonal matrix,
      ### we have to make sure everything is properly ordered by the cluster variable; otherwise, the 'meat'
      ### block-diagonal matrix is not in the same order as the rest; so we sort all relevant variables by
      ### the cluster variable (including the cluster variable itself)

      ocl <- order(cluster)
      cluster <- cluster[ocl]

      ### construct bread = (X'WX)^-1 X'W, where W is the weight matrix

      if (is.null(x$W)) {

         ### if no weights were specified, then vb = (X'WX)^-1, so we can use that part

         W <- try(chol2inv(chol(x$M[ocl,ocl])), silent=TRUE)

         if (inherits(W, "try-error"))
            stop(mstyle$stop("Cannot invert marginal var-cov matrix."))

         bread <- x$vb %*% crossprod(x$X[ocl,], W)

      } else {

         ### if weights were specified, then vb cannot be used

         A     <- x$W[ocl,ocl]
         stXAX <- chol2inv(chol(as.matrix(t(x$X[ocl,]) %*% A %*% x$X[ocl,]))) ### as.matrix() to avoid some issues with the matrix being not symmetric (when it must be)
         bread <- stXAX %*% crossprod(x$X[ocl,], A)

      }

      ### construct meat part

      ei <- c(x$yi - x$X %*% x$beta) ### use this instead of resid(), since this guarantees that the length is correct
      ei <- ei[ocl]

      cluster <- factor(cluster, levels=unique(cluster))

      if (x$sparse) {
         meat <- bdiag(lapply(split(ei, cluster), function(e) tcrossprod(e)))
      } else {
         meat <- bldiag(lapply(split(ei, cluster), function(e) tcrossprod(e)))
      }

      ### construct robust var-cov matrix

      vb <- bread %*% meat %*% t(bread)

      ### apply adjustments to vb as needed

      vbest <- "CR0"

      ### suggested in Hedges, Tipton, & Johnson (2010) -- analogous to HC1 adjustment

      if (.isTRUE(adjust)) {
         vb <- (n / dfs) * vb
         vbest <- "CR1"
      }

      ### what Stata does

      if (is.character(adjust) && (adjust=="Stata" || adjust=="Stata1")) {
         vb <- (n / (n-1) * (x$k-1) / (x$k-x$p)) * vb ### when the model was fitted with regress
         vbest <- "CR1.S1"
      }

      if (is.character(adjust) && adjust=="Stata2") {
         vb <- (n / (n-1)) * vb                       ### when the model was fitted with mixed
         vbest <- "CR1.S2"
      }

      ### dim(vb) is pxp and not sparse, so this won't blow up
      ### as.matrix() helps to avoid some issues with 'vb' appearing as non-symmetric (when it must be)

      if (x$sparse)
         vb <- as.matrix(vb)

      ### prepare results

      beta <- x$beta
      se <- sqrt(diag(vb))
      names(se) <- NULL
      zval <- c(beta/se)
      pval <- 2*pt(abs(zval), df=dfs, lower.tail=FALSE)
      crit <- qt(level/2, df=dfs, lower.tail=FALSE)
      ci.lb <- c(beta - crit * se)
      ci.ub <- c(beta + crit * se)

      QM <- try(as.vector(t(beta)[x$btt] %*% chol2inv(chol(vb[x$btt,x$btt])) %*% beta[x$btt]), silent=TRUE)

      if (inherits(QM, "try-error"))
         QM <- NA

      QM   <- QM / x$m ### note: m is the number of coefficients in btt, not the number of clusters
      QMdf <- c(x$m, dfs)
      QMp  <- pf(QM, df1=x$m, df2=dfs, lower.tail=FALSE)

      ### don't need this anymore at the moment

      #x$meat <- matrix(NA_real_, nrow=nrow(meat), ncol=ncol(meat))
      #x$meat[ocl,ocl] <- as.matrix(meat)

   }

   #########################################################################

   ### table of cluster variable
   tcl <- table(cluster)

   x$digits <- digits

   ### replace elements with robust results

   x$ddf   <- dfs
   x$dfs   <- dfs
   x$vb    <- vb
   x$se    <- se
   x$zval  <- zval
   x$pval  <- pval
   x$ci.lb <- ci.lb
   x$ci.ub <- ci.ub
   x$QM    <- QM
   x$QMdf  <- QMdf
   x$QMp   <- QMp
   x$n     <- n
   x$tcl   <- tcl
   x$test  <- "t"
   x$vbest <- vbest
   x$robumethod <- ifelse(clubSandwich, "clubsw", "default")

   class(x) <- c("robust.rma", "rma", "rma.mv")
   return(x)

}
