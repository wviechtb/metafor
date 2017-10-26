robust.rma.mv <- function(x, cluster, adjust=TRUE, digits, ...) {

   if (!inherits(x, "rma.mv"))
      stop("Argument 'x' must be an object of class \"rma.mv\".")

   if (missing(cluster))
      stop("Need to specify 'cluster' variable.")

   if (missing(digits))
      digits <- x$digits

   level <- ifelse(x$level > 1, (100-x$level)/100, ifelse(x$level > .5, 1-x$level, x$level))

   #########################################################################

   ### process cluster variable
   ### note: cluster variable is assumed to be of the same length as the original data passed to the model fitting function
   ###       so we have to apply the same subsetting (if necessary) and removing of missings as done during model fitting

   if (!is.null(x$subset))
      cluster <- cluster[x$subset]

   cluster <- cluster[x$not.na]

   ### checks on cluster variable

   if (anyNA(cluster))
      stop("No missing values allowed in 'cluster' variable.")

   if (length(cluster) != x$k)
      stop("Length of variable specified via 'cluster' does not match length of data.")

   ### number of clusters

   n <- length(unique(cluster))

   ### compute degrees of freedom
   ### note: Stata with vce(robust) also uses n-p as the dfs, but with vce(cluster <var>) always uses n-1 (which seems inconsistent)

   dfs <- n - x$p

   ### check if dfs are positive (note: this also handles the case where there is a single cluster)

   if (dfs <= 0)
      stop(paste0("Number of clusters (", n, ") must be larger than the number of fixed effects (", x$p, ")."))

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
         stop("Cannot invert marginal var-cov matrix.")

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

   ### suggested in Hedges, Tipton, & Johnson (2010) -- analogous to HC1 adjustment

   if (is.logical(adjust) && adjust)
      vb <- (n / dfs) * vb

   ### what Stata does

   if (is.character(adjust) && (adjust=="Stata" || adjust=="Stata1"))
      vb <- (n / (n-1) * (x$k-1) / (x$k-x$p)) * vb ### when the model was fitted with regress

   if (is.character(adjust) && adjust=="Stata2")
      vb <- (n / (n-1)) * vb                       ### when the model was fitted with mixed

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

   QM <- QM / x$m ### careful: m is the number of coefficients in btt, not the number of clusters
   QMp <- pf(QM, df1=x$m, df2=dfs, lower.tail=FALSE)

   #########################################################################

   ### table of cluster variable
   tcl <- table(cluster)

   res <- x
   res$digits <- digits

   ### replace elements with robust results
   res$dfs   <- dfs
   res$vb    <- vb
   res$se    <- se
   res$zval  <- zval
   res$pval  <- pval
   res$ci.lb <- ci.lb
   res$ci.ub <- ci.ub
   res$QM    <- QM
   res$QMp   <- QMp
   res$n     <- n
   res$tcl   <- tcl
   res$test  <- "t"
   res$meat <- matrix(NA_real_, nrow=nrow(meat), ncol=ncol(meat))
   res$meat[ocl,ocl] <- as.matrix(meat)

   class(res) <- c("robust.rma", "rma", "rma.mv")
   return(res)

}
