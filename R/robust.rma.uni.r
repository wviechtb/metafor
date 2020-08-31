robust.rma.uni <- function(x, cluster, adjust=TRUE, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.uni"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.uni\"."))

   if (missing(cluster))
      stop(mstyle$stop("Need to specify 'cluster' variable."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   level <- ifelse(x$level == 0, 1, ifelse(x$level >= 1, (100-x$level)/100, ifelse(x$level > .5, 1-x$level, x$level)))

   #########################################################################

   ### process cluster variable
   ### note: cluster variable is assumed to be of the same length as the original data passed to the model fitting function
   ###       so we have to apply the same subsetting (if necessary) and removing of missings as done during model fitting

   if (!is.null(x$subset))
      cluster <- cluster[x$subset]

   cluster <- cluster[x$not.na]

   ### checks on cluster variable

   if (anyNA(cluster))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster) != x$k)
      stop(mstyle$stop(paste0("Length of variable specified via 'cluster' (", length(cluster), ") does not match length of data (", x$k, ").")))

   ### number of clusters

   n <- length(unique(cluster))

   ### compute degrees of freedom
   ### note: Stata with vce(robust) also uses n-p as the dfs, but with vce(cluster <var>) always uses n-1 (which seems inconsistent)

   dfs <- n - x$p

   ### check if dfs are positive (note: this also handles the case where there is a single cluster)

   if (dfs <= 0)
      stop(mstyle$stop(paste0("Number of clusters (", n, ") must be larger than the number of fixed effects (", x$p, ").")))

   ### note: since we use split() below and then put things back together into a block-diagonal matrix,
   ### we have to make sure everything is properly ordered by the cluster variable; otherwise, the 'meat'
   ### block-diagonal matrix is not in the same order as the rest; so we sort all relevant variables by
   ### the cluster variable (including the cluster variable itself)

   ocl <- order(cluster)
   cluster <- cluster[ocl]

   ### construct bread = (X'WX)^-1 X'W, where W is the weight matrix

   if (x$weighted) {

      ### for weighted analysis

      if (is.null(x$weights)) {

         ### if no weights were specified, then vb = (X'WX)^-1, so we can use that part

         wi    <- 1/(x$vi + x$tau2)
         wi    <- wi[ocl]
         W     <- diag(wi, nrow=x$k, ncol=x$k)
         bread <- x$vb %*% crossprod(x$X[ocl,], W)

      } else {

         ### if weights were specified, then vb cannot be used

         A     <- diag(x$weights[ocl], nrow=x$k, ncol=x$k)
         stXAX <- .invcalc(X=x$X[ocl,], W=A, k=x$k)
         bread <- stXAX %*% crossprod(x$X[ocl,], A)

      }

   } else {

      ### for unweighted analysis

      stXX  <- .invcalc(X=x$X[ocl,], W=diag(x$k), k=x$k)
      bread <- stXX %*% t(x$X[ocl,])

   }

   ### construct meat part

   ei <- c(x$yi - x$X %*% x$beta) ### use this instead of resid(), since this guarantees that the length is correct
   ei <- ei[ocl]

   cluster <- factor(cluster, levels=unique(cluster))

   meat <- bldiag(lapply(split(ei, cluster), function(e) tcrossprod(e)))

   ### construct robust var-cov matrix

   vb <- bread %*% meat %*% t(bread)

   ### apply adjustments to vb as needed

   ### suggested in Hedges, Tipton, & Johnson (2010) -- analogous to HC1 adjustment

   if (.isTRUE(adjust))
      vb <- (n / dfs) * vb

   ### what Stata does

   if (is.character(adjust) && (adjust=="Stata" || adjust=="Stata1"))
      vb <- (n / (n-1) * (x$k-1) / (x$k-x$p)) * vb ### when the model was fitted with regress

   if (is.character(adjust) && adjust=="Stata2")
      vb <- (n / (n-1)) * vb                       ### when the model was fitted with mixed

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

   QM <- QM / x$m ### note: m is the number of coefficients in btt, not the number of unique clusters
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
   res$s2w   <- 1 ### just in case test="knha" originally
   res$meat <- matrix(NA_real_, nrow=nrow(meat), ncol=ncol(meat))
   res$meat[ocl,ocl] <- meat

   class(res) <- c("robust.rma", "rma", "rma.uni")
   return(res)

}
