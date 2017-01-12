robust.rma.uni <- function(x, cluster, adjust=TRUE, digits, ...) {

   if (!inherits(x, "rma.uni"))
      stop("Argument 'x' must be an object of class \"rma.uni\".")

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

   if (x$weighted) {

      ### for weighted analysis

      if (is.null(x$weights)) {

         ### if no weights were specified, then vb = (X'WX)^-1, so we can use that part

         wi <- 1/(x$vi[ocl] + x$tau2)
         W  <- diag(wi, nrow=x$k, ncol=x$k)
         bread <- x$vb %*% crossprod(x$X[ocl,], W)

      } else {

         ### if weights were specified, then vb cannot be used

         A     <- diag(x$weights[ocl], nrow=x$k, ncol=x$k)
         stXAX <- .invcalc(X=x$X[ocl,], W=A, k=x$k)
         bread <- stXAX %*% crossprod(x$X[ocl,], A)

      }

   } else {

      ### for unweighted analysis

      stXX <- .invcalc(X=x$X[ocl,], W=diag(x$k), k=x$k)
      bread <- stXX %*% t(x$X[ocl,])

   }

   ### construct meat part

   ei <- c(x$yi - x$X %*% x$b) ### use this instead of resid(), since this guarantees that the length is correct
   ei <- ei[ocl]

   cluster <- factor(cluster, levels=unique(cluster))

   meat <- bldiag(lapply(split(ei, cluster), function(e) tcrossprod(e)))

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

   ### prepare results

   b <- x$b
   se <- sqrt(diag(vb))
   names(se) <- NULL
   tval <- c(b/se)
   pval <- 2*pt(abs(tval), df=dfs, lower.tail=FALSE)
   crit <- qt(level/2, df=dfs, lower.tail=FALSE)
   ci.lb <- c(b - crit * se)
   ci.ub <- c(b + crit * se)

   QM <- try(as.vector(t(b)[x$btt] %*% chol2inv(chol(vb[x$btt,x$btt])) %*% b[x$btt]), silent=TRUE)

   if (inherits(QM, "try-error"))
      QM <- NA

   QM <- QM / x$m ### careful: m is the number of coefficients in btt, not the number of unique clusters
   QMp <- pf(QM, df1=x$m, df2=dfs, lower.tail=FALSE)

   #########################################################################

   ### table of cluster variable
   tcl <- table(cluster)

   res <- list(b=b, se=se, tval=tval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, vb=vb,
               k=x$k, k.f=x$k.f, p=x$p, m=x$m, n=n, dfs=dfs, tcl=tcl, QM=QM, QMp=QMp, yi.f=x$yi.f, vi.f=x$vi.f, X=x$X, X.f=x$X.f, method=x$method,
               int.only=x$int.only, int.incl=x$int.incl, test="t", btt=x$btt, intercept=x$intercept, digits=digits, level=x$level, tau2=x$tau2, slab=x$slab,
               slab.null=x$slab.null, not.na=x$not.na,
               fit.stats=x$fit.stats, k.eff=x$k.eff, p.eff=x$p.eff, parms=x$parms, measure=x$measure)

   class(res) <- c("robust.rma", "rma", "rma.uni")
   return(res)

}
