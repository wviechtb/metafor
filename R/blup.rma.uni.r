blup.rma.uni <- function(x, level, digits, transf, targs, ...) {

   if (!inherits(x, "rma.uni"))
      stop("Argument 'x' must be an object of class \"rma.uni\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (missing(level))
      level <- x$level

   if (missing(digits))
      digits <- x$digits

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   level <- ifelse(level > 1, (100-level)/100, ifelse(level > .5, 1-level, level))

   if (is.element(x$test, c("knha","adhoc","t"))) {
      crit <- qt(level/2, df=x$dfs, lower.tail=FALSE)
   } else {
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   #########################################################################

   pred  <- rep(NA_real_, x$k.f)
   vpred <- rep(NA_real_, x$k.f)

   ### see Appendix in: Raudenbush, S. W., & Bryk, A. S. (1985). Empirical
   ### Bayes meta-analysis. Journal of Educational Statistics, 10(2), 75-98

   ### TODO: is this still correct if user has specified arbitrary weights?

   li <- x$tau2 / (x$tau2 + x$vi.f)

   for (i in seq_len(x$k.f)[x$not.na]) { ### note: skipping NA cases
      Xi <- matrix(x$X.f[i,], nrow=1)
      pred[i]  <- li[i] * x$yi.f[i] + (1 - li[i])   * Xi %*% x$beta
      vpred[i] <- li[i] * x$vi.f[i] + (1 - li[i])^2 * Xi %*% tcrossprod(x$vb,Xi)
   }

   se <- sqrt(vpred)
   pi.lb <- pred - crit * se
   pi.ub <- pred + crit * se

   #########################################################################

   ### if requested, apply transformation function to 'pred' and interval bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         pred  <- sapply(pred, transf)
         se    <- rep(NA,x$k.f)
         pi.lb <- sapply(pi.lb, transf)
         pi.ub <- sapply(pi.ub, transf)
      } else {
         pred  <- sapply(pred, transf, targs)
         se    <- rep(NA,x$k.f)
         pi.lb <- sapply(pi.lb, transf, targs)
         pi.ub <- sapply(pi.ub, transf, targs)
      }
      transf <- TRUE
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(pi.lb, pi.ub)
   pi.lb <- tmp[,1]
   pi.ub <- tmp[,2]

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(pred=pred[x$not.na], se=se[x$not.na], pi.lb=pi.lb[x$not.na], pi.ub=pi.ub[x$not.na])
      out$slab <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(pred=pred, se=se, pi.lb=pi.lb, pi.ub=pi.ub)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop("Missing values in results.")

   #########################################################################

   out$digits <- digits
   out$transf <- transf

   class(out) <- "list.rma"
   return(out)

}
