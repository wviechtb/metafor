blup.rma.uni <- function(x, level, digits, transf, targs, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma.uni", notav=c("rma.uni.selmodel", "rma.gen"))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (is.null(x$X.f) || is.null(x$yi.f))
      stop(mstyle$stop("Information needed to compute the BLUPs is not available in the model object."))

   if (missing(level))
      level <- x$level

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   level <- .level(level)

   if (is.element(x$test, c("knha","adhoc","t"))) {
      crit <- qt(level/2, df=x$ddf, lower.tail=FALSE)
   } else {
      crit <- qnorm(level/2, lower.tail=FALSE)
   }

   ### TODO: check computations for user-defined weights

   if (!is.null(x$weights) || !x$weighted)
      stop(mstyle$stop("Extraction of random effects not available for models with non-standard weights."))

   ddd <- list(...)

   .chkdots(ddd, c("code1", "code2"))

   if (!is.null(ddd[["code1"]]))
      eval(expr = parse(text = ddd[["code1"]]))

   #########################################################################

   pred  <- rep(NA_real_, x$k.f)
   vpred <- rep(NA_real_, x$k.f)

   ### see Appendix in: Raudenbush, S. W., & Bryk, A. S. (1985). Empirical
   ### Bayes meta-analysis. Journal of Educational Statistics, 10(2), 75-98

   x$tau2.f <- .expand1(x$tau2.f, x$k.f)

   li <- ifelse(is.infinite(x$tau2.f), 1, x$tau2.f / (x$tau2.f + x$vi.f))

   for (i in seq_len(x$k.f)[x$not.na]) { # note: skipping NA cases

      if (!is.null(ddd[["code2"]]))
         eval(expr = parse(text = ddd[["code2"]]))

      Xi <- matrix(x$X.f[i,], nrow=1)

      pred[i]  <- li[i] * x$yi.f[i] + (1 - li[i])   * Xi %*% x$beta

      if (li[i] == 1) {
         vpred[i] <- li[i] * x$vi.f[i]
      } else {
         vpred[i] <- li[i] * x$vi.f[i] + (1 - li[i])^2 * Xi %*% tcrossprod(x$vb,Xi)
      }

   }

   se <- sqrt(vpred)
   pi.lb <- pred - crit * se
   pi.ub <- pred + crit * se

   #########################################################################

   ### if requested, apply transformation function to 'pred' and interval bounds

   if (is.function(transf)) {
      if (is.null(targs)) {
         pred  <- sapply(pred, transf)
         se    <- rep(NA_real_, x$k.f)
         pi.lb <- sapply(pi.lb, transf)
         pi.ub <- sapply(pi.ub, transf)
      } else {
         if (!is.primitive(transf) && !is.null(targs) && length(formals(transf)) == 1L)
            stop(mstyle$stop("Function specified via 'transf' does not appear to have an argument for 'targs'."))
         pred  <- sapply(pred, transf, targs)
         se    <- rep(NA_real_, x$k.f)
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
      stop(mstyle$stop("Missing values in results."))

   #########################################################################

   out$digits <- digits
   out$transf <- transf

   class(out) <- "list.rma"
   return(out)

}
