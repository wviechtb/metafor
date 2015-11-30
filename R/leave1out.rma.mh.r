leave1out.rma.mh <- function(x, digits, transf, targs, ...) {

   if (!is.element("rma.mh", class(x)))
      stop("Argument 'x' must be an object of class \"rma.mh\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   if (!x$int.only)
      stop("Method only applicable for models without moderators.")

   if (x$k == 1)
      stop("Stopped because k = 1.")

   if (missing(digits))
      digits <- x$digits

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   #########################################################################

   b     <- rep(NA_real_, x$k.f)
   se    <- rep(NA_real_, x$k.f)
   zval  <- rep(NA_real_, x$k.f)
   pval  <- rep(NA_real_, x$k.f)
   ci.lb <- rep(NA_real_, x$k.f)
   ci.ub <- rep(NA_real_, x$k.f)
   QE    <- rep(NA_real_, x$k.f)
   QEp   <- rep(NA_real_, x$k.f)
   #tau2 <- rep(NA_real_, x$k.f)
   #I2   <- rep(NA_real_, x$k.f)
   #H2   <- rep(NA_real_, x$k.f)

   ### note: skipping NA cases

   for (i in seq_len(x$k.f)[x$not.na]) {

      if (is.element(x$measure, c("RR","OR","RD"))) {
         res <- try(suppressWarnings(rma.mh(ai=x$ai.f[-i], bi=x$bi.f[-i], ci=x$ci.f[-i], di=x$di.f[-i], measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct)), silent=TRUE)
      } else {
         res <- try(suppressWarnings(rma.mh(x1i=x$x1i.f[-i], x2i=x$x2i.f[-i], t1i=x$t1i.f[-i], t2i=x$t2i.f[-i], measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct)), silent=TRUE)
      }

      if (inherits(res, "try-error"))
         next

      b[i]     <- res$b
      se[i]    <- res$se
      zval[i]  <- res$zval
      pval[i]  <- res$pval
      ci.lb[i] <- res$ci.lb
      ci.ub[i] <- res$ci.ub
      QE[i]    <- res$QE
      QEp[i]   <- res$QEp
      #tau2[i] <- res$tau2
      #I2[i]   <- res$I2
      #H2[i]   <- res$H2

   }

   #########################################################################

   ### if requested, apply transformation function

   if (is.logical(transf) && transf && is.element(x$measure, c("OR","RR","IRR"))) ### if transf=TRUE, apply exp transformation to ORs, RRs, and IRRs
      transf <- exp

   if (is.function(transf)) {
      if (is.null(targs)) {
         b     <- sapply(b, transf)
         se    <- rep(NA,x$k.f)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         b     <- sapply(b, transf, targs)
         se    <- rep(NA,x$k.f)
         ci.lb <- sapply(ci.lb, transf, targs)
         ci.ub <- sapply(ci.ub, transf, targs)
      }
      transf <- TRUE
   }

   ### make sure order of intervals is always increasing

   tmp <- .psort(ci.lb, ci.ub)
   ci.lb <- tmp[,1]
   ci.ub <- tmp[,2]

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(estimate=b[x$not.na], se=se[x$not.na], zval=zval[x$not.na], pval=pval[x$not.na], ci.lb=ci.lb[x$not.na], ci.ub=ci.ub[x$not.na], Q=QE[x$not.na], Qp=QEp[x$not.na])
      out$slab <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(estimate=b, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, Q=QE, Qp=QEp)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop("Missing values in results.")

   #out <- out[-c(2,3,4)]
   out$digits <- digits
   out$transf <- transf

   class(out) <- "list.rma"
   return(out)

}
