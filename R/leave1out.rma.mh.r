leave1out.rma.mh <- function(x, digits, transf, targs, progbar=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.mh")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (!x$int.only)
      stop(mstyle$stop("Method only applicable for models without moderators."))

   if (x$k == 1)
      stop(mstyle$stop("Stopped because k = 1."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   ddd <- list(...)

   .chkdots(ddd, c("time"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   #########################################################################

   beta  <- rep(NA_real_, x$k.f)
   se    <- rep(NA_real_, x$k.f)
   zval  <- rep(NA_real_, x$k.f)
   pval  <- rep(NA_real_, x$k.f)
   ci.lb <- rep(NA_real_, x$k.f)
   ci.ub <- rep(NA_real_, x$k.f)
   QE    <- rep(NA_real_, x$k.f)
   QEp   <- rep(NA_real_, x$k.f)
   #tau2 <- rep(NA_real_, x$k.f)
   I2   <- rep(NA_real_, x$k.f)
   H2   <- rep(NA_real_, x$k.f)

   ### elements that need to be returned

   outlist <- "beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, QE=QE, QEp=QEp, tau2=tau2, I2=I2, H2=H2"

   ### note: skipping NA cases

   if (progbar)
      pbar <- pbapply::startpb(min=0, max=x$k.f)

   for (i in seq_len(x$k.f)) {

      if (progbar)
         pbapply::setpb(pbar, i)

      if (!x$not.na[i])
         next

      if (is.element(x$measure, c("RR","OR","RD"))) {
         args <- list(ai=x$outdat.f$ai, bi=x$outdat.f$bi, ci=x$outdat.f$ci, di=x$outdat.f$di, measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct, level=x$level, subset=-i, outlist=outlist)
      } else {
         args <- list(x1i=x$outdat.f$x1i, x2i=x$outdat.f$x2i, t1i=x$outdat.f$t1i, t2i=x$outdat.f$t2i, measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct, level=x$level, subset=-i, outlist=outlist)
      }
      res <- try(suppressWarnings(.do.call(rma.mh, args)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      beta[i]  <- res$beta
      se[i]    <- res$se
      zval[i]  <- res$zval
      pval[i]  <- res$pval
      ci.lb[i] <- res$ci.lb
      ci.ub[i] <- res$ci.ub
      QE[i]    <- res$QE
      QEp[i]   <- res$QEp
      #tau2[i] <- res$tau2
      I2[i]   <- res$I2
      H2[i]   <- res$H2

   }

   if (progbar)
      pbapply::closepb(pbar)

   #########################################################################

   ### if requested, apply transformation function

   if (.isTRUE(transf) && is.element(x$measure, c("OR","RR","IRR"))) ### if transf=TRUE, apply exp transformation to ORs, RRs, and IRRs
      transf <- exp

   if (is.function(transf)) {
      if (is.null(targs)) {
         beta  <- sapply(beta, transf)
         se    <- rep(NA,x$k.f)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         beta  <- sapply(beta, transf, targs)
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
      out <- list(estimate=beta[x$not.na], se=se[x$not.na], zval=zval[x$not.na], pval=pval[x$not.na], ci.lb=ci.lb[x$not.na], ci.ub=ci.ub[x$not.na], Q=QE[x$not.na], Qp=QEp[x$not.na], I2=I2[x$not.na], H2=H2[x$not.na])
      out$slab <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(estimate=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, Q=QE, Qp=QEp, I2=I2, H2=H2)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   out$digits <- digits
   out$transf <- transf

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(out) <- "list.rma"
   return(out)

}
