leave1out.rma.uni <- function(x, digits, transf, targs, progbar=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.uni", notav=c("robust.rma", "rma.ls", "rma.uni.selmodel"))

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
   tau2  <- rep(NA_real_, x$k.f)
   I2    <- rep(NA_real_, x$k.f)
   H2    <- rep(NA_real_, x$k.f)

   ### note: skipping NA cases
   ### also: it is possible that model fitting fails, so that generates more NAs (these NAs will always be shown in output)

   if (progbar)
      pbar <- pbapply::startpb(min=0, max=x$k.f)

   for (i in seq_len(x$k.f)) {

      if (progbar)
         pbapply::setpb(pbar, i)

      if (!x$not.na[i])
         next

      args <- list(yi=x$yi.f, vi=x$vi.f, weights=x$weights.f, intercept=TRUE, method=x$method, weighted=x$weighted, test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=-i, skipr2=TRUE)
      args <- args[!sapply(args, is.null)]
      res <- try(suppressWarnings(do.call(rma.uni, args)), silent=TRUE)

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
      tau2[i]  <- res$tau2
      I2[i]    <- res$I2
      H2[i]    <- res$H2

   }

   if (progbar)
      pbapply::closepb(pbar)

   #########################################################################

   ### if requested, apply transformation function

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
      out <- list(estimate=beta[x$not.na], se=se[x$not.na], zval=zval[x$not.na], pval=pval[x$not.na], ci.lb=ci.lb[x$not.na], ci.ub=ci.ub[x$not.na], Q=QE[x$not.na], Qp=QEp[x$not.na], tau2=tau2[x$not.na], I2=I2[x$not.na], H2=H2[x$not.na])
      out$slab <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(estimate=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, Q=QE, Qp=QEp, tau2=tau2, I2=I2, H2=H2)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   if (is.element(x$test, c("knha","adhoc","t")))
      names(out)[3] <- "tval"

   ### remove tau2 for FE/EE/CE models

   if (is.element(x$method, c("FE","EE","CE")))
      out <- out[-9]

   out$digits <- digits
   out$transf <- transf

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(out) <- "list.rma"
   return(out)

}
