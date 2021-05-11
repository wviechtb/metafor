cumul.rma.uni <- function(x, order, digits, transf, targs, progbar=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.uni", notav=c("robust.rma", "rma.ls", "rma.uni.selmodel"))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (!x$int.only)
      stop(mstyle$stop("Method only applicable for models without moderators."))

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

   .chkdots(ddd, c("time", "decreasing"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   if (is.null(ddd$decreasing)) {
      decreasing <- FALSE
   } else {
      decreasing <- ddd$decreasing
   }

   #########################################################################

   if (grepl("^order\\(", deparse(substitute(order))))
      warning(mstyle$warning("Use of order() in 'order' argument is probably erroneous."), call.=FALSE)

   if (missing(order))
      order <- seq_len(x$k.all)

   if (length(order) != x$k.all)
      stop(mstyle$stop(paste0("Length of the 'order' argument (", length(order), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   ### note: order variable is assumed to be of the same length as the size of the
   ###       original dataset passed to the model fitting function and so we apply
   ###       the same subsetting (if necessary) as was done during model fitting

   if (!is.null(x$subset))
      order <- order[x$subset]

   order <- order(order, decreasing=decreasing)

   yi.f      <- x$yi.f[order]
   vi.f      <- x$vi.f[order]
   X.f       <- cbind(x$X.f[order,])
   weights.f <- x$weights.f[order]
   not.na    <- x$not.na[order]
   slab      <- x$slab[order]
   ids       <- x$ids[order]

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

      if (!not.na[i])
         next

      res <- try(suppressWarnings(rma.uni(yi.f, vi.f, weights=weights.f, intercept=TRUE, method=x$method, weighted=x$weighted, test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=seq_len(i))), silent=TRUE)

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
      out <- list(estimate=beta[not.na], se=se[not.na], zval=zval[not.na], pvals=pval[not.na], ci.lb=ci.lb[not.na], ci.ub=ci.ub[not.na], Q=QE[not.na], Qp=QEp[not.na], tau2=tau2[not.na], I2=I2[not.na], H2=H2[not.na])
      out$slab <- slab[not.na]
      out$ids  <- ids[not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(estimate=beta, se=se, zval=zval, pvals=pval, ci.lb=ci.lb, ci.ub=ci.ub, Q=QE, Qp=QEp, tau2=tau2, I2=I2, H2=H2)
      out$slab <- slab
      out$ids  <- ids
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   if (is.element(x$test, c("knha","adhoc","t")))
      names(out)[3] <- "tval"

   ### remove tau2 for FE/EE/CE models

   if (is.element(x$method, c("FE","EE","CE")))
      out <- out[-9]

   out$digits    <- digits
   out$transf    <- transf
   out$slab.null <- x$slab.null
   out$level     <- x$level
   out$measure   <- x$measure
   out$test      <- x$test

   attr(out$estimate, "measure") <- x$measure

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(out) <- c("list.rma", "cumul.rma")
   return(out)

}
