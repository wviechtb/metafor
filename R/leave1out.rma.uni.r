leave1out.rma.uni <- function(x, digits, transf, targs, progbar=FALSE, ...) {

   if (!inherits(x, "rma.uni"))
      stop("Argument 'x' must be an object of class \"rma.uni\".")

   if (inherits(x, "robust.rma"))
      stop("Method not yet implemented for objects of class \"robust.rma\". Sorry!")

   if (inherits(x, "rma.ls"))
      stop("Method not yet implemented for objects of class \"rma.ls\". Sorry!")

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
      pbar <- txtProgressBar(min=0, max=x$k.f, style=3)

   for (i in seq_len(x$k.f)) {

      if (progbar)
         setTxtProgressBar(pbar, i)

      if (!x$not.na[i])
         next

      res <- try(suppressWarnings(rma.uni(x$yi.f, x$vi.f, weights=x$weights.f, intercept=TRUE, method=x$method, weighted=x$weighted, test=x$test, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=-i)), silent=TRUE)

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
      close(pbar)

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
      stop("Missing values in results.")

   if (x$method == "FE")
      out <- out[-c(9,10,11)]

   #out <- out[-c(2,3,4)]
   out$digits <- digits
   out$transf <- transf

   class(out) <- "list.rma"
   return(out)

}
