cumul.rma.uni <- function(x, order, digits, transf, targs, collapse=FALSE, progbar=FALSE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma.uni", notav=c("robust.rma", "rma.ls", "rma.gen", "rma.uni.selmodel"))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in data."))

   if (!x$int.only)
      stop(mstyle$stop("Method only applicable to models without moderators."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   if (missing(transf))
      transf <- FALSE

   if (missing(targs))
      targs <- NULL

   funlist <- lapply(list(transf.exp.int, transf.ilogit.int, transf.iprobit.int, transf.ztor.int, transf.iarcsin.int, transf.iahw.int, transf.iabt.int, transf.exp.mode, transf.ilogit.mode, transf.iprobit.mode, transf.ztor.mode, transf.iarcsin.mode, transf.iahw.mode, transf.iabt.mode), deparse)

   if (is.null(targs) && any(sapply(funlist, identical, deparse(transf))) && inherits(x, c("rma.uni","rma.glmm")) && length(x$tau2 == 1L))
      targs <- list(tau2=x$tau2)

   ddd <- list(...)

   .chkdots(ddd, c("time", "decreasing", "code1", "code2"))

   if (.isTRUE(ddd$time))
      time.start <- proc.time()

   decreasing <- .chkddd(ddd$decreasing, FALSE)

   #########################################################################

   if (grepl("^order\\(", deparse1(substitute(order))))
      warning(mstyle$warning("Use of order() in the 'order' argument is probably erroneous."), call.=FALSE)

   if (missing(order)) {

      orvar <- seq_len(x$k.all)
      collapse <- FALSE

   } else {

      mf <- match.call()
      orvar <- .getx("order", mf=mf, data=x$data)

      if (length(orvar) != x$k.all)
         stop(mstyle$stop(paste0("Length of the 'order' argument (", length(orvar), ") does not correspond to the size of the original dataset (", x$k.all, ").")))

   }

   ### note: order variable must be of the same length as the original dataset
   ###       so apply the same subsetting as was done during the model fitting

   orvar <- .getsubset(orvar, x$subset)

   ### order data by the order variable (NAs in order variable are dropped)

   order <- base::order(orvar, decreasing=decreasing, na.last=NA)

   yi      <- x$yi.f[order]
   vi      <- x$vi.f[order]
   weights <- x$weights.f[order]
   not.na  <- x$not.na[order]
   slab    <- x$slab[order]
   ids     <- x$ids[order]
   orvar   <- orvar[order]

   if (inherits(x$data, "environment")) {
      data <- NULL
   } else {
      data <- x$data[order,,drop=FALSE]
   }

   if (collapse) {
      uorvar <- unique(orvar)
   } else {
      uorvar <- orvar
   }

   k.o <- length(uorvar)

   if (!is.null(ddd[["code1"]]))
      eval(expr = parse(text = ddd[["code1"]]))

   k     <- rep(NA_integer_, k.o)
   beta  <- rep(NA_real_, k.o)
   se    <- rep(NA_real_, k.o)
   zval  <- rep(NA_real_, k.o)
   pval  <- rep(NA_real_, k.o)
   ci.lb <- rep(NA_real_, k.o)
   ci.ub <- rep(NA_real_, k.o)
   QE    <- rep(NA_real_, k.o)
   QEp   <- rep(NA_real_, k.o)
   tau2  <- rep(NA_real_, k.o)
   I2    <- rep(NA_real_, k.o)
   H2    <- rep(NA_real_, k.o)
   show  <- rep(TRUE, k.o)

   ### elements that need to be returned

   outlist <- "k=k, beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, QE=QE, QEp=QEp, tau2=tau2, I2=I2, H2=H2"

   if (progbar)
      pbar <- pbapply::startpb(min=0, max=k.o)

   for (i in seq_len(k.o)) {

      if (progbar)
         pbapply::setpb(pbar, i)

      if (!is.null(ddd[["code2"]]))
         eval(expr = parse(text = ddd[["code2"]]))

      if (collapse) {

         if (all(!not.na[is.element(orvar, uorvar[i])])) {
            if (na.act == "na.omit")
               show[i] <- FALSE # if all studies to be added are !not.na, don't show (but a fit failure is still shown)
            next
         }

         incl <- is.element(orvar, uorvar[1:i])

      } else {

         if (!not.na[i]) {
            if (na.act == "na.omit")
               show[i] <- FALSE # if study to be added is !not.na, don't show (but a fit failure is still shown)
            next
         }

         incl <- 1:i

      }

      args <- list(yi=yi, vi=vi, weights=weights, intercept=TRUE, method=x$method, weighted=x$weighted,
                   test=x$test, level=x$level, tau2=ifelse(x$tau2.fix, x$tau2, NA), control=x$control, subset=incl, outlist=outlist)

      res <- try(suppressWarnings(.do.call(rma.uni, args)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      k[i]     <- res$k
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
         se    <- rep(NA_real_, k.o)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         if (!is.primitive(transf) && !is.null(targs) && length(formals(transf)) == 1L)
            stop(mstyle$stop("Function specified via 'transf' does not appear to have an argument for 'targs'."))
         beta  <- sapply(beta, transf, targs)
         se    <- rep(NA_real_, k.o)
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

   out <- list(k=k[show], estimate=beta[show], se=se[show], zval=zval[show], pval=pval[show], ci.lb=ci.lb[show], ci.ub=ci.ub[show], Q=QE[show], Qp=QEp[show], tau2=tau2[show], I2=I2[show], H2=H2[show])

   if (collapse) {
      out$slab <- uorvar[show]
      out$slab.null <- FALSE
   } else {
      out$slab <- slab[show]
      out$ids  <- ids[show]
      out$data <- data[show,,drop=FALSE]
      out$slab.null <- x$slab.null
   }

   out$order <- uorvar[show]

   if (is.element(x$test, c("knha","adhoc","t")))
      names(out)[4] <- "tval"

   ### remove tau2 for FE/EE/CE models

   if (is.element(x$method, c("FE","EE","CE")))
      out <- out[-10]

   out$digits <- digits
   out$transf <- transf
   out$level  <- x$level
   out$test   <- x$test

   if (!transf) {
      out$measure <- x$measure
      attr(out$estimate, "measure") <- x$measure
   }

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(out) <- c("list.rma", "cumul.rma")
   return(out)

}
