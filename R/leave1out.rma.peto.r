leave1out.rma.peto <- function(x, cluster, digits, transf, targs, progbar=FALSE, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma.peto")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   if (!x$int.only)
      stop(mstyle$stop("Method only applicable to models without moderators."))

   if (x$k == 1L)
      stop(mstyle$stop("Stopped because k = 1."))

   if (is.null(x$outdat.f))
      stop(mstyle$stop("Information needed to carry out a leave-one-out analysis is not available in the model object."))

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

   ### process cluster variable

   misscluster <- ifelse(missing(cluster), TRUE, FALSE)

   if (misscluster) {
      cluster <- seq_len(x$k.all)
   } else {
      mf <- match.call()
      cluster <- .getx("cluster", mf=mf, data=x$data)
   }

   ### note: cluster variable must be of the same length as the original dataset
   ###       so we have to apply the same subsetting (if necessary) and removing
   ###       of NAs as was done during model fitting

   if (length(cluster) != x$k.all)
      stop(mstyle$stop(paste0("Length of the variable specified via 'cluster' (", length(cluster), ") does not match the length of the data (", x$k.all, ").")))

   cluster <- .getsubset(cluster, x$subset)

   cluster.f <- cluster

   cluster <- cluster[x$not.na]

   ### checks on cluster variable

   if (anyNA(cluster.f))
      stop(mstyle$stop("No missing values allowed in 'cluster' variable."))

   if (length(cluster.f) == 0L)
      stop(mstyle$stop(paste0("Cannot find 'cluster' variable (or it has zero length).")))

   ### cluster ids and number of clusters

   ids <- unique(cluster)
   n <- length(ids)

   if (!misscluster)
      ids <- sort(ids)

   #########################################################################

   beta  <- rep(NA_real_, n)
   se    <- rep(NA_real_, n)
   zval  <- rep(NA_real_, n)
   pval  <- rep(NA_real_, n)
   ci.lb <- rep(NA_real_, n)
   ci.ub <- rep(NA_real_, n)
   QE    <- rep(NA_real_, n)
   QEp   <- rep(NA_real_, n)
   #tau2 <- rep(NA_real_, n)
   I2    <- rep(NA_real_, n)
   H2    <- rep(NA_real_, n)

   ### elements that need to be returned

   outlist <- "beta=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, QE=QE, QEp=QEp, tau2=tau2, I2=I2, H2=H2"

   if (progbar)
      pbar <- pbapply::startpb(min=0, max=n)

   for (i in seq_len(n)) {

      if (progbar)
         pbapply::setpb(pbar, i)

      args <- list(ai=x$outdat$ai, bi=x$outdat$bi, ci=x$outdat$ci, di=x$outdat$di,
                   add=x$add, to=x$to, drop00=x$drop00, level=x$level,
                   subset=ids[i]!=cluster, outlist=outlist)
      res <- try(suppressWarnings(.do.call(rma.peto, args)), silent=TRUE)

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
      I2[i]    <- res$I2
      H2[i]    <- res$H2

   }

   if (progbar)
      pbapply::closepb(pbar)

   #########################################################################

   ### if requested, apply transformation function

   if (.isTRUE(transf)) # if transf=TRUE, apply exp transformation to ORs
      transf <- exp

   if (is.function(transf)) {
      if (is.null(targs)) {
         beta  <- sapply(beta, transf)
         se    <- rep(NA_real_, n)
         ci.lb <- sapply(ci.lb, transf)
         ci.ub <- sapply(ci.ub, transf)
      } else {
         if (!is.primitive(transf) && !is.null(targs) && length(formals(transf)) == 1L)
            stop(mstyle$stop("Function specified via 'transf' does not appear to have an argument for 'targs'."))
         beta  <- sapply(beta, transf, targs)
         se    <- rep(NA_real_, n)
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

   out <- list(estimate=beta, se=se, zval=zval, pval=pval, ci.lb=ci.lb, ci.ub=ci.ub, Q=QE, Qp=QEp, I2=I2, H2=H2)

   if (na.act == "na.omit") {
      if (misscluster) {
         out$slab <- paste0("-", x$slab[x$not.na])
      } else {
         out$slab <- paste0("-", ids)
      }
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      if (misscluster) {
         out <- .expandna(out, x$not.na)
         out$slab <- paste0("-", x$slab)
      } else {
         out$slab <- paste0("-", ids)
      }
   }

   out$digits <- digits
   out$transf <- transf

   if (.isTRUE(ddd$time)) {
      time.end <- proc.time()
      .print.time(unname(time.end - time.start)[3])
   }

   class(out) <- "list.rma"
   return(out)

}
