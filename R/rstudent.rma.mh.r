rstudent.rma.mh <- function(model, digits, progbar=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(model, "rma.mh"))
      stop(mstyle$stop("Argument 'model' must be an object of class \"rma.mh\"."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- model

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   #########################################################################

   delpred  <- rep(NA_real_, x$k.f)
   vdelpred <- rep(NA_real_, x$k.f)

   ### note: skipping NA tables

   if (progbar)
      pbar <- txtProgressBar(min=0, max=x$k.f, style=3)

   for (i in seq_len(x$k.f)) {

      if (progbar)
         setTxtProgressBar(pbar, i)

      if (!x$not.na[i])
         next

      if (is.element(x$measure, c("RR","OR","RD"))) {
         res <- try(suppressWarnings(rma.mh(ai=x$ai.f, bi=x$bi.f, ci=x$ci.f, di=x$di.f, measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct, subset=-i)), silent=TRUE)
      } else {
         res <- try(suppressWarnings(rma.mh(x1i=x$x1i.f, x2i=x$x2i.f, t1i=x$t1i.f, t2i=x$t2i.f, measure=x$measure, add=x$add, to=x$to, drop00=x$drop00, correct=x$correct, subset=-i)), silent=TRUE)
      }

      if (inherits(res, "try-error"))
         next

      delpred[i]  <- res$beta
      vdelpred[i] <- res$vb

   }

   if (progbar)
      close(pbar)

   resid <- x$yi.f - delpred
   resid[abs(resid) < 100 * .Machine$double.eps] <- 0
   #resid[abs(resid) < 100 * .Machine$double.eps * median(abs(resid), na.rm=TRUE)] <- 0 ### see lm.influence
   seresid <- sqrt(x$vi.f + vdelpred)
   stresid <- resid / seresid

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(resid=resid[x$not.na.yivi], se=seresid[x$not.na.yivi], z=stresid[x$not.na.yivi])
      out$slab <- x$slab[x$not.na.yivi]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(resid=resid, se=seresid, z=stresid)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na.yivi))
      stop(mstyle$stop("Missing values in results."))

   out$digits <- digits

   class(out) <- "list.rma"
   return(out)

}
