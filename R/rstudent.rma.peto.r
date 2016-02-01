rstudent.rma.peto <- function(model, digits, ...) {

   if (!inherits(model, "rma.peto"))
      stop("Argument 'model' must be an object of class \"rma.peto\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   x <- model

   if (missing(digits))
      digits <- x$digits

   #########################################################################

   delpred  <- rep(NA_real_, x$k.f)
   vdelpred <- rep(NA_real_, x$k.f)

   ### note: skipping NA tables

   for (i in seq_len(x$k.f)[x$not.na]) {

      res <- try(suppressWarnings(rma.peto(ai=x$ai.f, bi=x$bi.f, ci=x$ci.f, di=x$di.f, add=x$add, to=x$to, drop00=x$drop00, subset=-i)), silent=TRUE)

      if (inherits(res, "try-error"))
         next

      delpred[i]  <- res$b
      vdelpred[i] <- res$vb

   }

   delresid <- x$yi.f - delpred
   delresid[abs(delresid) < 100 * .Machine$double.eps] <- 0
   #delresid[abs(delresid) < 100 * .Machine$double.eps * median(abs(delresid), na.rm=TRUE)] <- 0 ### see lm.influence
   sedelresid <- sqrt(x$vi.f + vdelpred)
   standelres <- delresid / sedelresid

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(resid=delresid[x$not.na.yivi], se=sedelresid[x$not.na.yivi], z=standelres[x$not.na.yivi])
      out$slab <- x$slab[x$not.na.yivi]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(resid=delresid, se=sedelresid, z=standelres)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na.yivi))
      stop("Missing values in results.")

   out$digits <- digits

   class(out) <- "list.rma"
   return(out)

}
