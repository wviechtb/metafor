rstandard.rma.peto <- function(model, digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(model), must="rma.peto")

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

   resid <- c(x$yi.f - x$beta)

   resid[abs(resid) < 100 * .Machine$double.eps] <- 0
   #resid[abs(resid) < 100 * .Machine$double.eps * median(abs(resid), na.rm=TRUE)] <- 0 # see lm.influence

   ### note: these are like Pearson (or semi-standardized) residuals

   seresid <- sqrt(x$vi.f)
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
