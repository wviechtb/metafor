rstandard.rma.peto <- function(model, digits, ...) {

   if (!inherits(model, "rma.peto"))
      stop("Argument 'model' must be an object of class \"rma.peto\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   x <- model

   if (missing(digits))
      digits <- x$digits

   #########################################################################

   ei <- c(x$yi.f - x$b)
   ei[abs(ei) < 100 * .Machine$double.eps] <- 0
   #ei[abs(ei) < 100 * .Machine$double.eps * median(abs(ei), na.rm=TRUE)] <- 0 ### see lm.influence
   sei <- sqrt(x$vi.f)
   zi <- ei / sei

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(resid=ei[x$not.na.yivi], se=sei[x$not.na.yivi], z=zi[x$not.na.yivi])
      out$slab <- x$slab[x$not.na.yivi]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(resid=ei, se=sei, z=zi)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na.yivi))
      stop("Missing values in results.")

   out$digits <- digits

   class(out) <- "list.rma"
   return(out)

}
