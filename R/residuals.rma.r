residuals.rma <- function(object, ...) {

   if (!is.element("rma", class(object)))
      stop("Argument 'object' must be an object of class \"rma\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   #########################################################################

   x <- object

   ### note: can calculate this even if vi is missing

   out <- c(x$yi.f - x$X.f %*% x$b)
   out[abs(out) < 100 * .Machine$double.eps] <- 0

   #########################################################################

   names(out) <- x$slab

   not.na <- !is.na(out)

   if (na.act == "na.omit")
      out <- out[not.na]

   if (na.act == "na.fail" && any(!not.na))
      stop("Missing values in results.")

   return(out)

}
