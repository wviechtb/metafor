fitted.rma <- function(object, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   ### note: fitted values can be calculated for all studies including those that
   ### have NAs on yi/vi; but if NA on X's, then the fitted value will also be NA

   out <- object$X.f %*% object$b

   names(out) <- object$slab

   not.na <- !is.na(out)

   if (na.act == "na.omit")
      out <- out[not.na]

   if (na.act == "na.fail" && any(!not.na))
      stop("Missing values in results.")

   return(out)

}
