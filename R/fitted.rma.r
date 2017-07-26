fitted.rma <- function(object, ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   ### note: fitted values can be calculated for all studies including those that
   ### have NAs on yi/vi; but if NA on X's, then the fitted value will also be NA

   out <- c(object$X.f %*% object$beta)
   names(out) <- object$slab

   #not.na <- !is.na(out)

   if (na.act == "na.omit")
      out <- out[object$not.na]
   if (na.act == "na.exclude")
      out[!object$not.na] <- NA
   if (na.act == "na.fail" && any(!object$not.na))
      stop("Missing values in results.")

   if (inherits(object, "rma.ls")) {

      out <- list(location = out)
      out$scale <- c(object$Z.f %*% object$alpha)

      names(out$scale) <- object$slab

      #not.na <- !is.na(out$scale)

      if (na.act == "na.omit")
         out$scale <- out$scale[object$not.na]
      if (na.act == "na.exclude")
         out$scale[!object$not.na] <- NA
      if (na.act == "na.fail" && any(!object$not.na))
         stop("Missing values in results.")

   }

   return(out)

}
