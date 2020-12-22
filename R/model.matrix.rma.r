model.matrix.rma <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   ### note: lm() always returns X (never the full model matrix, even with na.exclude or na.pass)
   ### but it seems a bit more logical to actually return X.f in that case

   if (na.act == "na.omit")
      out <- object$X

   if (na.act == "na.exclude" || na.act == "na.pass")
      out <- object$X.f

   if (na.act == "na.fail" && any(!object$not.na))
      stop(mstyle$stop("Missing values in results."))

   if (inherits(object, "rma.ls")) {

      out <- list(location = out)

      if (na.act == "na.omit")
         out$scale <- object$Z

      if (na.act == "na.exclude" || na.act == "na.pass")
         out$scale <- object$Z.f

      if (na.act == "na.fail" && any(!object$not.na))
         stop(mstyle$stop("Missing values in results."))

   }

   return(out)

}
