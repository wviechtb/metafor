vcov.rma <- function(object, type="fixed", ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   type <- match.arg(type, c("fixed", "obs", "fitted"))

   #########################################################################

   if (type=="fixed") {

      out <- object$vb

      if (inherits(object, "rma.ls"))
         out <- list(location = object$vb, scale = object$vb.alpha)

      return(out)

   }

   #########################################################################

   if (type=="obs") {

      if (inherits(object, c("rma.uni","rma.mv"))) {

         if (na.act == "na.omit")
            return(object$M)

         if (na.act == "na.pass" || na.act == "na.exclude") {
            M <- matrix(NA_real_, nrow=object$k.f, ncol=object$k.f)
            M[object$not.na,object$not.na] <- object$M
            return(M)
         }

         if (na.act == "na.fail" && any(!object$not.na))
            stop("Missing values in data.")

      } else {

         stop("Extraction of marginal var-cov matrix not available for objects of this class.")

      }

   }

   #########################################################################

   if (type=="fitted") {

      out <- object$X.f %*% object$vb %*% t(object$X.f)

      #rownames(out) <- colnames(out) <- object$slab

      if (na.act == "na.omit")
         out <- out[object$not.na, object$not.na]

      if (na.act == "na.exclude") {
         out[!object$not.na,] <- NA
         out[,!object$not.na] <- NA
      }

      return(out)

   }

   #########################################################################

}
