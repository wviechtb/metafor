# Note: works with "robust.rma" objects, but only type="fixed", not type="obs"
# (to make that work, would have to save M in the object returned by robust()
# and that seems unnecessary).

vcov.rma <- function(object, type="fixed", ...) {

   if (!is.element("rma", class(object)))
      stop("Argument 'object' must be an object of class \"rma\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   type <- match.arg(type, c("fixed", "obs"))

   #########################################################################

   if (type=="fixed")
      return(object$vb)

   #########################################################################

   if (type=="obs") {

      if (any(is.element(c("rma.uni","rma.mv"), class(object)))) {

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

}
