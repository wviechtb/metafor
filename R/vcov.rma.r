vcov.rma <- function(object, type="fixed", ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma")

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("fixed", "obs", "fitted", "resid"))

   #########################################################################

   if (type=="fixed") {

      out <- object$vb

      if (inherits(object, "rma.ls"))
         out <- list(beta = object$vb, alpha = object$va)

      return(out)

   }

   #########################################################################

   if (type=="obs") {

      if (inherits(object, c("rma.uni","rma.mv"))) {

         out <- matrix(NA_real_, nrow=object$k.f, ncol=object$k.f)
         out[object$not.na, object$not.na] <- as.matrix(object$M) # as.matrix() needed when sparse=TRUE

         rownames(out) <- colnames(out) <- object$slab

         if (na.act == "na.omit")
            out <- out[object$not.na, object$not.na]

         if (na.act == "na.fail" && any(!object$not.na))
            stop(mstyle$stop("Missing values in data."))

         return(out)

      } else {

         stop(mstyle$stop("Extraction of marginal var-cov matrix not available for objects of this class."))

      }

   }

   #########################################################################

   if (type=="fitted") {

      out <- object$X.f %*% object$vb %*% t(object$X.f)

      rownames(out) <- colnames(out) <- object$slab

      if (na.act == "na.omit")
         out <- out[object$not.na, object$not.na]

      if (na.act == "na.exclude" || na.act == "na.pass") {
         out[!object$not.na,] <- NA
         out[,!object$not.na] <- NA
      }

      return(out)

   }

   #########################################################################

   if (type=="resid") {

      ### the SEs of the residuals cannot be estimated consistently for "robust.rma" objects

      if (inherits(object, c("robust.rma", "rma.gen")))
         stop(mstyle$stop("Extraction of var-cov matrix of the residuals not available for objects of this type."))

      options(na.action="na.omit")
      H <- hatvalues(object, type="matrix")
      options(na.action = na.act)

      ImH <- diag(object$k) - H

      if (inherits(object, "robust.rma")) {
         ve <- ImH %*% tcrossprod(object$meat,ImH)
      } else {
         ve <- ImH %*% tcrossprod(as.matrix(object$M),ImH) # as.matrix() needed when sparse=TRUE
      }

      if (na.act == "na.omit") {
         out <- ve
         rownames(out) <- colnames(out) <- object$slab[object$not.na]
      }

      if (na.act == "na.exclude" || na.act == "na.pass") {
         out <- matrix(NA_real_, nrow=object$k.f, ncol=object$k.f)
         out[object$not.na, object$not.na] <- ve
         rownames(out) <- colnames(out) <- object$slab
      }

      return(out)

   }

   #########################################################################

}
