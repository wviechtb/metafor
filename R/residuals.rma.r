residuals.rma <- function(object, type="response", ...) {

   if (!inherits(object, "rma"))
      stop("Argument 'object' must be an object of class \"rma\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   type <- match.arg(type, c("response", "rstandard", "rstudent", "pearson", "cholesky"))

   ### for objects of class "rma.mh" and "rma.peto", use rstandard() to get the Pearson residuals

   if (inherits(object, c("rma.mh", "rma.peto")) && type == "pearson")
      type <- "rstandard"

   #########################################################################

   if (type == "rstandard") {
      tmp <- rstandard(object)
      out <- c(tmp$z)
      names(out) <- tmp$slab
   }

   if (type == "rstudent") {
      tmp <- rstudent(object)
      out <- c(tmp$z)
      names(out) <- tmp$slab
   }

   #########################################################################

   if (type == "response") {

      ### note: can calculate this even if vi is missing

      out <- c(object$yi.f - object$X.f %*% object$beta)
      out[abs(out) < 100 * .Machine$double.eps] <- 0

   }

   if (type == "pearson") {

      if (inherits(object, c("rma.mh", "rma.peto", "rma.glmm")))
         stop("Extraction of Pearson residuals not implemented for objects of class \"rma.mh\", \"rma.peto\", or \"rma.glmm\".")

      out <- c(object$yi.f - object$X.f %*% object$beta)
      out[abs(out) < 100 * .Machine$double.eps] <- 0

      se <- rep(NA_real_, object$k.f)
      se[object$not.na] <- sqrt(diag(object$M))

      out <- out / se

   }

   if (type == "cholesky") {

      ### note: Cholesky residuals depend on the data order

      if (inherits(object, c("rma.mh", "rma.peto", "rma.glmm")))
         stop("Extraction of Cholesky residuals not implemented for objects of class \"rma.mh\", \"rma.peto\", or \"rma.glmm\".")

      out <- c(object$yi - object$X %*% object$beta)
      out[abs(out) < 100 * .Machine$double.eps] <- 0

      L <- try(chol(chol2inv(chol(object$M))))
      tmp <- L %*% out

      out <- rep(NA_real_, object$k.f)
      out[object$not.na] <- tmp

   }

   if (is.element(type, c("response", "pearson", "cholesky"))) {

      names(out) <- object$slab

      not.na <- !is.na(out)

      if (na.act == "na.omit")
         out <- out[not.na]

      if (na.act == "na.fail" && any(!not.na))
         stop("Missing values in results.")

   }

   #########################################################################

   return(out)

}
