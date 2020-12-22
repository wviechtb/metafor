residuals.rma <- function(object, type="response", ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

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
         stop(mstyle$stop("Extraction of Pearson residuals not available for objects of class \"rma.mh\", \"rma.peto\", or \"rma.glmm\"."))

      out <- c(object$yi.f - object$X.f %*% object$beta)
      out[abs(out) < 100 * .Machine$double.eps] <- 0

      se <- rep(NA_real_, object$k.f)
      se[object$not.na] <- sqrt(diag(object$M))

      out <- out / se

   }

   if (type == "cholesky") {

      ### note: Cholesky residuals depend on the data order
      ### but only for the Cholesky residuals is QE = sum(residuals(res, type="cholesky)^2) for models where M (or rather: V) is not diagonal

      if (inherits(object, c("rma.mh", "rma.peto", "rma.glmm")))
         stop(mstyle$stop("Extraction of Cholesky residuals not available for objects of class \"rma.mh\", \"rma.peto\", or \"rma.glmm\"."))

      out <- c(object$yi - object$X %*% object$beta)
      out[abs(out) < 100 * .Machine$double.eps] <- 0

      L <- try(chol(chol2inv(chol(object$M))))

      if (inherits(L, "try-error"))
         stop(mstyle$stop("Could not take Cholesky decomposition of the marginal var-cov matrix."))

      tmp <- L %*% out

      out <- rep(NA_real_, object$k.f)
      out[object$not.na] <- tmp

   }

   if (is.element(type, c("response", "pearson", "cholesky"))) {

      names(out) <- object$slab

      #not.na <- !is.na(out)

      if (na.act == "na.omit")
         out <- out[object$not.na]

      if (na.act == "na.exclude")
         out[!object$not.na] <- NA

      if (na.act == "na.fail" && any(!object$not.na))
         stop(mstyle$stop("Missing values in results."))

   }

   #########################################################################

   return(out)

}
