weights.rma.uni <- function(object, type="diagonal", ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="rma.uni", notav="rma.uni.selmodel")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("diagonal", "matrix"))

   x <- object

   #########################################################################

   if (x$weighted) {
      if (is.null(x$weights)) {
         W <- diag(1/(x$vi + x$tau2), nrow=x$k, ncol=x$k)
      } else {
         W <- diag(x$weights, nrow=x$k, ncol=x$k)
      }
   } else {
      W <- diag(1/x$k, nrow=x$k, ncol=x$k)
   }

   #########################################################################

   if (type == "diagonal") {

      wi <- as.vector(diag(W))
      weight <- rep(NA_real_, x$k.f)
      weight[x$not.na] <- wi/sum(wi) * 100
      names(weight) <- x$slab

      if (na.act == "na.omit")
         weight <- weight[x$not.na]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in weights."))

      return(weight)

   }

   if (type == "matrix") {

      Wfull <- matrix(NA_real_, nrow=x$k.f, ncol=x$k.f)
      Wfull[x$not.na, x$not.na] <- W

      rownames(Wfull) <- x$slab
      colnames(Wfull) <- x$slab

      if (na.act == "na.omit")
         Wfull <- Wfull[x$not.na, x$not.na, drop=FALSE]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in results."))

      return(Wfull)

   }

}
