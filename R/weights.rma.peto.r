weights.rma.peto <- function(object, type="diagonal", ...) {

   if (!inherits(object, "rma.peto"))
      stop("Argument 'object' must be an object of class \"rma.peto\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   type <- match.arg(type, c("diagonal", "matrix"))

   x <- object

   #########################################################################

   n1i <- x$ai + x$bi
   n2i <- x$ci + x$di
   Ni  <- x$ai + x$bi + x$ci + x$di
   xt  <- x$ai + x$ci
   yt  <- x$bi + x$di

   wi  <- xt * yt * (n1i/Ni) * (n2i/Ni) / (Ni - 1)

   #########################################################################

   if (type == "diagonal") {

      weight <- rep(NA_real_, x$k.f)
      weight[x$not.na] <- wi/sum(wi) * 100
      names(weight) <- x$slab

      if (na.act == "na.omit")
         weight <- weight[x$not.na]

      if (na.act == "na.fail" && any(!x$not.na))
         stop("Missing values in weights.")

      return(weight)

   }

   if (type == "matrix") {

      Wfull <- matrix(NA_real_, nrow=x$k.f, ncol=x$k.f)
      Wfull[x$not.na, x$not.na] <- diag(wi)

      rownames(Wfull) <- x$slab
      colnames(Wfull) <- x$slab

      if (na.act == "na.omit")
         Wfull <- Wfull[x$not.na, x$not.na, drop=FALSE]

      if (na.act == "na.fail" && any(!x$not.na))
         stop("Missing values in results.")

      return(Wfull)

   }

}
