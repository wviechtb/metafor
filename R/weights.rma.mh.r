weights.rma.mh <- function(object, type="diagonal", ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma.mh")

   if (is.null(object$outdat))
      stop(mstyle$stop("Information needed to compute the weights is not available in the model object."))

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("diagonal", "matrix"))

   x <- object

   #########################################################################

   if (is.element(x$measure, c("RR","OR","RD"))) {
      Ni <- with(x$outdat, ai + bi + ci + di)
   } else {
      Ti <- with(x$outdat, t1i + t2i)
   }

   if (x$measure == "OR")
      wi <- with(x$outdat, (bi / Ni) * ci)

   if (x$measure == "RR")
      wi <- with(x$outdat, (ci / Ni) * (ai+bi))

   if (x$measure == "RD")
      wi <- with(x$outdat, ((ai+bi) / Ni) * (ci+di))

   if (x$measure == "IRR")
      wi <- with(x$outdat, (x2i / Ti) * t1i)

   if (x$measure == "IRD")
      wi <- with(x$outdat, (t1i / Ti) * t2i)

   #########################################################################

   if (type == "diagonal") {

      weight <- rep(NA_real_, x$k.f)
      weight[x$not.na] <- wi / sum(wi) * 100
      names(weight) <- x$slab

      if (na.act == "na.omit")
         weight <- weight[x$not.na]

      if (na.act == "na.fail" && any(!x$not.na))
         stop(mstyle$stop("Missing values in weights."))

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
         stop(mstyle$stop("Missing values in results."))

      return(Wfull)

   }

}
