weights.rma.mh <- function(object, type="diagonal", ...) {

   if (!is.element("rma.mh", class(object)))
      stop("Argument 'object' must be an object of class \"rma.mh\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop("Unknown 'na.action' specified under options().")

   type <- match.arg(type, c("diagonal", "matrix"))

   x <- object

   #########################################################################

   if (is.element(x$measure, c("RR","OR","RD"))) {
      Ni <- x$ai + x$bi + x$ci + x$di
   } else {
      Ti <- x$t1i + x$t2i
   }

   if (x$measure == "OR")
      wi <- (x$bi / Ni) * x$ci

   if (x$measure == "RR")
      wi <- (x$ci / Ni) * (x$ai+x$bi)

   if (x$measure == "RD")
      wi <- ((x$ai+x$bi) / Ni) * (x$ci+x$di)

   if (x$measure == "IRR")
      wi <- (x$x2i / Ti) * x$t1i

   if (x$measure == "IRD")
      wi <- (x$t1i / Ti) * x$t2i

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
