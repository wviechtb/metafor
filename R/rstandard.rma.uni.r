rstandard.rma.uni <- function(model, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(model, "rma.uni"))
      stop("Argument 'model' must be an object of class \"rma.uni\".")

   na.act <- getOption("na.action")

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   x <- model

   if (missing(digits))
      digits <- x$digits

   #########################################################################

   options(na.action="na.omit")
   H <- hatvalues(x, type="matrix")
   options(na.action = na.act)

   #########################################################################

   ImH <- diag(x$k) - H
   #ei <- ImH %*% cbind(x$yi)
   ei <- c(x$yi - x$X %*% x$beta)

   ei[abs(ei) < 100 * .Machine$double.eps] <- 0
   #ei[abs(ei) < 100 * .Machine$double.eps * median(abs(ei), na.rm=TRUE)] <- 0 ### see lm.influence

   if (inherits(x, "robust.rma")) {
      ve <- ImH %*% tcrossprod(x$meat,ImH)
   } else {
      ve <- ImH %*% tcrossprod(x$M,ImH)
   }

   #ve <- x$M + x$X %*% x$vb %*% t(x$X) - 2*H%*%x$M
   sei <- sqrt(diag(ve))

   resid   <- rep(NA_real_, x$k.f)
   seresid <- rep(NA_real_, x$k.f)
   stresid <- rep(NA_real_, x$k.f)

   resid[x$not.na]   <- ei
   seresid[x$not.na] <- sei
   stresid[x$not.na] <- ei / sei

   #########################################################################

   if (na.act == "na.omit") {
      out <- list(resid=resid[x$not.na], se=seresid[x$not.na], z=stresid[x$not.na])
      out$slab <- x$slab[x$not.na]
   }

   if (na.act == "na.exclude" || na.act == "na.pass") {
      out <- list(resid=resid, se=seresid, z=stresid)
      out$slab <- x$slab
   }

   if (na.act == "na.fail" && any(!x$not.na))
      stop(mstyle$stop("Missing values in results."))

   out$digits <- digits

   class(out) <- "list.rma"
   return(out)

}
