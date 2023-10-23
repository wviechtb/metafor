rstandard.rma.uni <- function(model, digits, type="marginal", ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(model), must="rma.uni", notav=c("robust.rma", "rma.gen", "rma.uni.selmodel"))

   na.act <- getOption("na.action")
   on.exit(options(na.action=na.act), add=TRUE)

   if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", "na.pass")))
      stop(mstyle$stop("Unknown 'na.action' specified under options()."))

   type <- match.arg(type, c("marginal", "conditional"))

   x <- model

   if (type == "conditional" && (!is.null(x$weights) || !x$weighted))
      stop(mstyle$stop("Extraction of conditional residuals not available for models with non-standard weights."))

   #if (type == "conditional" & inherits(x, "robust.rma"))
   #   stop(mstyle$stop("Extraction of conditional residuals not available for objects of class \"robust.rma\"."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   #########################################################################

   options(na.action="na.omit")
   H <- hatvalues(x, type="matrix")
   options(na.action = na.act)

   #########################################################################

   ImH <- diag(x$k) - H
   #ei <- ImH %*% cbind(x$yi)

   if (type == "marginal") {

      ei <- c(x$yi - x$X %*% x$beta)

      ei[abs(ei) < 100 * .Machine$double.eps] <- 0
      #ei[abs(ei) < 100 * .Machine$double.eps * median(abs(ei), na.rm=TRUE)] <- 0 ### see lm.influence

      ### don't allow this; the SEs of the residuals cannot be estimated consistently for "robust.rma" objects
      #if (inherits(x, "robust.rma")) {
      #   ve <- ImH %*% tcrossprod(x$meat,ImH)
      #} else {
         #ve <- ImH %*% tcrossprod(x$M,ImH)
      #}

      ve <- ImH %*% tcrossprod(x$M,ImH)

      #ve <- x$M + x$X %*% x$vb %*% t(x$X) - 2*H%*%x$M
      sei <- sqrt(diag(ve))

   }

   if (type == "conditional") {

      li <- x$tau2 / (x$tau2 + x$vi)

      pred  <- rep(NA_real_, x$k)

      for (i in seq_len(x$k)) {
         Xi <- matrix(x$X[i,], nrow=1)
         pred[i] <- li[i] * x$yi[i] + (1 - li[i]) * Xi %*% x$beta
      }

      ei <- x$yi - pred
      sei <- sqrt(x$vi^2 * 1/(x$vi + x$tau2) * (1 - diag(H)))

   }

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
