vif.rma <- function(x, intercept=FALSE, table=FALSE, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

   if (inherits(x, "robust.rma"))
      stop(mstyle$stop("Function not applicable to objects of class \"robust.rma\"."))

   if (x$int.only)
      stop(mstyle$stop("VIF not applicable for intercept-only models."))

   if (missing(digits))
      digits <- x$digits

   #########################################################################

   vb <- vcov(x)

   if (inherits(x, "rma.ls"))
      vb <- vb$vb

   ### remove intercept row/colum from vb if model includes one and intercept=FALSE
   if (x$intercept && !intercept)
     vb <- vb[-1,-1,drop=FALSE]

   ### rescale vb to correlation matrix
   rb <- cov2cor(vb)

   ### try computing the VIFs
   vif <- try(diag(chol2inv(chol(rb))), silent=TRUE)

   if (inherits(vif, "try-error"))
      stop(mstyle$stop("Cannot invert var-cov matrix to compute VIFs."))

   ### add NA for intercept if model includes one and intercept=FALSE
   if (x$intercept && !intercept && table)
      vif <- c(NA, vif)

   if (table) {
      vif <- cbind(coef(summary(x)), vif=vif)
   } else {
      names(vif) <- rownames(vb)
   }

   vif <- round(vif, digits=digits)

   return(vif)

}
