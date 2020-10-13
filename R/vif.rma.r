vif.rma <- function(x, btt, intercept=FALSE, table=FALSE, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma\"."))

   if (inherits(x, "robust.rma"))
      stop(mstyle$stop("Method not available for objects of class \"robust.rma\"."))

   if (inherits(x, "rma.uni.selmodel"))
      stop(mstyle$stop("Method not available for objects of class \"rma.uni.selmodel\"."))

   if (x$int.only)
      stop(mstyle$stop("VIF not applicable for intercept-only models."))

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   #########################################################################

   vb <- vcov(x)

   if (inherits(x, "rma.ls"))
      vb <- vb$vb

   if (missing(btt)) {

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
         vif <- cbind(coef(summary(x)), vif=vif, sif=sqrt(vif))
      } else {
         names(vif) <- rownames(vb)
      }

      res <- list(vif=vif, digits=digits, table=table, test=x$test)

   } else {

      btt <- .set.btt(btt, x$p, x$int.incl, colnames(x$X))

      if (x$intercept && !intercept) {
         vb <- vb[-1,-1,drop=FALSE]
         btt <- btt - 1
         if (any(btt < 1))
            warning(mstyle$warning("Intercept term not included in GVIF computation."), call.=FALSE)
         btt <- btt[btt > 0]
      }
      m <- length(btt)

      rb <- cov2cor(vb)
      detrv <- det(rb)
      gvif <- det(as.matrix(rb[btt,btt])) * det(as.matrix(rb[-btt,-btt])) / detrv
      gsif <- gvif^(1/(2*m))

      if (x$intercept && !intercept)
         btt <- btt + 1

      res <- list(gvif=gvif, gsif=gsif, btt=btt, m=m, digits=digits)

   }

   class(res) <- "vif.rma"
   return(res)

}
