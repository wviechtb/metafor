vif.rma <- function(x, btt, intercept=FALSE, table=FALSE, digits, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma", notav=c("rma.uni.selmodel"))

   # allow for 'robust.rma' based on the same principle as used for standard 'rma.uni' models

   vif.loc   <- TRUE
   vif.scale <- TRUE

   if (inherits(x, "rma.ls")) {
      if (x$int.only)
         vif.loc <- FALSE
      if (x$Z.int.only)
         vif.scale <- FALSE
      if (!vif.loc && !vif.scale)
         stop(mstyle$stop("VIF not applicable for intercept-only models."))
   } else {
      if (x$int.only)
         stop(mstyle$stop("VIF not applicable for intercept-only models."))
   }

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("att"))

   #########################################################################

   if (vif.loc) {

      vcov <- vcov(x)

      if (inherits(x, "rma.ls"))
         vcov <- vcov$beta

      if (missing(btt) || is.null(btt)) {

         ### remove intercept row/colum from vcov if model includes one and intercept=FALSE
         if (x$intercept && !intercept)
            vcov <- vcov[-1,-1,drop=FALSE]

         ### rescale vcov to correlation matrix
         rb <- cov2cor(vcov)

         ### try computing the VIFs
         vif <- try(diag(chol2inv(chol(rb))), silent=TRUE)

         if (inherits(vif, "try-error"))
            stop(mstyle$stop("Cannot invert var-cov matrix to compute VIFs."))

         ### add NA for intercept if model includes one and intercept=FALSE
         if (x$intercept && !intercept && table)
            vif <- c(NA, vif)

         if (table) {
            tab <- coef(summary(x))
            if (inherits(x, "rma.ls"))
               tab <- tab$beta
            vif <- cbind(tab, vif=vif, sif=sqrt(vif))
         } else {
            names(vif) <- rownames(vcov)
         }

         res <- list(vif=vif, digits=digits, table=table, test=x$test)

      } else {

         btt <- .set.btt(btt, x$p, x$int.incl, colnames(x$X))

         if (x$intercept && !intercept) {
            vcov <- vcov[-1,-1,drop=FALSE]
            btt <- btt - 1
            if (any(btt < 1))
               warning(mstyle$warning("Intercept term not included in GVIF computation."), call.=FALSE)
            btt <- btt[btt > 0]
         }
         m <- length(btt)

         rb <- cov2cor(vcov)
         detrv <- det(rb)
         gvif <- det(as.matrix(rb[btt,btt])) * det(as.matrix(rb[-btt,-btt])) / detrv
         gsif <- gvif^(1/(2*m))

         if (x$intercept && !intercept)
            btt <- btt + 1

         res <- list(gvif=gvif, gsif=gsif, btt=btt, m=m, digits=digits)

      }

      class(res) <- "vif.rma"

   } else {

      res <- NULL

   }

   #########################################################################

   if (inherits(x, "rma.ls") && vif.scale) {

      res.loc <- res

      vcov <- vcov(x)$alpha

      if (is.null(ddd$att)) {

         ### remove intercept row/colum from vcov if model includes one and intercept=FALSE
         if (x$Z.intercept && !intercept)
            vcov <- vcov[-1,-1,drop=FALSE]

         ### rescale vcov to correlation matrix
         rb <- cov2cor(vcov)

         ### try computing the VIFs
         vif <- try(diag(chol2inv(chol(rb))), silent=TRUE)

         if (inherits(vif, "try-error"))
            stop(mstyle$stop("Cannot invert var-cov matrix to compute VIFs for the scale model."))

         ### add NA for intercept if model includes one and intercept=FALSE
         if (x$Z.intercept && !intercept && table)
            vif <- c(NA, vif)

         if (table) {
            tab <- coef(summary(x))$alpha
            vif <- cbind(tab, vif=vif, sif=sqrt(vif))
         } else {
            names(vif) <- rownames(vcov)
         }

         res.scale <- list(vif=vif, digits=digits, table=table, test=x$test)

      } else {

         att <- .set.btt(ddd$att, x$q, x$Z.int.incl, colnames(x$Z))

         if (x$Z.intercept && !intercept) {
            vcov <- vcov[-1,-1,drop=FALSE]
            att <- att - 1
            if (any(att < 1))
               warning(mstyle$warning("Intercept term not included in GVIF computation."), call.=FALSE)
            att <- att[att > 0]
         }
         m <- length(att)

         rb <- cov2cor(vcov)
         detrv <- det(rb)
         gvif <- det(as.matrix(rb[att,att])) * det(as.matrix(rb[-att,-att])) / detrv
         gsif <- gvif^(1/(2*m))

         if (x$Z.intercept && !intercept)
            att <- att + 1

         res.scale <- list(gvif=gvif, gsif=gsif, btt=att, m=m, digits=digits)

      }

      class(res.scale) <- "vif.rma"

      if (vif.loc) {
         if (vif.scale) {
            res <- list(beta=res.loc, alpha=res.scale)
         } else {
            res <- res.loc
         }
      } else {
         res <- res.scale
      }

   }

   #########################################################################

   return(res)

}
