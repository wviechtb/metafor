coef.rma <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "rma"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"rma\"."))

   coefs <- c(object$beta)
   names(coefs) <- rownames(object$beta)

   if (inherits(object, "rma.ls")) {
      coefs <- list(beta=coefs)
      coefs$alpha <- c(object$alpha)
      names(coefs$alpha) <- rownames(object$alpha)
   }

   if (inherits(object, "rma.uni.selmodel")) {
      coefs <- list(beta=coefs)
      coefs$delta <- c(object$delta)
      if (length(object$delta) == 1L) {
         names(coefs$delta) <- "delta"
      } else {
         names(coefs$delta) <- paste0("delta.", 1:length(object$delta))
      }
   }

   return(coefs)

}
