coef.rma <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   ddd <- list(...)

   coefs <- c(object$beta)
   names(coefs) <- rownames(object$beta)

   if (isTRUE(ddd$type=="beta"))
      return(coefs)

   if (inherits(object, "rma.ls")) {
      coefs <- list(beta=coefs)
      coefs$alpha <- c(object$alpha)
      names(coefs$alpha) <- rownames(object$alpha)
      if (isTRUE(ddd$type=="alpha"))
         return(coefs$alpha)
   }

   if (inherits(object, "rma.uni.selmodel")) {
      coefs <- list(beta=coefs)
      coefs$delta <- c(object$delta)
      if (length(object$delta) == 1L) {
         names(coefs$delta) <- "delta"
      } else {
         names(coefs$delta) <- paste0("delta.", seq_along(object$delta))
      }
      if (isTRUE(ddd$type=="delta"))
         return(coefs$delta)
   }

   return(coefs)

}
