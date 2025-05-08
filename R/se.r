se <- function(object, ...)
   UseMethod("se")

se.default <- function(object, ...) {

   vb <- try(vcov(object, ...), silent=TRUE)

   if (inherits(vb, "try-error") || !is.matrix(vb) || !.is.square(vb))
      stop("Default method for extracting the standard errors does not work for such model objects.")

   return(sqrt(diag(vb)))

}

se.rma <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="rma")

   ddd <- list(...)

   ses <- c(object$se)
   names(ses) <- rownames(object$beta)

   if (isTRUE(ddd$type=="beta"))
      return(ses)

   if (inherits(object, "rma.ls")) {
      ses <- list(beta=ses)
      ses$alpha <- c(object$se.alpha)
      names(ses$alpha) <- rownames(object$alpha)
      if (isTRUE(ddd$type=="alpha"))
         return(ses$alpha)
   }

   if (inherits(object, "rma.uni.selmodel")) {
      ses <- list(beta=ses)
      ses$delta <- c(object$se.delta)
      if (length(object$delta) == 1L) {
         names(ses$delta) <- "delta"
      } else {
         names(ses$delta) <- paste0("delta.", seq_along(object$delta))
      }
      if (isTRUE(ddd$type=="delta"))
         return(ses$delta)
   }

   return(ses)

}
