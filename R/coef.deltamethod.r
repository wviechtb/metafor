coef.deltamethod <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="deltamethod")

   coefs <- c(object$tab$coef)
   names(coefs) <- rownames(object$tab)

   return(coefs)

}
