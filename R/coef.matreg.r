coef.matreg <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="matreg")

   coefs <- c(object$tab$beta)
   names(coefs) <- rownames(object$tab)

   return(coefs)

}
