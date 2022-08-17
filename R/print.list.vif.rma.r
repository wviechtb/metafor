print.list.vif.rma <- function(x, digits=x[[1]]$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="list.vif.rma")

   digits <- .get.digits(digits=digits, xdigits=x[[1]]$digits, dmiss=FALSE)

   .space()

   cat(mstyle$section(paste0("Location Coefficients:\n")))

   print(x[[1]], digits=digits, ...)

   .space(FALSE)

   cat(mstyle$section(paste0("Scale Coefficients:\n")))

   print(x[[2]], digits=digits, ...)

   invisible()

}
