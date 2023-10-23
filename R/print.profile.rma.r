print.profile.rma <- function(x, ...) {

   #########################################################################

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="profile.rma")

   #########################################################################

   if (x$comps == 1) {

      res <- data.frame(x[1], x[2])
      print(res)

   } else {

      x$comps <- NULL
      print(lapply(x, function(x) data.frame(x[1], x[2])))

   }

}
