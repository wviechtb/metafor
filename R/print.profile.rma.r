print.profile.rma <- function(x, ...) {

   #########################################################################

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "profile.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"profile.rma\"."))

   #########################################################################

   if (x$comps == 1) {

      res <- data.frame(x[1], x[2])
      print(res)

   } else {

      x$comps <- NULL
      print(lapply(x, function(x) data.frame(x[1], x[2])))

   }

}
