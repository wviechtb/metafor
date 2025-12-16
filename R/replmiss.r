replmiss <- function(x, y, data) {

   mstyle <- .get.mstyle()

   # check if the 'data' argument was specified

   if (missing(data))
      data <- NULL

   if (is.null(data)) {
      data <- sys.frame(sys.parent())
   } else {
      if (!is.data.frame(data))
         data <- data.frame(data)
   }

   mf <- match.call()

   x <- .getx("x", mf=mf, data=data, checknull=FALSE)
   y <- .getx("y", mf=mf, data=data, checknull=FALSE)

   if (length(y) == 0L)
      return(x)

   if (length(x) == 0L)
      x <- rep(NA_real_, length(y))

   # in case the user specified a constant for 'y' to use for the replacement

   y <- .expand1(y, length(x))

   # check that 'x' and 'y' are of the same length

   if (length(x) != length(y))
      stop(mstyle$stop("Length of 'x' and 'y' are not the same."))

   #x <- ifelse(is.na(x), y, x) # this is quite a bit slower than the following
   is.na.x <- is.na(x)
   x[is.na.x] <- y[is.na.x]

   return(x)

}
