replmiss <- function(x, y) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (length(y) == 0L)
      return(x)

   if (length(x) == 0L)
      x <- rep(NA, length(y))

   ### in case user specifies a constant to use for replacement

   if (length(y) == 1L)
      y <- rep(y, length(x))

   ### check that x and y are of the same length

   if (length(x) != length(y))
      stop(mstyle$stop("Length of 'x' and 'y' is not the same."))

   #x <- ifelse(is.na(x), y, x) # this is quite a bit slower than the following
   is.na.x <- is.na(x)
   x[is.na.x] <- y[is.na.x]

   return(x)

}
