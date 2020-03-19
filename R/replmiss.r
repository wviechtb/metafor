replmiss <- function(x, y) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (length(y) == 0L)
      y <- NA

   ### catch cases where x is of length 0

   if (length(x) == 0L)
      return(x)

   ### in case user specifies a constant to use for replacement

   if (length(y) == 1L)
      y <- rep(y, length(x))

   ### check that x and y are of the same length

   if (length(x) != length(y))
      stop(mstyle$stop("Length of 'x' and 'y' is not the same."))

   #x <- ifelse(is.na(x), y, x) # this is quite a bit slower than the following
   x[is.na(x)] <- y[is.na(x)]

   return(x)

}
