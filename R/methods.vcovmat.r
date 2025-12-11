print.vcovmat <- function(x, digits=4, tol, zero=".", na="", ...) {

   d <- dim(x)
   if (any(d == 0)) {
      cat("< table of extent", paste(d, collapse = " x "), ">\n")
      return(invisible(x))
   }

   if (missing(tol))
      tol <- 10 * .Machine$double.eps

   xx <- formatC(unclass(x), format="f", digits=digits, flag=" ")

   if (any(ina <- is.na(x)))
      xx[ina] <- na

   if (zero != "0" && any(i0 <- !ina & abs(x) <= tol))
      xx[i0] <- zero

   print(xx, quote=FALSE, ...)

   invisible(x)

}

"[.vcovmat" <- function(x, i, j, ...) {

   out <- NextMethod("[")

   if (inherits(out, "matrix"))
      class(out) <- class(x)

   return(out)

}
