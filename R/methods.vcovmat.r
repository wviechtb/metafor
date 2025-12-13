print.vcovmat <- function(x, digits=4, tol, zero=".", na="NA", ...) {

   mstyle <- .get.mstyle()

   d <- dim(x)

   if (any(d == 0)) {
      cat("< table of extent", paste(d, collapse = " x "), ">\n")
      return(invisible(x))
   }

   if (missing(tol))
      tol <- 10 * .Machine$double.eps

   xx <- formatC(unclass(x), format="f", digits=digits, flag=if (any(x < 0, na.rm=TRUE)) " " else "")

   if (any(ina <- is.na(x)))
      xx[ina] <- na

   if (zero != "0" && any(i0 <- !ina & abs(x) <= tol))
      xx[i0] <- zero

   if (is.null(colnames(xx)))
      colnames(xx) <- 1:ncol(xx)

   if (is.null(rownames(xx)))
      rownames(xx) <- 1:ncol(xx)

   #print(xx, quote=FALSE, right=TRUE, ...)

   .space()

   tmp <- capture.output(print(xx, quote=FALSE, right=TRUE, ...))
   .print.vcovmat(tmp, mstyle)

   .space()

   invisible()

}

"[.vcovmat" <- function(x, i, j, ...) {

   out <- NextMethod("[")

   if (inherits(out, "matrix"))
      class(out) <- class(x)

   return(out)

}
