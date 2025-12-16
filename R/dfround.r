dfround <- function(x, digits, drop0=TRUE) {

   mstyle <- .get.mstyle()

   if (inherits(x, "matrix") && length(dim(x)) == 2L)
      x <- data.frame(x, check.names=FALSE)

   .chkclass(class(x), must="data.frame")

   p <- ncol(x)

   if (missing(digits))
      digits <- 0

   digits <- .expand1(digits, p)
   drop0  <- .expand1(drop0, p)

   if (p != length(digits))
      stop(mstyle$stop(paste0("Number of columns in 'x' (", p, ") does not match the length of 'digits' (", length(digits), ").")))

   if (p != length(drop0))
      stop(mstyle$stop(paste0("Number of columns in 'x' (", p, ") does not match the length of 'drop0' (", length(drop0), ").")))

   if (!is.numeric(digits))
      stop(mstyle$stop("Argument 'digits' must be a numeric vector."))

   if (!is.logical(drop0))
      stop(mstyle$stop("Argument 'drop0' must be a logical vector."))

   for (i in seq_len(p)) {
      if (!is.numeric(x[[i]]))
         next
      if (drop0[i]) {
         x[[i]] <- round(x[[i]], digits[i])
      } else {
         x[[i]] <- formatC(x[[i]], format="f", digits=digits[i])
      }
   }

   return(x)

}
