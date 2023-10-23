vec2mat <- function(x, diag=FALSE, corr=!diag, dimnames) {

   mstyle <- .get.mstyle()

   p <- length(x)

   dims <- sqrt(2*p + 1/4) + ifelse(diag, -1/2, 1/2)

   if (abs(dims - round(dims)) >= .Machine$double.eps^0.5)
      stop(mstyle$stop("Length of 'x' does not correspond to a square matrix."))

   dims <- round(dims)

   R <- matrix(NA_real_, nrow=dims, ncol=dims)

   if (!missing(dimnames)) {
      if (length(dimnames) != dims)
         stop(mstyle$stop(paste0("Length of 'dimnames' (", length(dimnames), ") does not correspond to the dimensions of the matrix (", dims, ").")))
      rownames(R) <- colnames(R) <- dimnames
   }

   R[lower.tri(R, diag=diag)] <- x
   R[upper.tri(R, diag=diag)] <- t(R)[upper.tri(R, diag=diag)]

   if (corr)
      diag(R) <- 1

   return(R)

}
