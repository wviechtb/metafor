############################################################################

"[.list.rma" <- function(x, i, ...) { ### removed j argument (see below), so can only select rows, not columns

   out <- x

   attr(out, "class") <- NULL

   slab.pos <- which(names(out) == "slab")

   if (!missing(i)) ### for X element
      out[seq_len(slab.pos-1)] <- lapply(out[seq_len(slab.pos-1)], function(r) if (class(r) == "matrix") r[i,] else r[i])

   ### catch cases where user selects values outside 1:k

   if (length(out[[1]]) == 0L)
      return(NULL)

   #out <- out[j] ### this causes all kinds of problems, so left out for now

   out$slab <- x$slab[i]

   ### slab can only contain NAs if user selects values outside 1:k

   if (anyNA(out$slab))
      return(NULL)

   out$digits <- x$digits
   out$transf <- x$transf
   out$method <- x$method

   class(out) <- "list.rma"
   return(out)

}

############################################################################

as.data.frame.list.rma <- function(x, ...) {

   attr(x, "class") <- NULL

   ### turn all vectors before the slab vector into a data frame

   slab.pos <- which(names(x) == "slab")
   out <- x[seq_len(slab.pos-1)]
   out <- data.frame(out, row.names=x$slab)

   ### in case all values were NA and have been omitted

   if (nrow(out) == 0L)
      return(data.frame())

   ### if transf exists and is TRUE, set SEs to NULL so that column is omitted from the output

   if (exists("transf", where=x, inherits=FALSE) && x$transf)
      out$se <- NULL

   return(out)

}

############################################################################

as.matrix.list.rma <- function(x, ...) {

   attr(x, "class") <- NULL

   ### turn all vectors before the slab vector into a matrix

   slab.pos <- which(names(x) == "slab")
   out <- x[seq_len(slab.pos-1)]
   out <- do.call(cbind, out)
   rownames(out) <- x$slab

   ### if transf exists and is TRUE, set SEs to NULL so that column is omitted from the output

   if (exists("transf", where=x, inherits=FALSE) && x$transf)
      out <- out[,-which(colnames(out) == "se")]

   return(out)

}

############################################################################

### like utils:::head.data.frame and utils:::tail.data.frame,
### but with nrow(x) replaced by length(x[[1]])

head.list.rma <- function (x, n = 6L, ...) {

   stopifnot(length(n) == 1L)

   n <- if (n < 0L)  {
      max(length(x[[1]]) + n, 0L)
   } else {
      min(n, length(x[[1]]))
   }

   x[seq_len(n), , drop = FALSE]

}

tail.list.rma <- function (x, n = 6L, ...) {

   stopifnot(length(n) == 1L)

   nrx <- length(x[[1]])

   n <- if (n < 0L) {
      max(nrx + n, 0L)
   } else {
      min(n, nrx)
   }

   x[seq.int(to = nrx, length.out = n), , drop = FALSE]

}

############################################################################

`$<-.list.rma` <- function(x, name, value) {

   slab.pos <- which(names(x) == "slab")

   out <- list()

   for (i in seq_len(slab.pos-1)) {
      out[[i]] <- x[[i]]
   }

   names(out) <- names(x)[seq_len(slab.pos-1)]

   out[[name]] <- value

   for (i in (slab.pos:length(x))) {
      out[[i+1]] <- x[[i]]
   }

   names(out)[(slab.pos+1):(length(x)+1)] <- names(x)[slab.pos:length(x)]

   class(out) <- class(x)
   return(out)

}

############################################################################
