############################################################################

print.list.rma <- function(x, digits, ...) {

   if (!inherits(x, "list.rma"))
      stop("Argument 'x' must be an object of class \"list.rma\".")

   if (missing(digits))
      digits <- x$digits

   attr(x, "class") <- NULL

   ### turn all vectors before the slab vector into a data frame

   slab.pos <- which(names(x) == "slab")
   out <- x[seq_len(slab.pos-1)]
   out <- data.frame(out, row.names=x$slab)

   ### in case all values were NA and have been omitted

   if (nrow(out) == 0L)
      stop("All values are NA.", call.=FALSE)

   ### if transf exists and is TRUE, set SEs to NULL so that column is omitted from the output

   transf.true <- 0

   if (exists("transf", where=x, inherits=FALSE) && x$transf) {
      transf.true <- 1
      out$se <- NULL
   }

   ### objects created by predict.rma() have a 'method' element
   ### properly format columns 1-4 (for FE models) or columns 1-6 (for RE/ME models)
   ### leave element tau2.level, gamma2.level, and/or element X untouched

   if (exists("method", where=x, inherits=FALSE)) {
      min.pos <- slab.pos - is.element("tau2.level", names(x)) - is.element("gamma2.level", names(x)) - is.element("X", names(x)) - transf.true
   } else {
      min.pos <- slab.pos - transf.true
   }

   sav <- out[,seq_len(min.pos-1)]

   for (i in 1:(min.pos-1)) {
      if (inherits(out[,i], "integer")) { ### do not apply formating to integers
         out[,i] <- out[,i]
      } else {
         out[,i] <- formatC(out[,i], digits=digits, format="f")
      }
   }

   print(out, quote=FALSE, right=TRUE)

   invisible(sav)

}

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

   if (!inherits(x, "list.rma"))
      stop("Argument 'x' must be an object of class \"list.rma\".")

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

   for (i in 1:(slab.pos-1)) {
      out[[i]] <- x[[i]]
   }

   names(out) <- names(x)[1:(slab.pos-1)]

   out[[name]] <- value

   for (i in (slab.pos:length(x))) {
      out[[i+1]] <- x[[i]]
   }

   names(out)[(slab.pos+1):(length(x)+1)] <- names(x)[slab.pos:length(x)]

   class(out) <- class(x)
   return(out)

}

############################################################################
