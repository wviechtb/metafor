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
