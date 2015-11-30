print.list.rma <- function(x, digits, ...) {

   if (!is.element("list.rma", class(x)))
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

   out[,seq_len(min.pos-1)] <- apply(out[,seq_len(min.pos-1), drop=FALSE], 2, formatC, digits=digits, format="f")

   print(out, quote=FALSE, right=TRUE)

   invisible(sav)

}
