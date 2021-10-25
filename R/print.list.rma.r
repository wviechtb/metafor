print.list.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="list.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   attr(x, "class") <- NULL

   ### remove cr.lb and cr.ub elements (if they are there)

   x$cr.lb <- NULL
   x$cr.ub <- NULL

   ### turn all vectors before the slab vector into a data frame

   slab.pos <- which(names(x) == "slab")
   out <- x[seq_len(slab.pos-1)]
   out <- data.frame(out, row.names=x$slab, stringsAsFactors=FALSE)

   ### in case all values were NA and have been omitted

   if (nrow(out) == 0L)
      stop(mstyle$stop("All values are NA."), call.=FALSE)

   ### in case there is a select element, apply it

   if (exists("select", where=x, inherits=FALSE))
      out <- out[x$select,]

   if (nrow(out) == 0L) {
      message(mstyle$message("No values to print."))
      return(invisible())
   }

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
      min.pos <- slab.pos - is.element("tau2.level", names(x)) - is.element("gamma2.level", names(x)) - is.element("X", names(x)) - is.element("Z", names(x)) - transf.true
   } else {
      min.pos <- slab.pos - transf.true
   }

   sav <- out[,seq_len(min.pos-1)]

   for (i in seq_len(min.pos-1)) {
      if (inherits(out[,i], c("integer","logical","factor","character"))) { ### do not apply formating to these classes
         out[,i] <- out[,i]
      } else {
         if (names(out)[i] %in% c("pred", "resid"))
            out[,i] <- .fcf(out[,i], digits[["est"]])
         if (names(out)[i] %in% c("se"))
            out[,i] <- .fcf(out[,i], digits[["se"]])
         if (names(out)[i] %in% c("ci.lb", "ci.ub", "cr.lb", "cr.ub", "pi.lb", "pi.ub"))
            out[,i] <- .fcf(out[,i], digits[["ci"]])
         if (names(out)[i] %in% c("zval", "Q", "z", "X2"))
            out[,i] <- .fcf(out[,i], digits[["test"]])
         if (names(out)[i] %in% c("pval", "Qp"))
            out[,i] <- .fcf(out[,i], digits[["pval"]])
         if (names(out)[i] %in% c("I2", "H2"))
            out[,i] <- .fcf(out[,i], digits[["het"]])
         if (names(out)[i] %in% c("tau2"))
            out[,i] <- .fcf(out[,i], digits[["var"]])
         # if (names(out)[i] == "rstudent")
         #    out[,i] <- .fcf(out[,i], digits[["test"]])
         # if (names(out)[i] == "dffits")
         #    out[,i] <- .fcf(out[,i], digits[["test"]])
         # if (names(out)[i] == "cook.d")
         #    out[,i] <- .fcf(out[,i], digits[["test"]])
         # if (names(out)[i] == "cov.r")
         #    out[,i] <- .fcf(out[,i], digits[["test"]])
         # if (names(out)[i] == "tau2.del")
         #    out[,i] <- .fcf(out[,i], digits[["var"]])
         # if (names(out)[i] == "QE.del")
         #    out[,i] <- .fcf(out[,i], digits[["test"]])
         # if (names(out)[i] == "hat")
         #    out[,i] <- .fcf(out[,i], digits[["test"]])
         # if (names(out)[i] == "weight")
         #    out[,i] <- .fcf(out[,i], digits[["test"]])
         # if (names(out)[i] == "dfbs")
         #    out[,i] <- .fcf(out[,i], digits[["est"]])
         if (!is.character(out[,i]))
            out[,i] <- .fcf(out[,i], digits[["est"]])
      }
   }

   .space()

   tmp <- capture.output(print(out, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   if (is.null(attr(x, ".rmspace"))) .space()

   invisible(sav)

}
