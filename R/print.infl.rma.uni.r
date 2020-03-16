print.infl.rma.uni <- function(x, digits=x$digits, infonly=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "infl.rma.uni"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"infl.rma.uni\"."))

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   if (x$p == 1) {

      out <- list(rstudent=x$inf$rstudent, dffits=x$inf$dffits, cook.d=x$inf$cook.d, cov.r=x$inf$cov.r,
                tau2.del=x$inf$tau2.del, QE.del=x$inf$QE.del, hat=x$inf$hat, weight=x$inf$weight,
                dfbs=x$dfbs[[1]], inf=x$inf$inf, slab=x$inf$slab, digits=digits)
      class(out) <- "list.rma"

      if (infonly)
         out[["select"]] <- !is.na(x$is.infl) & x$is.infl

   } else {

      out <- x[1:2]
      out$inf[["digits"]]  <- digits
      out$dfbs[["digits"]] <- digits
      attr(out$inf,  ".rmspace") <- TRUE
      attr(out$dfbs, ".rmspace") <- TRUE

      if (infonly) {
         out$inf[["select"]] <- !is.na(x$is.infl) & x$is.infl
         out$dfbs[["select"]] <- !is.na(x$is.infl) & x$is.infl
      }

   }

   print(out)

}
