print.confint.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="confint.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   if (names(x)[1] == "fixed") {

      res.fixed <- cbind(fmtx(x$fixed[,1,drop=FALSE], digits[["est"]]), fmtx(x$fixed[,2:3,drop=FALSE], digits[["ci"]]))
      tmp <- capture.output(print(res.fixed, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

   }

   if (is.element("random", names(x))) {

      if (names(x)[1] == "fixed")
         cat("\n")

      res.random <- fmtx(x$random, digits[["var"]])
      res.random[,2] <- paste0(x$lb.sign, res.random[,2])
      res.random[,3] <- paste0(x$ub.sign, res.random[,3])
      tmp <- capture.output(print(res.random, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

      # this can only (currently) happen for 'rma.uni' models

      if (x$ci.null)
         message(mstyle$message(paste0("\nThe upper and lower CI bounds for tau^2 both fall below ", round(x$tau2.min,4), ".\nThe CIs are therefore equal to the null/empty set.")))

   }

   .space()

   invisible()

}
