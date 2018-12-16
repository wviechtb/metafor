print.confint.rma <- function(x, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "confint.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"confint.rma\"."))

   if (missing(digits))
      digits <- x$digits

   if (!exists(".rmspace"))
      cat("\n")

   if (names(x)[1] == "fixed") {

      res.fixed <- formatC(x$fixed, digits=digits, format="f")
      tmp <- capture.output(print(res.fixed, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

   }

   if (is.element("random", names(x))) {

      if (names(x)[1] == "fixed")
         cat("\n")

      res.random <- formatC(x$random, digits=digits, format="f")
      res.random[,2] <- paste0(x$lb.sign, res.random[,2])
      res.random[,3] <- paste0(x$ub.sign, res.random[,3])
      tmp <- capture.output(print(res.random, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

      ### this can only (currently) happen for 'rma.uni' models

      if (x$ci.null)
         message(mstyle$message(paste0("\nThe upper and lower CI bounds for tau^2 both fall below ", x$tau2.min, ".\nThe CIs are therefore equal to the null/empty set.")))

   }

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
