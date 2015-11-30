print.confint.rma <- function(x, digits, ...) {

   if (!is.element("confint.rma", class(x)))
      stop("Argument 'x' must be an object of class \"confint.rma\".")

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   if (names(x)[1] == "fixed") {

      res.fixed <- formatC(x$fixed, digits=digits, format="f")
      print(res.fixed, quote=FALSE, right=TRUE)

   }

   if (is.element("random", names(x))) {

      if (names(x)[1] == "fixed")
         cat("\n")

      res.random <- formatC(x$random, digits=digits, format="f")
      res.random[,2] <- paste0(x$lb.sign, res.random[,2])
      res.random[,3] <- paste0(x$ub.sign, res.random[,3])
      print(res.random, quote=FALSE, right=TRUE)

      ### this can only (currently) happen for 'rma.uni' models

      if (x$ci.null)
         message("\nThe upper and lower CI bounds for tau^2 both fall below ", x$tau2.min, ".\nThe CIs are therefore equal to the null/empty set.", sep="")

   }

   cat("\n")
   invisible()

}
