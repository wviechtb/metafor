print.fsn <- function(x, digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "fsn"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"fsn\"."))

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   cat(mstyle$section(paste("Fail-safe N Calculation Using the", x$type, "Approach")))
   cat("\n\n")

   if (x$type == "Rosenthal") {
      cat(mstyle$text("Observed Significance Level: "))
      cat(mstyle$result(.pval(x$pval, digits=digits)))
      cat("\n")
      cat(mstyle$text("Target Significance Level:   "))
      cat(mstyle$result(x$alpha))
   }

   if (x$type == "Orwin") {
      cat(mstyle$text("Average Effect Size: "))
      cat(mstyle$result(formatC(x$meanes, digits=digits, format="f")))
      cat("\n")
      cat(mstyle$text("Target Effect Size:  "))
      cat(mstyle$result(formatC(x$target, digits=digits, format="f")))
   }

   if (x$type == "Rosenberg") {
      cat(mstyle$text("Average Effect Size:        "))
      cat(mstyle$result(formatC(x$meanes, digits=digits, format="f")))
      cat("\n")
      cat(mstyle$text("Observed Significance Level: "))
      cat(mstyle$result(.pval(x$pval, digits=digits)))
      cat("\n")
      cat(mstyle$text("Target Significance Level:   "))
      cat(mstyle$result(x$alpha))
   }

   cat("\n\n")
   cat(mstyle$text("Fail-safe N: "))
   cat(mstyle$result(x$fsnum))
   cat("\n\n")

   invisible()

}
