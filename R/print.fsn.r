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
      cat("Observed Significance Level:", .pval(x$pval, digits=digits), "\n")
      cat("Target Significance Level:  ", x$alpha, "\n\n")
   }

   if (x$type == "Orwin") {
      cat("Average Effect Size:", formatC(x$meanes, digits=digits, format="f"), "\n")
      cat("Target Effect Size: ", formatC(x$target, digits=digits, format="f"), "\n\n")
   }

   if (x$type == "Rosenberg") {
      cat("Average Effect Size:        ", formatC(x$meanes, digits=digits, format="f"), "\n")
      cat("Observed Significance Level:", .pval(x$pval, digits=digits), "\n")
      cat("Target Significance Level:  ", x$alpha, "\n\n")
   }

   cat("Fail-safe N:", x$fsnum, "\n\n")

   invisible()

}
