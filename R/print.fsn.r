print.fsn <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="fsn")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   cat(mstyle$section(paste("Fail-safe N Calculation Using the", x$type, "Approach")))
   cat("\n\n")

   if (x$type == "Rosenthal") {
      cat(mstyle$text("Observed Significance Level: "))
      cat(mstyle$result(.pval(x$pval, digits[["pval"]])))
      cat("\n")
      cat(mstyle$text("Target Significance Level:   "))
      cat(mstyle$result(x$alpha))
   }

   if (x$type == "Orwin") {
      cat(mstyle$text("Average Effect Size: "))
      cat(mstyle$result(.fcf(x$meanes, digits[["est"]])))
      cat("\n")
      cat(mstyle$text("Target Effect Size:  "))
      cat(mstyle$result(.fcf(x$target, digits[["est"]])))
   }

   if (x$type == "Rosenberg") {
      cat(mstyle$text("Average Effect Size:         "))
      cat(mstyle$result(.fcf(x$meanes, digits[["est"]])))
      cat("\n")
      cat(mstyle$text("Observed Significance Level: "))
      cat(mstyle$result(.pval(x$pval, digits[["pval"]])))
      cat("\n")
      cat(mstyle$text("Target Significance Level:   "))
      cat(mstyle$result(x$alpha))
   }

   if (x$type == "REM") {
      cat(mstyle$text("Average Effect Size: "))
      cat(mstyle$result(.fcf(x$meanes, digits[["est"]])))
      cat("\n")
      cat(mstyle$text("Target Effect Size:  "))
      cat(mstyle$result(.fcf(x$target, digits[["est"]])))
      #cat("\n")
      #cat(mstyle$text("Rejection Rate:      "))
      #cat(mstyle$result(.pval(x$rejrate, digits[["pval"]])))
   }

   cat("\n\n")
   cat(mstyle$text("Fail-safe N: "))
   cat(mstyle$result(x$fsnum))
   cat("\n")

   .space()

   invisible()

}
