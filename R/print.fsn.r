print.fsn <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="fsn")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   cat(mstyle$section(paste("Fail-safe N Calculation Using the", x$type, "Approach")))
   cat("\n\n")

   if (x$type == "Rosenthal" || x$type == "Binomial") {
      cat(mstyle$text("Observed Significance Level: "))
      cat(mstyle$result(fmtp(x$pval, digits[["pval"]])))
      cat("\n")
      cat(mstyle$text("Target Significance Level:   "))
      cat(mstyle$result(round(x$alpha, digits[["pval"]])))
   }

   if (x$type == "Orwin") {
      cat(mstyle$text("Average Effect Size: "))
      cat(mstyle$result(fmtx(x$est, digits[["est"]])))
      cat("\n")
      cat(mstyle$text("Target Effect Size:  "))
      cat(mstyle$result(fmtx(x$target, digits[["est"]])))
   }

   if (x$type == "Rosenberg") {
      flag.left  <- ifelse(isTRUE(x$est < 0), " ", "")
      cat(mstyle$text("Average Effect Size:         "))
      cat(mstyle$result(fmtx(x$est, digits[["est"]], flag=flag.left)))
      cat("\n")
      cat(mstyle$text("Observed Significance Level: "))
      cat(flag.left)
      cat(mstyle$result(fmtp(x$pval, digits[["pval"]])))
      cat("\n")
      cat(mstyle$text("Target Significance Level:   "))
      cat(flag.left)
      cat(mstyle$result(round(x$alpha, digits[["pval"]])))
   }

   if (x$type == "General") {
      flag.left  <- ifelse(isTRUE(x$est < 0), " ", "")
      flag.right <- ifelse(isTRUE(x$est.fsn < 0), " ", "")
      cat(mstyle$text("Average Effect Size:         "))
      cat(mstyle$result(fmtx(x$est, digits[["est"]], flag=flag.left)))
      if (x$fsnum > 0) {
         cat(mstyle$text(" (with file drawer: "))
         cat(mstyle$result(fmtx(x$est.fsn, digits[["est"]], flag=flag.right)))
         cat(mstyle$text(")"))
      }
      cat("\n")
      if (!is.element(x$method, c("FE","EE","CE"))) {
         cat(mstyle$text("Amount of Heterogeneity:     "))
         cat(mstyle$result(fmtx(x$tau2, digits[["var"]], flag=flag.left)))
         if (x$fsnum > 0) {
            cat(mstyle$text(" (with file drawer: "))
            cat(mstyle$result(fmtx(x$tau2.fsn, digits[["var"]], flag=flag.right)))
            cat(mstyle$text(")"))
         }
         cat("\n")
      }
      cat(mstyle$text("Observed Significance Level: "))
      cat(flag.left)
      cat(mstyle$result(fmtp(x$pval, digits[["pval"]])))
      if (x$fsnum > 0) {
         cat(mstyle$text(" (with file drawer: "))
         cat(flag.right)
         cat(mstyle$result(fmtp(x$pval.fsn, digits[["pval"]])))
         cat(mstyle$text(")"))
      }
      cat("\n")
      if (is.na(x$target)) {
         cat(mstyle$text("Target Significance Level:   "))
         cat(flag.left)
         cat(mstyle$result(round(x$alpha, digits[["pval"]])))
      } else {
         cat(mstyle$text("Target Effect Size:          "))
         cat(mstyle$result(fmtx(x$target, digits[["est"]], , flag=flag.left)))
      }
   }

   cat("\n\n")
   cat(mstyle$text("Fail-safe N: "))
   cat(mstyle$result(paste0(x$ub.sign, x$fsnum)))
   cat("\n")

   .space()

   invisible()

}
