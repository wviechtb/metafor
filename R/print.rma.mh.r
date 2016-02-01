print.rma.mh <- function(x, digits, showfit=FALSE, ...) {

   if (!inherits(x, "rma.mh"))
      stop("Argument 'x' must be an object of class \"rma.mh\".")

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   cat("Fixed-Effects Model (k = ", x$k, ")", sep="")

   if (showfit) {
      cat("\n")
      if (anyNA(x$fit.stats$ML)) {
         fs <- x$fit.stats$ML
      } else {
         fs <- c(formatC(x$fit.stats$ML, digits=digits, format="f"))
      }
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      print(fs, quote=FALSE, print.gap=2)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (!is.na(x$QE)) {
      cat("Test for Heterogeneity: \n")
      cat("Q(df = ", ifelse(x$k.yi-1 >= 0, x$k.yi-1, 0), ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "), sep="")
   }

   if (x$measure == "OR" || x$measure == "RR" || x$measure == "IRR") {

      res.table     <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
      res.table.exp <- c(exp(x$b), exp(x$ci.lb), exp(x$ci.ub))

      if (!is.na(x$b)) {
         res.table    <- formatC(res.table, digits=digits, format="f")
         res.table[4] <- .pval(x$pval, digits=digits)
      }

      if (!is.na(x$b))
         res.table.exp <- formatC(res.table.exp, digits=digits, format="f")

      names(res.table)     <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
      names(res.table.exp) <- c("estimate", "ci.lb", "ci.ub")

      cat("\n\nModel Results (log scale):\n\n")
      print(res.table, quote=FALSE, right=TRUE)

      cat("\nModel Results (", x$measure, " scale):", "\n\n", sep="")
      print(res.table.exp, quote=FALSE, right=TRUE)
      cat("\n")

      if (x$measure == "OR") {
         MH <- ifelse(is.na(x$MH), NA, formatC(x$MH, digits=digits, format="f"))
         TA <- ifelse(is.na(x$TA), NA, formatC(x$TA, digits=digits, format="f"))
         width <- max(nchar(MH), nchar(TA))
         if (is.na(MH)) {
            cat("Cochran-Mantel-Haenszel Test:    test value not computable for these data \n", sep="")
         } else {
            cat("Cochran-Mantel-Haenszel Test:    CMH = ", formatC(MH, width=width), ", df = 1,", paste(rep(" ", nchar(x$k.pos)-1, collapse="")), " p-val ", .pval(x$MHp, digits=digits, showeq=TRUE, sep=" ", add0=TRUE), "\n", sep="")
         }
         if (is.na(TA)) {
            cat("Tarone's Test for Heterogeneity: test value not computable for these data \n\n", sep="")
         } else {
            cat("Tarone's Test for Heterogeneity: X^2 = ", formatC(TA, width=width), ", df = ", x$k.pos-1, ", p-val ", .pval(x$TAp, digits=digits, showeq=TRUE, sep=" ", add0=TRUE), "\n\n", sep="")
         }
      }

      if (x$measure == "IRR") {
         if (is.na(x$MH)) {
            cat("Mantel-Haenszel Test: test value not computable for these data \n", sep="")
         } else {
            cat("Mantel-Haenszel Test: MH = ", formatC(x$MH, digits, format="f"), ", df = 1, p-val ", .pval(x$MHp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
         }
      }

   } else {

      res.table <- c(x$b, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)

      if (!is.na(x$b)) {
         res.table    <- formatC(res.table, digits=digits, format="f")
         res.table[4] <- .pval(x$pval, digits=digits)
      }

      names(res.table) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

      cat("\n\nModel Results:\n\n")
      print(res.table, quote=FALSE, right=TRUE)
      cat("\n")

   }

   invisible()

}
