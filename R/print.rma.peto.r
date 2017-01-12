print.rma.peto <- function(x, digits, showfit=FALSE, ...) {

   if (!inherits(x, "rma.peto"))
      stop("Argument 'x' must be an object of class \"rma.peto\".")

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
      cat("Q(df = ", x$k.pos-1, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "), sep="")
   }

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
   .print.out(res.table)
   #print(res.table, quote=FALSE, right=TRUE)

   cat("\nModel Results (OR scale):\n\n")
   .print.out(res.table.exp)
   #print(res.table.exp, quote=FALSE, right=TRUE)
   cat("\n")

   invisible()

}
