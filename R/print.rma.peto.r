print.rma.peto <- function(x, digits, showfit=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.peto"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.peto\"."))

   if (missing(digits))
      digits <- x$digits

   if (!exists(".rmspace"))
      cat("\n")

   cat(mstyle$section("Fixed-Effects Model"))
   cat(mstyle$section(paste0(" (k = ", x$k, ")")))

   if (showfit) {
      cat("\n")
      if (anyNA(x$fit.stats$ML)) {
         fs <- x$fit.stats$ML
      } else {
         fs <- c(formatC(x$fit.stats$ML, digits=digits, format="f"))
      }
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      .print.table(tmp, mstyle)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (!is.na(x$QE)) {
      cat(mstyle$section("Test for Heterogeneity:"), "\n")
      cat(mstyle$result(paste0("Q(df = ", x$k.pos-1, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "))))
   }

   res.table     <- c(x$beta, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
   res.table.exp <- c(exp(x$beta), exp(x$ci.lb), exp(x$ci.ub))

   if (!is.na(x$beta)) {
      res.table    <- formatC(res.table, digits=digits, format="f")
      res.table[4] <- .pval(x$pval, digits=digits)
   }

   if (!is.na(x$beta))
      res.table.exp <- formatC(res.table.exp, digits=digits, format="f")

   names(res.table)     <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
   names(res.table.exp) <- c("estimate", "ci.lb", "ci.ub")

   cat("\n\n")
   cat(mstyle$section("Model Results (log scale):"))
   cat("\n\n")
   tmp <- capture.output(.print.vector(res.table))
   .print.table(tmp, mstyle)

   cat("\n")
   cat(mstyle$section("Model Results (OR scale):"))
   cat("\n\n")
   tmp <- capture.output(.print.vector(res.table.exp))
   .print.table(tmp, mstyle)

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
