print.rma.peto <- function(x, digits, showfit=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.peto")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   .space()

   cat(mstyle$section("Equal-Effects Model"))
   cat(mstyle$section(paste0(" (k = ", x$k, ")")))

   cat("\n")

   if (showfit) {
      fs <- .fcf(x$fit.stats$ML, digits[["fit"]])
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      .print.table(tmp, mstyle)
   }

   cat("\n")

   if (!is.na(x$I2)) {
      cat(mstyle$text("I^2 (total heterogeneity / total variability):  "))
      cat(mstyle$result(paste0(ifelse(is.na(x$I2), NA, .fcf(x$I2, 2)), "%")))
      cat("\n")
   }
   if (!is.na(x$H2)) {
      cat(mstyle$text("H^2 (total variability / sampling variability): "))
      cat(mstyle$result(paste0(ifelse(is.na(x$H2), NA, .fcf(x$H2, 2)))))
      cat("\n")
   }

   if (!is.na(x$QE)) {
      cat("\n")
      cat(mstyle$section("Test for Heterogeneity:"), "\n")
      cat(mstyle$result(paste0("Q(df = ", x$k.pos-1, ") = ", .fcf(x$QE, digits[["test"]]), ", p-val ", .pval(x$QEp, digits[["pval"]], showeq=TRUE, sep=" "))))
   }

   if (any(!is.na(c(x$I2, x$H2, x$QE))))
      cat("\n\n")

   res.table <- c(estimate=.fcf(unname(x$beta), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), pval=.pval(x$pval, digits[["pval"]]), ci.lb=.fcf(x$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$ci.ub, digits[["ci"]]))
   res.table.exp <- c(estimate=.fcf(exp(unname(x$beta)), digits[["est"]]), ci.lb=.fcf(exp(x$ci.lb), digits[["ci"]]), ci.ub=.fcf(exp(x$ci.ub), digits[["ci"]]))

   cat(mstyle$section("Model Results (log scale):"))
   cat("\n\n")
   tmp <- capture.output(.print.vector(res.table))
   .print.table(tmp, mstyle)

   cat("\n")
   cat(mstyle$section("Model Results (OR scale):"))
   cat("\n\n")
   tmp <- capture.output(.print.vector(res.table.exp))
   .print.table(tmp, mstyle)

   .space()

   invisible()

}
