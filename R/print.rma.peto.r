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
      fs <- fmtx(x$fit.stats$ML, digits[["fit"]])
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      #tmp[1] <- paste0(tmp[1], "\u200b")
      .print.table(tmp, mstyle)
   }

   cat("\n")

   if (!is.na(x$I2)) {
      cat(mstyle$text("I^2 (total heterogeneity / total variability):  "))
      cat(mstyle$result(paste0(fmtx(x$I2, 2), "%")))
      cat("\n")
   }
   if (!is.na(x$H2)) {
      cat(mstyle$text("H^2 (total variability / sampling variability): "))
      cat(mstyle$result(fmtx(x$H2, 2)))
      cat("\n")
   }

   if (!is.na(x$QE)) {
      cat("\n")
      cat(mstyle$section("Test for Heterogeneity:"), "\n")
      cat(mstyle$result(fmtt(x$QE, "Q", df=x$k.pos-1, pval=x$QEp, digits=digits)))
   }

   if (any(!is.na(c(x$I2, x$H2, x$QE))))
      cat("\n\n")

   res.table <- c(estimate=fmtx(unname(x$beta), digits[["est"]]), se=fmtx(x$se, digits[["se"]]), zval=fmtx(x$zval, digits[["test"]]), pval=fmtp(x$pval, digits[["pval"]]), ci.lb=fmtx(x$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$ci.ub, digits[["ci"]]))
   res.table.exp <- c(estimate=fmtx(exp(unname(x$beta)), digits[["est"]]), ci.lb=fmtx(exp(x$ci.lb), digits[["ci"]]), ci.ub=fmtx(exp(x$ci.ub), digits[["ci"]]))

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
