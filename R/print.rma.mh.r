print.rma.mh <- function(x, digits, showfit=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.mh")

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
      cat(mstyle$result(fmtt(x$QE, "Q", df=ifelse(x$k.yi-1 >= 0, x$k.yi-1, 0), pval=x$QEp, digits=digits)))
   }

   if (any(!is.na(c(x$I2, x$H2, x$QE))))
      cat("\n\n")

   if (is.element(x$measure, c("OR","RR","IRR"))) {

      res.table <- c(estimate=fmtx(unname(x$beta), digits[["est"]]), se=fmtx(x$se, digits[["se"]]), zval=fmtx(x$zval, digits[["test"]]), pval=fmtp(x$pval, digits[["pval"]]), ci.lb=fmtx(x$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$ci.ub, digits[["ci"]]))
      res.table.exp <- c(estimate=fmtx(exp(unname(x$beta)), digits[["est"]]), ci.lb=fmtx(exp(x$ci.lb), digits[["ci"]]), ci.ub=fmtx(exp(x$ci.ub), digits[["ci"]]))

      cat(mstyle$section("Model Results (log scale):"))
      cat("\n\n")
      tmp <- capture.output(.print.vector(res.table))
      .print.table(tmp, mstyle)

      cat("\n")
      cat(mstyle$section(paste0("Model Results (", x$measure, " scale):")))
      cat("\n\n")
      tmp <- capture.output(.print.vector(res.table.exp))
      .print.table(tmp, mstyle)

      if (x$measure == "OR") {
         cat("\n")
         MH <- fmtx(x$MH, digits[["test"]])
         TA <- fmtx(x$TA, digits[["test"]])
         if (is.na(MH) && is.na(TA)) {
            width <- 1
         } else {
            width <- max(nchar(MH), nchar(TA), na.rm=TRUE)
         }
         cat(mstyle$text("Cochran-Mantel-Haenszel Test:    "))
         if (is.na(MH)) {
            cat(mstyle$result("test value not computable for these data"))
            cat("\n")
         } else {
            cat(mstyle$result(paste0("CMH = ", formatC(MH, width=width), ", df = 1,", paste(rep(" ", nchar(x$k.pos)-1L), collapse=""), " p-val ", fmtp(x$MHp, digits[["pval"]], equal=TRUE, sep=TRUE, add0=TRUE))))
            cat("\n")
         }
         cat(mstyle$text("Tarone's Test for Heterogeneity: "))
         if (is.na(TA)) {
            cat(mstyle$result("test value not computable for these data"))
         } else {
            cat(mstyle$result(paste0("X^2 = ", formatC(TA, width=width), ", df = ", x$k.pos-1, ", p-val ", fmtp(x$TAp, digits[["pval"]], equal=TRUE, sep=TRUE, add0=TRUE))))
         }
         cat("\n")
      }

      if (x$measure == "IRR") {
         cat("\n")
         cat(mstyle$text("Mantel-Haenszel Test: "))
         if (is.na(x$MH)) {
            cat(mstyle$result("test value not computable for these data"))
         } else {
            cat(mstyle$result(paste0("MH = ", fmtx(x$MH, digits[["test"]]), ", df = 1, p-val ", fmtp(x$MHp, digits[["pval"]], equal=TRUE, sep=TRUE))))
         }
         cat("\n")
      }

   } else {

      res.table <- c(estimate=fmtx(unname(x$beta), digits[["est"]]), se=fmtx(x$se, digits[["se"]]), zval=fmtx(x$zval, digits[["test"]]), pval=fmtp(x$pval, digits[["pval"]]), ci.lb=fmtx(x$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$ci.ub, digits[["ci"]]))

      cat(mstyle$section("Model Results:"))
      cat("\n\n")
      tmp <- capture.output(.print.vector(res.table))
      .print.table(tmp, mstyle)

   }

   .space()

   invisible()

}
