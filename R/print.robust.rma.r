print.robust.rma <- function(x, digits=x$digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="robust.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   if (!exists(".rmspace"))
      cat("\n")

   cat(mstyle$text("Number of outcomes:   "))
   cat(mstyle$result(x$k))
   cat("\n")
   cat(mstyle$text("Number of clusters:   "))
   cat(mstyle$result(x$n))
   cat("\n")

   cat(mstyle$text("Outcomes per cluster: "))
   if (all(x$tcl[1] == x$tcl)) {
      cat(mstyle$result(x$tcl[1]))
   } else {
      cat(mstyle$result(paste0(min(x$tcl), "-", max(x$tcl), " (mean: ", .fcf(mean(x$tcl), digits=2), ", median: ", round(median(x$tcl), digits=2), ")")))
   }
   cat("\n\n")

   if (x$p > 1L && !is.na(x$QM)) {
      cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      cat("\n")
      cat(mstyle$result(paste0("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
      cat("\n\n")
   }

   res.table <- data.frame(estimate=.fcf(c(x$beta), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), pval=.pval(x$pval, digits[["pval"]]), ci.lb=.fcf(x$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$ci.ub, digits[["ci"]]))
   rownames(res.table) <- rownames(x$beta)
   if (is.element(x$test, c("knha","adhoc","t")))
      colnames(res.table)[3] <- "tval"
   signif <- symnum(x$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(res.table, signif)
      colnames(res.table)[7] <- ""
   }

   if (x$int.only)
      res.table <- res.table[1,]

   cat(mstyle$section("Model Results:"))
   cat("\n\n")
   if (x$int.only) {
      tmp <- capture.output(.print.vector(res.table))
   } else {
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   }
   .print.table(tmp, mstyle)

   if (signif.legend) {
      cat("\n")
      cat(mstyle$legend("---\nSignif. codes: "), mstyle$legend(attr(signif, "legend")))
      cat("\n")
   }

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
