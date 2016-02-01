print.robust.rma <- function(x, digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   if (!inherits(x, "robust.rma"))
      stop("Argument 'x' must be an object of class \"robust.rma\".")

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   cat("Number of outcomes:  ", x$k, "\n")
   cat("Number of clusters:  ", x$n, "\n")

   if (all(x$tcl[1] == x$tcl)) {
      cat("Outcomes per cluster:", x$tcl[1], "\n")
   } else {
      cat("Outcomes per cluster: ", min(x$tcl), "-", max(x$tcl), " (mean: ", formatC(mean(x$tcl), format="f", digits=2), ", median: ", median(x$tcl), ")\n", sep="")
   }
   cat("\n")

   if (x$p > 1) {
      cat("Test of Moderators (coefficient(s) ", paste(x$btt, collapse=","),"): \n", sep="")
      cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
   }

   res.table <- cbind(estimate=c(x$b), se=x$se, tval=x$tval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   rownames(res.table) <- rownames(x$b)
   signif <- symnum(x$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(formatC(res.table, digits=digits, format="f"), signif)
      colnames(res.table)[7] <- ""
   } else {
      res.table <- formatC(res.table, digits=digits, format="f")
   }
   res.table[,4] <- .pval(x$pval, digits=digits)

   if (x$int.only)
      res.table <- res.table[1,]

   cat("Model Results:\n\n")
   if (x$int.only) {
      print(res.table, quote=FALSE, right=TRUE)
   } else {
      print(res.table, quote=FALSE, right=TRUE, print.gap=2)
   }
   cat("\n")
   if (signif.legend)
      cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")

   invisible()

}
