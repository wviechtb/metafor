print.permutest.rma.uni <- function(x, digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "permutest.rma.uni"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"permutest.rma.uni\"."))

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   if (!x$int.only) {
      cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val* ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val* ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   res.table <- cbind(estimate=c(x$beta), se=x$se, zval=x$zval, "pval*"=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   if (x$permci)
      colnames(res.table)[5:6] <- c("ci.lb*", "ci.ub*")
   rownames(res.table) <- rownames(x$beta)
   if (is.element(x$test, c("knha","adhoc","t")))
      colnames(res.table)[3] <- "tval"
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

   cat(mstyle$section("Model Results:"))
   cat("\n\n")
   if (x$int.only) {
      tmp <- capture.output(.print.vector(res.table))
   } else {
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   }
   .print.table(tmp, mstyle)

   cat("\n")
   if (signif.legend)
      cat(mstyle$legend("---\nSignif. codes: "), mstyle$legend(attr(signif, "legend")))
   cat("\n\n")

   invisible()

}
