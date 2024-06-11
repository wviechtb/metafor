print.list.anova.rma <- function(x, digits=x[[1]]$digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="list.anova.rma")

   digits <- .get.digits(digits=digits, xdigits=x[[1]]$digits, dmiss=FALSE)

   .space()

   res.table <- as.data.frame(x)

   if ("QM" %in% names(res.table))
      res.table$QM <- fmtx(res.table$QM, digits[["test"]])
   if ("QS" %in% names(res.table))
      res.table$QS <- fmtx(res.table$QS, digits[["test"]])
   if ("Fval" %in% names(res.table))
      res.table$Fval <- fmtx(res.table$Fval, digits[["test"]])

   signif <- symnum(res.table$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

   res.table$pval <- fmtp(res.table$pval, digits[["pval"]])

   if (getOption("show.signif.stars")) {
      res.table <- cbind(res.table, signif)
      colnames(res.table)[ncol(res.table)] <- ""
   }

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   .space()

   invisible()

}
