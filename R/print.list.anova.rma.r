print.list.anova.rma <- function(x, digits=x[[1]]$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

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

   res.table$pval <- fmtp(res.table$pval, digits[["pval"]])

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   .space()

   invisible()

}
