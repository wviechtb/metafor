print.list.anova.rma <- function(x, digits=x[[1]]$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="list.anova.rma")

   digits <- .get.digits(digits=digits, xdigits=x[[1]]$digits, dmiss=FALSE)

   .space()

   res.table <- as.data.frame(x)

   if ("QM" %in% names(res.table))
      res.table$QM <- .fcf(res.table$QM, digits[["test"]])
   if ("Fval" %in% names(res.table))
      res.table$Fval <- .fcf(res.table$Fval, digits[["test"]])
   res.table$pval <- .pval(res.table$pval, digits[["pval"]])

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   .space()

   invisible()

}
