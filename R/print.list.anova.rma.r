print.list.anova.rma <- function(x, digits=x[[1]]$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="list.anova.rma")

   digits <- .get.digits(digits=digits, xdigits=x[[1]]$digits, dmiss=FALSE)

   .space()

   res.table <- data.frame(QM   = sapply(x, function(x) .fcf(x$QM, digits[["test"]])),
                           df   = sapply(x, function(x) round(x$QMdf[1], 2)),
                           pval = sapply(x, function(x) .pval(x$QMp, digits[["pval"]])))

   if (is.element(x[[1]]$test, c("knha","adhoc","t"))) {
      names(res.table)[1:2] <- c("F-val", "df1")
      res.table <- cbind(res.table[1:2], df2 = round(x[[1]]$QMdf[2], 2), res.table[3])
   }

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   .space()

   invisible()

}
