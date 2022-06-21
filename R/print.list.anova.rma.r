print.list.anova.rma <- function(x, digits=x[[1]]$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="list.anova.rma")

   digits <- .get.digits(digits=digits, xdigits=x[[1]]$digits, dmiss=FALSE)

   .space()

   res.table <- data.frame(spec  = names(x),
                           coefs = sapply(x, function(x) .format.btt(x$btt)),
                           QM    = sapply(x, function(x) .fcf(x$QM, digits[["test"]])),
                           df    = sapply(x, function(x) round(x$QMdf[1], 2)),
                           pval  = sapply(x, function(x) .pval(x$QMp, digits[["pval"]])))

   if (is.element(x[[1]]$test, c("knha","adhoc","t"))) {
      names(res.table)[3:4] <- c("Fval", "df1")
      res.table <- cbind(res.table[1:4], df2 = round(x[[1]]$QMdf[2], 2), res.table[5])
   }

   # if all btt specifications were numeric, remove the 'spec' column
   if (all(substr(res.table$spec, 1, 1) %in% as.character(1:9)))
       res.table$spec <- NULL

   rownames(res.table) <- NULL

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   .space()

   invisible()

}
