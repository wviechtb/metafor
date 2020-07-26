print.matreg <- function(x, digits=x$digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "matreg"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"matreg\"."))

   if (!exists(".rmspace"))
      cat("\n")

   res.table <- x$tab
   res.table <- round(res.table, digits)
   signif <- symnum(x$tab$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(res.table, signif)
      colnames(res.table)[7] <- ""
   }

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
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
