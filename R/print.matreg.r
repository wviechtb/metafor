print.matreg <- function(x, digits=x$digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="matreg")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   .space()

   if (x$test == "t") {
      res.table <- data.frame(estimate=fmtx(c(x$tab$beta), digits[["est"]]), se=fmtx(x$tab$se, digits[["se"]]), tval=fmtx(x$tab$tval, digits[["test"]]), df=round(x$tab$df,2), pval=fmtp(x$tab$pval, digits[["pval"]]), ci.lb=fmtx(x$tab$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$tab$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
      res.table$df[is.na(x$tab$beta)] <- NA_real_
   } else {
      res.table <- data.frame(estimate=fmtx(c(x$tab$beta), digits[["est"]]), se=fmtx(x$tab$se, digits[["se"]]), zval=fmtx(x$tab$zval, digits[["test"]]), pval=fmtp(x$tab$pval, digits[["pval"]]), ci.lb=fmtx(x$tab$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$tab$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
   }

   rownames(res.table) <- rownames(x$tab)

   signif <- symnum(x$tab$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

   if (signif.stars) {
      res.table <- cbind(res.table, signif)
      colnames(res.table)[ncol(res.table)] <- ""
   }

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   #tmp[1] <- paste0(tmp[1], "\u200b")
   .print.table(tmp, mstyle)

   if (signif.legend) {
      cat("\n")
      cat(mstyle$legend("---"))
      cat("\n")
      cat(mstyle$legend("Signif. codes: "), mstyle$legend(attr(signif, "legend")))
      cat("\n")
   }

   .space()

   invisible()

}

summary.matreg <- function(object, digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="matreg")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=object$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=object$digits, dmiss=FALSE)
   }

   object$digits <- digits

   class(object) <- c("summary.matreg", class(object))
   return(object)

}

print.summary.matreg <- function(x, digits=x$digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="summary.matreg")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   # strip summary.matreg class from object (otherwise get recursion)
   class(x) <- class(x)[-1]

   print(x, digits=digits, signif.stars=signif.stars, signif.legend=signif.legend, ...)

   .space(FALSE)

   if (x$test == "t") {

      cat(mstyle$text("Residual standard error: "))
      cat(mstyle$result(fmtx(sqrt(x$sigma2.reml), digits[["se"]])))
      cat(mstyle$text(paste0(" on ", x$Fdf[2], " degrees of freedom\n")))

      cat(mstyle$text("Multiple R-squared: "))
      cat(mstyle$result(fmtx(x$R2, digits[["het"]])))
      cat(mstyle$text(",  Adjusted R-squared: "))
      cat(mstyle$result(fmtx(x$R2adj, digits[["het"]])))

      cat("\n")

      cat(mstyle$text("F-statistic: "))
      cat(mstyle$result(fmtx(x$F[["value"]], digits[["test"]])))
      cat(mstyle$text(paste0(" on ", x$Fdf[1], " and ", x$Fdf[2], " DF,  p-value: ")))
      cat(mstyle$result(fmtp(x$Fp, digits[["pval"]], equal=FALSE, sep=FALSE)))

   } else {

      cat(mstyle$result("R^2: "))
      cat(mstyle$result(fmtx(x$R2, digits[["het"]])))
      cat(mstyle$result(", "))
      cat(mstyle$result(fmtt(x$QM, "QM", df=x$QMdf[1], pval=x$QMp, digits=digits)))

   }

   cat("\n")

   .space()

   invisible()

}
