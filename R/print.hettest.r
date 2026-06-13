print.hettest <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="hettest")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   cat(mstyle$section("Test for Heteroscedastic Heterogeneity"))
   cat("\n\n")
   cat(mstyle$text(paste0("Estimation method: "), ifelse(x$method == "ML", "Maximum likelihood", "Restricted maximum likelihood")))
   cat("\n")
   cat(mstyle$text(paste0("Test type:         "), sapply(x$test, switch, "lrt"="Likelihood ratio test", "wald"="Wald-type test", "score"="Score test", USE.NAMES=FALSE)))
   cat("\n")
   cat(mstyle$text(paste0("Bootstrapping:     "), ifelse(x$boot, "Yes", "No"), ifelse(x$boot, paste0(" (", sum(!is.na(x$x2.boot)), "/", x$iter, " iterations)"), "")))
   cat("\n\n")
   cat(mstyle$result(paste0("X^2(df = ", x$df, ") = ", fmtx(x$x2, digits[["test"]]), ", p ", fmtp(x$pval, digits[["pval"]], equal=TRUE, sep=TRUE))))
   cat("\n")

   .space()

   invisible()

}
