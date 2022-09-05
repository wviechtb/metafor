print.ranktest <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="ranktest")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   cat(mstyle$section("Rank Correlation Test for Funnel Plot Asymmetry"))
   cat("\n\n")
   cat(mstyle$result(paste0("Kendall's tau = ", fmtx(x$tau, digits[["est"]]), ", p ", fmtp(x$pval, digits[["pval"]], equal=TRUE, sep=TRUE))))
   cat("\n")
   #cat("H0: true tau is equal to 0\n\n")

   .space()

   invisible()

}
