print.tes <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="tes")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   cat(mstyle$section(paste("Test of Excess Significance")))
   cat("\n\n")

   cat(mstyle$text("Observed Number of Significant Findings: "))
   cat(mstyle$result(x$O))
   cat(mstyle$result(paste0(" (out of ", x$k, ")")))
   cat("\n")
   cat(mstyle$text("Expected Number of Significant Findings: "))
   cat(mstyle$result(fmtx(x$E, digits[["est"]])))
   cat("\n")
   cat(mstyle$text("Observed Number / Expected Number:       "))
   cat(mstyle$result(fmtx(x$OEratio, digits[["est"]])))
   cat("\n\n")

   if (length(x$theta) == 1L) {
      cat(mstyle$text("Estimated Power of Tests (based on theta = "))
      cat(mstyle$result(fmtx(x$theta, digits[["est"]])))
      cat(mstyle$text(")"))
   } else {
      cat(mstyle$text("Estimated Power of Tests: "))
   }

   cat("\n\n")
   if (x$k > 5L) {
      power <- quantile(x$power)
      names(power) <- c("min", "q1", "median", "q3", "max")
   } else {
      power <- x$power
      names(power) <- seq_len(x$k)
   }
   tmp <- capture.output(.print.vector(fmtx(power, digits[["pval"]])))
   .print.table(tmp, mstyle)
   cat("\n")

   cat(mstyle$text("Test of Excess Significance: "))
   cat(mstyle$result(paste0("p ", fmtp(x$pval, digits[["pval"]], equal=TRUE, sep=TRUE))))
   if (x$test == "chi2") {
      cat(mstyle$result(paste0(" (X^2 = ", fmtx(x$X2, digits[["test"]]), ", df = 1)")))
   }
   if (x$test == "binom") {
      cat(mstyle$result(" (binomial test)"))
   }
   if (x$test == "exact") {
      cat(mstyle$result(" (exact test)"))
   }
   cat("\n")

   if (!is.null(x$theta.lim)) {
      cat(mstyle$text(paste0("Limit Estimate (theta_lim):  ")))
      if (is.na(x$theta.lim[1])) {
         cat(mstyle$result("not estimable"))
      } else {
         cat(mstyle$result(fmtx(x$theta.lim[1], digits[["est"]])))
      }
      if (length(x$theta.lim) == 2L) {
         cat(mstyle$result(", "))
         if (is.na(x$theta.lim[2])) {
            cat(mstyle$result("not estimable"))
         } else {
            cat(mstyle$result(fmtx(x$theta.lim[2], digits[["est"]])))
         }
      }

      if (any(!is.na(x$theta.lim)))
         cat(mstyle$result(paste0(" (where p = ", ifelse(x$tes.alternative == "two.sided", x$tes.alpha/2, x$tes.alpha), ")")))

      cat("\n")

   }

   .space()

   invisible()

}
