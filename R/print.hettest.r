print.hettest <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="hettest")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   for (j in 1:length(x$test)) {

      .space()

      cat(mstyle$section("Test for Heteroscedastic Heterogeneity"))
      cat("\n\n")
      cat(mstyle$text(paste0("Estimation method: "), ifelse(x$method == "ML", "Maximum likelihood", "Restricted maximum likelihood")))
      cat("\n")
      cat(mstyle$text(paste0("Test type:         "), sapply(x$test[j], switch, "lrt"="Likelihood ratio test", "wald"="Wald-type test", "score"="Score test", "ksn"="Kolmogorov-Smirnov test (normal)", "ksx2"="Kolmogorov-Smirnov test (chi-squared)", "adn"="Anderson-Darling test (normal)", "adx2"="Anderson-Darling test (chi-squared)", USE.NAMES=FALSE)))
      cat("\n")
      cat(mstyle$text(paste0("Bootstrapping:     "), ifelse(x$boot[j], "Yes", "No"), ifelse(x$boot[j], paste0(" (", sum(!is.na(x$statistic.boot[,j])), "/", x$iter, " iterations)"), "")))
      cat("\n\n")
      if (is.element(x$test[j], c("lrt", "wald", "score"))) {
         cat(mstyle$result(paste0("X^2(df = ", x$df[j], ") = ", fmtx(x$statistic[j], digits[["test"]]), ", p ", fmtp(x$pval[j], digits[["pval"]], equal=TRUE, sep=TRUE))))
      } else {
         cat(mstyle$result(paste0("statistic = ", fmtx(x$statistic[j], digits[["test"]]), ", p ", fmtp(x$pval[j], digits[["pval"]], equal=TRUE, sep=TRUE))))
      }
      cat("\n")

      .space()

   }

   invisible()

}
