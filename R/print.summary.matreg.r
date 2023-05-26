print.summary.matreg <- function(x, digits=x$digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="summary.matreg")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   ### strip summary.matreg class from object (otherwise get recursion)

   class(x) <- class(x)[-1]

   ### print with showfit=TRUE

   print(x, digits=digits, signif.stars=signif.stars, signif.legend=signif.legend, ...)

   .space(FALSE)

   if (x$test == "t") {

      cat(mstyle$text("Residual standard error: "))
      cat(mstyle$result(fmtx(sqrt(x$mse), digits[["se"]])))
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
