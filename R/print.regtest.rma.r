print.regtest.rma <- function(x, digits, ret.fit, ...) {

   if (class(x) != "regtest.rma")
      stop("Argument 'x' must be an object of class \"regtest.rma\".")

   if (missing(digits))
      digits <- x$digits

   if (missing(ret.fit))
      ret.fit <- x$ret.fit

   cat("\n")
   cat("Regression Test for Funnel Plot Asymmetry\n\n")
   if (x$model == "lm") {
      cat("model:     weighted regression with multiplicative dispersion\n")
   } else {
      cat("model:    ", ifelse(x$method=="FE", "fixed-effects", "mixed-effects"), "meta-regression model\n")
   }
   if (x$predictor == "sei")
      cat("predictor: standard error\n")
   if (x$predictor == "vi")
      cat("predictor: sampling variance\n")
   if (x$predictor == "ni")
      cat("predictor: sample size\n")
   if (x$predictor == "ninv")
      cat("predictor: inverse of the sample size\n")
   if (x$predictor == "sqrtni")
      cat("predictor: square-root sample size\n")
   if (x$predictor == "sqrtninv")
      cat("predictor: inverse of the square-root sample size\n")

   if (ret.fit) {
      print(x$fit)
   } else {
      cat("\n")
   }

   if (is.na(x$dfs)) {
      cat("test for funnel plot asymmetry: z = ", formatC(x$zval, digits=digits, format="f"), ", p ", .pval(x$pval, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
   } else {
      cat("test for funnel plot asymmetry: t = ", formatC(x$zval, digits=digits, format="f"), ", df = ", x$dfs, ", p ", .pval(x$pval, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
   }
   #cat("H0: coefficient for predictor is equal to 0\n\n")

   invisible()

}
