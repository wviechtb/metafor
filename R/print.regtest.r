print.regtest <- function(x, digits=x$digits, ret.fit=x$ret.fit, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="regtest")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   if (!exists(".rmspace"))
      cat("\n")

   cat(mstyle$section("Regression Test for Funnel Plot Asymmetry"))
   cat("\n\n")
   if (x$model == "lm") {
      cat(mstyle$text("Model:     weighted regression with multiplicative dispersion"))
   } else {
      cat(mstyle$text(paste("Model:    ", ifelse(x$method=="FE", "fixed-effects", "mixed-effects"), "meta-regression model")))
   }
   cat("\n")
   if (x$predictor == "sei")
      cat(mstyle$text("Predictor: standard error"))
   if (x$predictor == "vi")
      cat(mstyle$text("Predictor: sampling variance"))
   if (x$predictor == "ni")
      cat(mstyle$text("Predictor: sample size"))
   if (x$predictor == "ninv")
      cat(mstyle$text("Predictor: inverse of the sample size"))
   if (x$predictor == "sqrtni")
      cat(mstyle$text("Predictor: square root sample size"))
   if (x$predictor == "sqrtninv")
      cat(mstyle$text("Predictor: inverse of the square root sample size"))

   cat("\n")

   if (ret.fit) {
      if (exists(".rmspace"))
         cat("\n")
      if (x$model == "lm") {
         print(summary(x$fit))
      } else {
         print(x$fit)
      }
      if (exists(".rmspace"))
         cat("\n")
   } else {
      cat("\n")
   }

   cat(mstyle$text("Test for Funnel Plot Asymmetry: "))
   if (is.na(x$dfs)) {
      cat(mstyle$result(paste0("z = ", .fcf(x$zval, digits[["test"]]), ", p ", .pval(x$pval, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
   } else {
      cat(mstyle$result(paste0("t = ", .fcf(x$zval, digits[["test"]]), ", df = ", x$dfs, ", p ", .pval(x$pval, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
   }
   cat("\n")

   if (!is.null(x$est)) {
      if (x$predictor == "sei")
         cat(mstyle$text("Limit Estimate (as sei -> 0):   "))
      if (x$predictor == "vi")
         cat(mstyle$text("Limit Estimate (as vi -> 0):    "))
      if (x$predictor %in% c("ninv", "sqrtninv"))
         cat(mstyle$text("Limit Estimate (as ni -> inf):  "))
      cat(mstyle$result(paste0("b = ", .fcf(x$est, digits[["est"]]), " (CI: ", .fcf(x$ci.lb, digits[["est"]]), ", ", .fcf(x$ci.ub, digits[["est"]]), ")")))
      cat("\n")
   }

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
