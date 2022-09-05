print.regtest <- function(x, digits=x$digits, ret.fit=x$ret.fit, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="regtest")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   cat(mstyle$section("Regression Test for Funnel Plot Asymmetry"))
   cat("\n\n")
   if (x$model == "lm") {
      cat(mstyle$text("Model:     weighted regression with multiplicative dispersion"))
   } else {
      cat(mstyle$text(paste("Model:    ", ifelse(is.element(x$method, c("FE","EE","CE")), "fixed-effects", "mixed-effects"), "meta-regression model")))
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
      .space(FALSE)
      if (x$model == "lm") {
         print(summary(x$fit))
      } else {
         print(x$fit)
      }
      .space(FALSE)
   } else {
      cat("\n")
   }

   cat(mstyle$text("Test for Funnel Plot Asymmetry: "))
   if (is.na(x$ddf)) {
      cat(mstyle$result(fmtt(x$zval, "z", pval=x$pval, pname="p", format=2, digits=digits)))
   } else {
      cat(mstyle$result(fmtt(x$zval, "t", df=x$ddf, pval=x$pval, pname="p", format=2, digits=digits)))
   }
   cat("\n")

   if (!is.null(x$est)) {
      if (x$predictor == "sei")
         cat(mstyle$text("Limit Estimate (as sei -> 0):   "))
      if (x$predictor == "vi")
         cat(mstyle$text("Limit Estimate (as vi -> 0):    "))
      if (x$predictor %in% c("ninv", "sqrtninv"))
         cat(mstyle$text("Limit Estimate (as ni -> inf):  "))
      cat(mstyle$result(paste0("b = ", fmtx(x$est, digits[["est"]]), " (CI: ", fmtx(x$ci.lb, digits[["est"]]), ", ", fmtx(x$ci.ub, digits[["est"]]), ")")))
      cat("\n")
   }

   .space()

   invisible()

}
