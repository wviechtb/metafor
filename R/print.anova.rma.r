print.anova.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "anova.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"anova.rma\"."))

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   if (!exists(".rmspace"))
      cat("\n")

   if (x$type == "Wald.b") {

      cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QM(df = ", x$m, ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
      }
      cat("\n")

   }

   if (x$type == "Wald.L") {

      if (x$m == 1) {
         cat(mstyle$section("Hypothesis:"))
      } else {
         cat(mstyle$section("Hypotheses:"))
      }

      tmp <- capture.output(print(x$hyp))
      .print.output(tmp, mstyle$text)

      cat("\n")
      cat(mstyle$section("Results:"))
      cat("\n")

      res.table <- data.frame(estimate=.fcf(c(x$Lb), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), pval=.pval(x$pval, digits=digits[["pval"]]))
      if (is.element(x$test, c("knha","adhoc","t")))
         colnames(res.table)[3] <- "tval"
      rownames(res.table) <- paste0(seq_len(x$m), ":")
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

      if (!is.na(x$QM)) {
         cat("\n")
         if (x$m == 1) {
            cat(mstyle$section("Test of Hypothesis:"))
         } else {
            cat(mstyle$section("Omnibus Test of Hypotheses:"))
         }
         cat("\n")
         if (is.element(x$test, c("knha","adhoc","t"))) {
            cat(mstyle$result(paste0("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
         } else {
            cat(mstyle$result(paste0("QM(df = ", x$m, ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
         }
         cat("\n")
      }

   }

   if (x$type == "LRT") {

      res.table <- rbind(
         c(x$p.f, .fcf(x$fit.stats.f[c("AIC","BIC","AICc","ll")], digits[["fit"]]), NA, NA, .fcf(x$QE.f, digits[["test"]]), .fcf(x$tau2.f, digits[["var"]]), NA),
         c(x$p.r, .fcf(x$fit.stats.r[c("AIC","BIC","AICc","ll")], digits[["fit"]]), .fcf(x$LRT, digits[["test"]]), .pval(x$pval, digits=digits[["pval"]]), .fcf(x$QE.r, digits[["test"]]), .fcf(x$tau2.r, digits[["var"]]), NA)
      )

      colnames(res.table) <- c("df", "AIC", "BIC", "AICc", "logLik", "LRT", "pval", "QE", "tau^2", "R^2")
      rownames(res.table) <- c("Full", "Reduced")

      res.table["Full",c("LRT","pval")] <- ""
      res.table["Full","R^2"] <- ""
      res.table["Reduced","R^2"] <- paste0(.fcf(x$R2, digits[["het"]]), "%")

      ### remove tau^2 and R^2 columns if full model is FE or if dealing with rma.mv models

      if (x$method == "FE" || is.element("rma.mv", x$class.f))
         res.table <- res.table[,seq_len(8)]

      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

   }

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
