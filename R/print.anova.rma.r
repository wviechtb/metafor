print.anova.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="anova.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   if (x$type == "Wald.btt") {

      if (is.element("rma.ls", x$class)) {
         cat(mstyle$section(paste0("Test of Location Coefficients (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      } else {
         cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      }
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$QMdf[1], ", df2 = ", round(x$QMdf[2], 2), ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits[["pval"]], showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QM(df = ", x$QMdf[1], ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits[["pval"]], showeq=TRUE, sep=" "))))
      }
      cat("\n")

   }

   if (x$type == "Wald.att") {

      cat(mstyle$section(paste0("Test of Scale Coefficients (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$att),"):")))
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$QSdf[1], ", df2 = ", round(x$QSdf[2], 2), ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits[["pval"]], showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QS(df = ", x$QSdf[1], ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits[["pval"]], showeq=TRUE, sep=" "))))
      }
      cat("\n")

   }

   if (x$type == "Wald.Xb") {

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

      if (is.element(x$test, c("knha","adhoc","t"))) {
         res.table <- data.frame(estimate=.fcf(c(x$Xb), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), tval=.fcf(x$zval, digits[["test"]]), df=round(x$ddf,2), pval=.pval(x$pval, digits[["pval"]]), stringsAsFactors=FALSE)
      } else {
         res.table <- data.frame(estimate=.fcf(c(x$Xb), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), pval=.pval(x$pval, digits[["pval"]]), stringsAsFactors=FALSE)
      }
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
            cat(mstyle$result(paste0("F(df1 = ", x$QMdf[1], ", df2 = ", round(x$QMdf[2], 2), ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits[["pval"]], showeq=TRUE, sep=" "))))
         } else {
            cat(mstyle$result(paste0("QM(df = ", x$QMdf[1], ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits[["pval"]], showeq=TRUE, sep=" "))))
         }
         cat("\n")
      }

   }

   if (x$type == "Wald.Za") {

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

      if (is.element(x$test, c("knha","adhoc","t"))) {
         res.table <- data.frame(estimate=.fcf(c(x$Za), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), tval=.fcf(x$zval, digits[["test"]]), df=round(x$ddf,2), pval=.pval(x$pval, digits[["pval"]]), stringsAsFactors=FALSE)
      } else {
         res.table <- data.frame(estimate=.fcf(c(x$Za), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), pval=.pval(x$pval, digits[["pval"]]), stringsAsFactors=FALSE)
      }
      rownames(res.table) <- paste0(seq_len(x$m), ":")
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

      if (!is.na(x$QS)) {
         cat("\n")
         if (x$m == 1) {
            cat(mstyle$section("Test of Hypothesis:"))
         } else {
            cat(mstyle$section("Omnibus Test of Hypotheses:"))
         }
         cat("\n")
         if (is.element(x$test, c("knha","adhoc","t"))) {
            cat(mstyle$result(paste0("F(df1 = ", x$QSdf[1], ", df2 = ", round(x$QSdf[2], 2), ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits[["pval"]], showeq=TRUE, sep=" "))))
         } else {
            cat(mstyle$result(paste0("QS(df = ", x$QSdf[1], ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits[["pval"]], showeq=TRUE, sep=" "))))
         }
         cat("\n")
      }

   }

   if (x$type == "LRT") {

      res.table <- data.frame(c(x$parms.f, x$parms.r),
                              c(.fcf(x$fit.stats.f["AIC"],  digits[["fit"]]), .fcf(x$fit.stats.r["AIC"],  digits[["fit"]])),
                              c(.fcf(x$fit.stats.f["BIC"],  digits[["fit"]]), .fcf(x$fit.stats.r["BIC"],  digits[["fit"]])),
                              c(.fcf(x$fit.stats.f["AICc"], digits[["fit"]]), .fcf(x$fit.stats.r["AICc"], digits[["fit"]])),
                              c(.fcf(x$fit.stats.f["ll"],   digits[["fit"]]), .fcf(x$fit.stats.r["ll"],   digits[["fit"]])),
                              c(NA, .fcf(x$LRT, digits[["test"]])),
                              c(NA, .pval(x$pval, digits[["pval"]])),
                              c(.fcf(x$QE.f, digits[["test"]]),  .fcf(x$QE.r, digits[["test"]])),
                              c(.fcf(x$tau2.f, digits[["var"]]), .fcf(x$tau2.r, digits[["var"]])),
                              c(NA, NA), stringsAsFactors=FALSE)

      colnames(res.table) <- c("df", "AIC", "BIC", "AICc", "logLik", "LRT", "pval", "QE", "tau^2", "R^2")
      rownames(res.table) <- c("Full", "Reduced")

      res.table["Full",c("LRT","pval")] <- ""
      res.table["Full","R^2"] <- ""
      res.table["Reduced","R^2"] <- paste0(.fcf(x$R2, digits[["het"]]), "%")

      ### remove tau^2 column if full model is a FE/EE/CE model

      if (is.element(x$method, c("FE","EE","CE")))
         res.table <- res.table[-which(names(res.table) == "tau^2")]

      ### remove R^2 column if full model is a rma.mv or rma.ls model

      if (is.element("rma.mv", x$class.f) || is.element("rma.ls", x$class.f))
         res.table <- res.table[-which(names(res.table) == "R^2")]

      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
      .print.table(tmp, mstyle)

   }

   .space()

   invisible()

}
