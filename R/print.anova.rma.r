print.anova.rma <- function(x, digits, ...) {

   if (!inherits(x, "anova.rma"))
      stop("Argument 'x' must be an object of class \"anova.rma\".")

   if (missing(digits))
      digits <- x$digits

   if (x$type == "Wald.b") {

      cat("\n")

      cat("Test of Moderators (coefficient(s) ", .format.btt(x$btt),"): \n", sep="")
      #cat("Test of Moderators (coefficient(s) ", paste(x$btt, collapse=","),"): \n", sep="")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      } else {
         cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      }

   }

   if (x$type == "Wald.L") {

      cat("\n")

      if (x$m == 1) {
         cat("Hypothesis:")
      } else {
         cat("Hypotheses:")
      }

      print(x$hyp)

      cat("\nResults:\n")

      res.table <- cbind(estimate=c(x$Lb), se=x$se, zval=x$zval, pval=x$pval)
      if (is.element(x$test, c("knha","adhoc","t")))
         colnames(res.table)[3] <- "tval"
      rownames(res.table) <- paste0(seq_len(x$m), ":")
      res.table <- formatC(res.table, digits=digits, format="f")
      res.table[,4] <- .pval(x$pval, digits=digits)
      print(res.table, quote=FALSE, right=TRUE)

      cat("\n")

      if (!is.null(x$QM)) {
         if (x$m == 1) {
            cat("Test of Hypothesis:\n")
         } else {
            cat("Omnibus Test of Hypotheses:\n")
         }
         if (is.element(x$test, c("knha","adhoc","t"))) {
            cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
         } else {
            cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
         }
      }

   }

   if (x$type == "LRT") {

      res.table <- rbind(
         c(x$p.f, x$fit.stats.f["AIC"], x$fit.stats.f["BIC"], x$fit.stats.f["AICc"], x$fit.stats.f["ll"], NA,    NA,     x$QE.f, x$tau2.f, NA),
         c(x$p.r, x$fit.stats.r["AIC"], x$fit.stats.r["BIC"], x$fit.stats.r["AICc"], x$fit.stats.r["ll"], x$LRT, x$pval, x$QE.r, x$tau2.r, NA)
      )

      res.table[,seq.int(from=2, to=10)] <- formatC(res.table[,seq.int(from=2, to=10)], digits=digits, format="f")
      colnames(res.table) <- c("df", "AIC", "BIC", "AICc", "logLik", "LRT", "pval", "QE", "tau^2", "R^2")
      rownames(res.table) <- c("Full", "Reduced")

      res.table["Reduced","pval"] <- .pval(x$pval, digits=digits)

      res.table["Full",c("LRT","pval")] <- ""
      res.table["Full","R^2"] <- ""
      res.table["Reduced","R^2"] <- paste0(ifelse(is.na(x$R2), NA, formatC(x$R2, format="f", digits=2)), "%")

      ### remove tau^2 and R^2 columns if full model is FE or if dealing with rma.mv models

      if (x$method == "FE" || is.element("rma.mv", x$class.f))
         res.table <- res.table[,seq_len(8)]

      print(res.table, quote=FALSE, right=TRUE)

   }

   invisible()

}
