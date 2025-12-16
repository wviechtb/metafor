print.rma.glmm <- function(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(x), must="rma.glmm")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   ddd <- list(...)

   .chkdots(ddd, c("num"))

   .space()

   if (is.element(x$method, c("FE","EE","CE"))) {
      if (x$int.only) {
         cat(mstyle$section(sapply(x$method, switch, "FE"="Fixed-Effects Model", "EE"="Equal-Effects Model", "CE"="Common-Effects Model", USE.NAMES=FALSE)))
      } else {
         cat(mstyle$section("Fixed-Effects with Moderators Model"))
      }
      cat(mstyle$section(paste0(" (k = ", x$k, ")")))
   } else {
      if (x$int.only) {
         cat(mstyle$section("Random-Effects Model"))
      } else {
         cat(mstyle$section("Mixed-Effects Model"))
      }
      cat(mstyle$section(paste0(" (k = ", x$k, "; ")))
      cat(mstyle$section(paste0("tau^2 estimator: ", x$method, ")")))
   }

   if (is.element(x$measure, c("OR","IRR"))) {
      cat("\n")
      if (x$model == "UM.FS")
         cat(mstyle$section("Model Type: Unconditional Model with Fixed Study Effects"))
      if (x$model == "UM.RS")
         cat(mstyle$section("Model Type: Unconditional Model with Random Study Effects"))
      if (x$model == "CM.AL")
         cat(mstyle$section("Model Type: Conditional Model with Approximate Likelihood"))
      if (x$model == "CM.EL")
         cat(mstyle$section("Model Type: Conditional Model with Exact Likelihood"))
   }

   if (showfit) {
      cat("\n")
      fs <- fmtx(x$fit.stats$ML, digits[["fit"]])
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      #tmp[1] <- paste0(tmp[1], "\u200b")
      .print.table(tmp, mstyle)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (!is.element(x$method, c("FE","EE","CE"))) {
      if (x$int.only) {
         cat(mstyle$text("tau^2 (estimated amount of total heterogeneity): "))
         cat(mstyle$result(paste0(fmtx(x$tau2, digits[["var"]], thresh=.Machine$double.eps*10), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , fmtx(x$se.tau2, digits[["sevar"]]), ")")))))
         cat("\n")
         cat(mstyle$text("tau (square root of estimated tau^2 value):      "))
         cat(mstyle$result(fmtx(.sqrt(x$tau2), digits[["var"]], thresh=.Machine$double.eps*10)))
         cat("\n")
         cat(mstyle$text("I^2 (total heterogeneity / total variability):   "))
         cat(mstyle$result(paste0(fmtx(x$I2, 2), "%")))
         cat("\n")
         cat(mstyle$text("H^2 (total variability / sampling variability):  "))
         cat(mstyle$result(fmtx(x$H2, 2)))
      } else {
         cat(mstyle$text("tau^2 (estimated amount of residual heterogeneity):     "))
         cat(mstyle$result(paste0(fmtx(x$tau2, digits[["var"]], thresh=.Machine$double.eps*10), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , fmtx(x$se.tau2, digits[["sevar"]]), ")")))))
         cat("\n")
         cat(mstyle$text("tau (square root of estimated tau^2 value):             "))
         cat(mstyle$result(fmtx(.sqrt(x$tau2), digits[["var"]], thresh=.Machine$double.eps*10)))
         cat("\n")
         cat(mstyle$text("I^2 (residual heterogeneity / unaccounted variability): "))
         cat(mstyle$result(paste0(fmtx(x$I2, 2), "%")))
         cat("\n")
         cat(mstyle$text("H^2 (unaccounted variability / sampling variability):   "))
         cat(mstyle$result(fmtx(x$H2, 2)))
      }
      cat("\n\n")
   }

   if (!is.na(x$sigma2)) {
      cat(mstyle$text("sigma^2 (estimated amount of study level variability): "))
      cat(mstyle$result(fmtx(x$sigma2, digits[["var"]], thresh=.Machine$double.eps*10)))
      cat("\n")
      cat(mstyle$text("sigma (square root of estimated sigma^2 value):        "))
      cat(mstyle$result(fmtx(.sqrt(x$sigma2), digits[["var"]], thresh=.Machine$double.eps*10)))
      cat("\n\n")
   }

   if (!is.na(x$QE.Wld) || !is.na(x$QE.LRT)) {
      QE.Wld <- fmtx(x$QE.Wld, digits[["test"]])
      QE.LRT <- fmtx(x$QE.LRT, digits[["test"]])
      nchar.Wld <- nchar(QE.Wld, keepNA=FALSE)
      nchar.LRT <- nchar(QE.LRT, keepNA=FALSE)

      if (nchar.Wld > nchar.LRT)
         QE.LRT <- paste0(paste(rep(" ", nchar.Wld - nchar.LRT), collapse=""), QE.LRT)
      if (nchar.LRT > nchar.Wld)
         QE.Wld <- paste0(paste(rep(" ", nchar.LRT - nchar.Wld), collapse=""), QE.Wld)

      if (x$int.only) {
         cat(mstyle$section("Tests for Heterogeneity:"))
      } else {
         cat(mstyle$section("Tests for Residual Heterogeneity:"))
      }
      cat("\n")
      cat(mstyle$result(fmtt(x$QE.Wld, "Wld", df=x$QE.df, pval=x$QEp.Wld, digits=digits)))
      cat("\n")
      cat(mstyle$result(fmtt(x$QE.LRT, "LRT", df=x$QE.df, pval=x$QEp.LRT, digits=digits)))
      cat("\n\n")
   }

   if (x$p > 1L && !is.na(x$QM)) {
      cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(fmtt(x$QM, "F", df1=x$QMdf[1], df2=x$QMdf[2], pval=x$QMp, digits=digits)))
      } else {
         cat(mstyle$result(fmtt(x$QM, "QM", df=x$QMdf[1], pval=x$QMp, digits=digits)))
      }
      cat("\n\n")
   }

   if (is.element(x$test, c("knha","adhoc","t"))) {
      res.table <- data.frame(estimate=fmtx(c(x$beta), digits[["est"]]), se=fmtx(x$se, digits[["se"]]), tval=fmtx(x$zval, digits[["test"]]), df=round(x$ddf,2), pval=fmtp(x$pval, digits[["pval"]]), ci.lb=fmtx(x$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
   } else {
      res.table <- data.frame(estimate=fmtx(c(x$beta), digits[["est"]]), se=fmtx(x$se, digits[["se"]]), zval=fmtx(x$zval, digits[["test"]]), pval=fmtp(x$pval, digits[["pval"]]), ci.lb=fmtx(x$ci.lb, digits[["ci"]]), ci.ub=fmtx(x$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
   }
   rownames(res.table) <- rownames(x$beta)
   signif <- symnum(x$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(res.table, signif)
      colnames(res.table)[ncol(res.table)] <- ""
   }

   if (isTRUE(ddd$num)) {
      width <- nchar(nrow(res.table))
      rownames(res.table) <- paste0(formatC(seq_len(nrow(res.table)), format="d", width=width), ") ", rownames(res.table))
   }

   if (x$int.only)
      res.table <- res.table[1,]

   cat(mstyle$section("Model Results:"))
   cat("\n\n")
   if (x$int.only) {
      tmp <- capture.output(.print.vector(res.table))
   } else {
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   }
   #tmp[1] <- paste0(tmp[1], "\u200b")
   .print.table(tmp, mstyle)

   if (signif.legend) {
      cat("\n")
      cat(mstyle$legend("---"))
      cat("\n")
      cat(mstyle$legend("Signif. codes: "), mstyle$legend(attr(signif, "legend")))
      cat("\n")
   }

   .space()

   invisible()

}
