print.rma.glmm <- function(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   if (!inherits(x, "rma.glmm"))
      stop("Argument 'x' must be an object of class \"rma.glmm\".")

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   if (x$method == "FE") {
      if (x$int.only) {
         cat("Fixed-Effects Model (k = ", x$k, ")", sep="")
      } else {
         cat("Fixed-Effects with Moderators Model (k = ", x$k, ")", sep="")
      }
   } else {
      if (x$int.only) {
         cat("Random-Effects Model (k = ", x$k, "; ", sep="")
      } else {
         cat("Mixed-Effects Model (k = ", x$k, "; ", sep="")
      }
      cat("tau^2 estimator: ", x$method, ")", sep="")
   }

   if (is.element(x$measure, c("OR","IRR"))) {
      cat("\n")
      if (x$model == "UM.FS")
         cat("Model Type: Unconditional Model with Fixed Study Effects")
      if (x$model == "UM.RS")
         cat("Model Type: Unconditional Model with Random Study Effects")
      if (x$model == "CM.AL")
         cat("Model Type: Conditional Model with Approximate Likelihood")
      if (x$model == "CM.EL")
         cat("Model Type: Conditional Model with Exact Likelihood")
   }

   if (showfit) {
      cat("\n")
      fs <- c(formatC(round(x$fit.stats$ML, digits=digits), digits=digits, format="f")) ### formatC(round()) - see comment below
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      print(fs, quote=FALSE, print.gap=2)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (x$method != "FE") {
      if (x$int.only) {
         cat("tau^2 (estimated amount of total heterogeneity): ", formatC(x$tau2, digits=ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits), format="f"), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , formatC(x$se.tau2, digits=digits, format="f"), ")") ), "\n", sep = "")
         cat("tau (square root of estimated tau^2 value):      ", ifelse(x$tau2>=0, formatC(sqrt(x$tau2), digits=ifelse(x$tau2 <= .Machine$double.eps*10,0,digits), format="f"), NA), "\n", sep = "")
         cat("I^2 (total heterogeneity / total variability):   ", ifelse(is.na(x$I2), NA, formatC(x$I2, digits=2, format="f")), "%", "\n", sep = "")
         cat("H^2 (total variability / sampling variability):  ", ifelse(is.na(x$H2), NA, formatC(x$H2, digits=2, format="f")), sep = "")
      } else {
         cat("tau^2 (estimated amount of residual heterogeneity):     ", formatC(x$tau2, digits=ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits), format="f"), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , formatC(x$se.tau2, digits=digits, format="f"), ")") ), "\n", sep = "")
         cat("tau (square root of estimated tau^2 value):             ", ifelse(x$tau2>=0, formatC(sqrt(x$tau2), digits=ifelse(x$tau2 <= .Machine$double.eps*10,0,digits), format="f"), NA), "\n", sep="")
         cat("I^2 (residual heterogeneity / unaccounted variability): ", ifelse(is.na(x$I2), NA, formatC(x$I2, digits=2, format="f")), "%", "\n", sep = "")
         cat("H^2 (unaccounted variability / sampling variability):   ", ifelse(is.na(x$H2), NA, formatC(x$H2, digits=2, format="f")), sep = "")
      }
      cat("\n\n")
   }

   if (!is.na(x$sigma2)) {
      cat("sigma^2 (estimated amount of study level variability): ", formatC(x$sigma2, digits=ifelse(abs(x$sigma2) <= .Machine$double.eps*10,0,digits), format="f"), "\n", sep = "")
      cat("sigma (square root of estimated sigma^2 value):        ", ifelse(x$sigma2>=0, formatC(sqrt(x$sigma2), digits=ifelse(x$sigma2 <= .Machine$double.eps*10,0,digits), format="f"), NA), "\n\n", sep = "")
   }

   if (!is.na(x$QE.Wld) || !is.na(x$QE.LRT)) {
      QE.Wld <- formatC(round(x$QE.Wld, digits=digits), digits=digits, format="f") ### formatC(round()) seems a bit redundant, but had problem with a "true" 0
      QE.LRT <- formatC(round(x$QE.LRT, digits=digits), digits=digits, format="f") ### (that even tested as TRUE with == 0) show up as -0.0000 (not sure why)

      if (nchar(QE.Wld) > nchar(QE.LRT))
         QE.LRT <- paste0(paste(rep(" ", nchar(QE.Wld) - nchar(QE.LRT)), collapse=""), QE.LRT)
      if (nchar(QE.LRT) > nchar(QE.Wld))
         QE.Wld <- paste0(paste(rep(" ", nchar(QE.LRT) - nchar(QE.Wld)), collapse=""), QE.Wld)

      if (x$int.only) {
         cat("Tests for Heterogeneity: \n")
         cat("Wld(df = ", x$QE.df, ") = ", QE.Wld, ", p-val ", .pval(x$QEp.Wld, digits=digits, showeq=TRUE, sep=" "), "\n", sep="")
         cat("LRT(df = ", x$QE.df, ") = ", QE.LRT, ", p-val ", .pval(x$QEp.LRT, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      } else {
         cat("Tests for Residual Heterogeneity: \n")
         cat("Wld(df = ", x$QE.df, ") = ", QE.Wld, ", p-val ", .pval(x$QEp.Wld, digits=digits, showeq=TRUE, sep=" "), "\n", sep="")
         cat("LRT(df = ", x$QE.df, ") = ", QE.LRT, ", p-val ", .pval(x$QEp.LRT, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      }
   }

   if (x$p > 1 && !is.na(x$QM)) {
      cat("Test of Moderators (coefficient(s) ", paste(x$btt, collapse=","),"): \n", sep="")
      if (x$knha) {
         cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      } else {
         cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      }
   }

   res.table <- cbind(estimate=c(x$b), se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   rownames(res.table) <- rownames(x$b)
   if (x$knha)
      colnames(res.table)[3] <- "tval"
   signif <- symnum(x$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(formatC(res.table, digits=digits, format="f"), signif)
      colnames(res.table)[7] <- ""
   } else {
      res.table <- formatC(res.table, digits=digits, format="f")
   }
   res.table[,4] <- .pval(x$pval, digits=digits)

   if (x$int.only)
      res.table <- res.table[1,]

   cat("Model Results:")
   cat("\n\n")
   if (x$int.only) {
      print(res.table, quote=FALSE, right=TRUE)
   } else {
      print(res.table, quote=FALSE, right=TRUE, print.gap=2)
   }
   cat("\n")
   if (signif.legend)
      cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")

   invisible()

}
