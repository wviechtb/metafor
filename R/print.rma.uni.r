print.rma.uni <- function(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   if (!inherits(x, "rma.uni"))
      stop("Argument 'x' must be an object of class \"rma.uni\".")

   if (missing(digits))
      digits <- x$digits

   if (inherits(x, "rma.uni.trimfill")) {
      cat("\nEstimated number of missing studies on the ", x$side, " side: ", x$k0, " (SE = ", ifelse(is.na(x$se.k0), NA, formatC(x$se.k0, digits=digits, format="f")), ")\n", sep="")
      if (x$k0.est == "R0")
         cat("Test of H0: no missing studies on the ", x$side, " side: p-val ", .pval(x$p.k0, digits=digits, showeq=TRUE, sep=" "), "\n", sep="")
   }

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
      if (x$tau2.fix) {
         cat("user-specified tau^2 value)", sep="")
      } else {
         cat("tau^2 estimator: ", x$method, ")", sep="")
      }
   }

   if (showfit) {
      cat("\n")
      if (x$method == "REML") {
         fs <- c(formatC(round(x$fit.stats$REML, digits=digits), digits=digits, format="f"))
      } else {
         fs <- c(formatC(round(x$fit.stats$ML, digits=digits), digits=digits, format="f"))
      }
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      print(fs, quote=FALSE, print.gap=2)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (x$method != "FE") {
      if (x$int.only) {
         if (x$model == "rma.uni") {
            cat("tau^2 (estimated amount of total heterogeneity): ", formatC(x$tau2, digits=ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits), format="f"), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , formatC(x$se.tau2, digits=digits, format="f"), ")") ), "\n", sep = "")
            cat("tau (square root of estimated tau^2 value):      ", ifelse(x$tau2>=0, formatC(sqrt(x$tau2), digits=ifelse(x$tau2 <= .Machine$double.eps*10,0,digits), format="f"), NA), "\n", sep = "")
         }
         if (!is.na(x$I2))
            cat("I^2 (total heterogeneity / total variability):   ", ifelse(is.na(x$I2), NA, formatC(x$I2, digits=2, format="f")), "%", "\n", sep = "")
         if (!is.na(x$H2))
            cat("H^2 (total variability / sampling variability):  ", ifelse(is.na(x$H2), NA, formatC(x$H2, digits=2, format="f")), "\n", sep = "")
      } else {
         if (x$model == "rma.uni") {
            cat("tau^2 (estimated amount of residual heterogeneity):     ", formatC(x$tau2, digits=ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits), format="f"), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , formatC(x$se.tau2, digits=digits, format="f"), ")") ), "\n", sep = "")
            cat("tau (square root of estimated tau^2 value):             ", ifelse(x$tau2>=0, formatC(sqrt(x$tau2), digits=ifelse(x$tau2 <= .Machine$double.eps*10,0,digits), format="f"), NA), "\n", sep="")
         }
         if (!is.na(x$I2))
            cat("I^2 (residual heterogeneity / unaccounted variability): ", ifelse(is.na(x$I2), NA, formatC(x$I2, digits=2, format="f")), "%", "\n", sep = "")
         if (!is.na(x$H2))
            cat("H^2 (unaccounted variability / sampling variability):   ", ifelse(is.na(x$H2), NA, formatC(x$H2, digits=2, format="f")), "\n", sep = "")
         if (!is.null(x$R2))
            cat("R^2 (amount of heterogeneity accounted for):            ", ifelse(is.na(x$R2), NA, formatC(x$R2, digits=2, format="f")), "%", "\n", sep = "")
      }
      cat("\n")
   }

   if (!is.na(x$QE)) {
      if (x$int.only) {
         cat("Test for Heterogeneity: \n")
         cat("Q(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      } else {
         cat("Test for Residual Heterogeneity: \n")
         cat("QE(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      }
   }

   if (x$p > 1) {
      cat("Test of Moderators (coefficient(s) ", .format.btt(x$btt),"): \n", sep="")
      #cat("Test of Moderators (coefficient(s) ", paste(x$btt, collapse=","),"): \n", sep="")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      } else {
         cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "), "\n\n", sep="")
      }
   }

   res.table <- cbind(estimate=c(x$b), se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   rownames(res.table) <- rownames(x$b)
   if (is.element(x$test, c("knha","adhoc","t")))
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

   if (x$model == "rma.uni") {
      cat("Model Results:\n\n")
   } else {
      cat("Model Results (Location):\n\n")
   }
   if (x$int.only) {
      print(res.table, quote=FALSE, right=TRUE)
   } else {
      print(res.table, quote=FALSE, right=TRUE, print.gap=2)
   }

   if (x$model == "rma.ls") {

      res.table <- cbind(estimate=c(x$b.tau2), se=x$se.tau2, zval=x$zval.tau2, pval=x$pval.tau2, ci.lb=x$ci.lb.tau2, ci.ub=x$ci.ub.tau2)
      rownames(res.table) <- rownames(x$b.tau2)
      signif <- symnum(x$pval.tau2, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
      if (signif.stars) {
         res.table <- cbind(formatC(res.table, digits=digits, format="f"), signif)
         colnames(res.table)[7] <- ""
      } else {
         res.table <- formatC(res.table, digits=digits, format="f")
      }
      res.table[,4] <- .pval(x$pval.tau2, digits=digits)

      if (length(x$b.tau2) == 1)
         res.table <- res.table[1,]

      cat("\nModel Results (Scale):\n\n")

      if (length(x$b.tau2) == 1) {
         print(res.table, quote=FALSE, right=TRUE)
      } else {
         print(res.table, quote=FALSE, right=TRUE, print.gap=2)
      }

   }

   cat("\n")
   if (signif.legend)
      cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")

   invisible()

}
