print.rma.uni <- function(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.uni"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.uni\"."))

   if (missing(digits))
      digits <- x$digits

   if (inherits(x, "rma.uni.trimfill")) {
      cat("\n")
      cat(mstyle$text(paste0("Estimated number of missing studies on the ", x$side, " side: ")))
      cat(mstyle$result(paste0(x$k0, " (SE = ", ifelse(is.na(x$se.k0), NA, formatC(x$se.k0, digits=digits, format="f")), ")")))
      cat("\n")
      if (x$k0.est == "R0") {
         cat(mstyle$text(paste0("Test of H0: no missing studies on the ", x$side, " side:     ")))
         cat(paste0(rep(" ", nchar(x$k0)), collapse=""))
         cat(mstyle$result(paste0("p-val ", .pval(x$p.k0, digits=digits, showeq=TRUE, sep=" "))))
         cat("\n")
      }
   }

   cat("\n")

   if (x$method == "FE") {
      if (x$int.only) {
         cat(mstyle$section("Fixed-Effects Model"))
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
      if (x$tau2.fix) {
         cat(mstyle$section("user-specified tau^2 value)"))
      } else {
         cat(mstyle$section(paste0("tau^2 estimator: ", x$method, ")")))
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
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      .print.table(tmp, mstyle)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (x$method != "FE") {
      if (x$int.only) {
         if (x$model == "rma.uni") {
            cat(mstyle$text("tau^2 (estimated amount of total heterogeneity): "))
            cat(mstyle$result(paste0(formatC(x$tau2, digits=ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits), format="f"), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , formatC(x$se.tau2, digits=digits, format="f"), ")")))))
            cat("\n")
            cat(mstyle$text("tau (square root of estimated tau^2 value):      "))
            cat(mstyle$result(paste0(ifelse(x$tau2>=0, formatC(sqrt(x$tau2), digits=ifelse(x$tau2 <= .Machine$double.eps*10,0,digits), format="f"), NA))))
            cat("\n")
         }
         if (!is.na(x$I2)) {
            cat(mstyle$text("I^2 (total heterogeneity / total variability):   "))
            cat(mstyle$result(paste0(ifelse(is.na(x$I2), NA, formatC(x$I2, digits=2, format="f")), "%")))
            cat("\n")
         }
         if (!is.na(x$H2)) {
            cat(mstyle$text("H^2 (total variability / sampling variability):  "))
            cat(mstyle$result(paste0(ifelse(is.na(x$H2), NA, formatC(x$H2, digits=2, format="f")))))
            cat("\n")
         }
      } else {
         if (x$model == "rma.uni") {
            cat(mstyle$text("tau^2 (estimated amount of residual heterogeneity):     "))
            cat(mstyle$result(paste0(formatC(x$tau2, digits=ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits), format="f"), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , formatC(x$se.tau2, digits=digits, format="f"), ")")))))
            cat("\n")
            cat(mstyle$text("tau (square root of estimated tau^2 value):             "))
            cat(mstyle$result(paste0(ifelse(x$tau2>=0, formatC(sqrt(x$tau2), digits=ifelse(x$tau2 <= .Machine$double.eps*10,0,digits), format="f"), NA))))
            cat("\n")
         }
         if (!is.na(x$I2)) {
            cat(mstyle$text("I^2 (residual heterogeneity / unaccounted variability): "))
            cat(mstyle$result(paste0(ifelse(is.na(x$I2), NA, formatC(x$I2, digits=2, format="f")), "%")))
            cat("\n")
         }
         if (!is.na(x$H2)) {
            cat(mstyle$text("H^2 (unaccounted variability / sampling variability):   "))
            cat(mstyle$result(paste0(ifelse(is.na(x$H2), NA, formatC(x$H2, digits=2, format="f")))))
            cat("\n")
         }
         if (!is.null(x$R2)) {
            cat(mstyle$text("R^2 (amount of heterogeneity accounted for):            "))
            cat(mstyle$result(paste0(ifelse(is.na(x$R2), NA, formatC(x$R2, digits=2, format="f")), "%")))
            cat("\n")
         }
      }
      cat("\n")
   }

   if (!is.na(x$QE)) {
      if (x$int.only) {
         cat(mstyle$section("Test for Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("Q(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$section("Test for Residual Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("QE(df = ", x$k-x$p, ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   if (x$p > 1 && !is.na(x$QM)) {
      cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QM(df = ", x$m, ") = ", formatC(x$QM, digits=digits, format="f"), ", p-val ", .pval(x$QMp, digits=digits, showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   res.table <- cbind(estimate=c(x$beta), se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   rownames(res.table) <- rownames(x$beta)
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
      cat(mstyle$section("Model Results:"))
   } else {
      cat(mstyle$section("Model Results (Location):"))
   }
   cat("\n\n")
   if (x$int.only) {
      tmp <- capture.output(.print.vector(res.table))
   } else {
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   }
   .print.table(tmp, mstyle)

   if (x$model == "rma.ls") {

      res.table <- cbind(estimate=c(x$alpha), se=x$se.alpha, zval=x$zval.alpha, pval=x$pval.alpha, ci.lb=x$ci.lb.alpha, ci.ub=x$ci.ub.alpha)
      rownames(res.table) <- rownames(x$alpha)
      if (is.element(x$test, c("t")))
         colnames(res.table)[3] <- "tval"
      signif <- symnum(x$pval.alpha, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
      if (signif.stars) {
         res.table <- cbind(formatC(res.table, digits=digits, format="f"), signif)
         colnames(res.table)[7] <- ""
      } else {
         res.table <- formatC(res.table, digits=digits, format="f")
      }
      res.table[,4] <- .pval(x$pval.alpha, digits=digits)

      if (length(x$alpha) == 1)
         res.table <- res.table[1,]

      cat("\n")
      cat(mstyle$section("Model Results (Scale):"))
      cat("\n\n")

      if (length(x$alpha) == 1) {
         tmp <- capture.output(.print.vector(res.table))
      } else {
         tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
      }
      .print.table(tmp, mstyle)

   }

   cat("\n")
   if (signif.legend)
      cat(mstyle$legend("---\nSignif. codes: "), mstyle$legend(attr(signif, "legend")))
   cat("\n\n")

   invisible()

}
