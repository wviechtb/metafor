print.rma.uni <- function(x, digits, showfit=FALSE, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="rma.uni")

   if (missing(digits)) {
      digits <- .get.digits(xdigits=x$digits, dmiss=TRUE)
   } else {
      digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)
   }

   footsym <- .get.footsym()

   ddd <- list(...)

   .chkdots(ddd, c("num", "legend"))

   if (is.null(ddd$legend)) {
      legend <- ifelse(inherits(x, "robust.rma"), TRUE, FALSE)
   } else {
      legend <- .isTRUE(ddd$legend)
   }

   if (inherits(x, "rma.uni.trimfill")) {
      .space()
      cat(mstyle$text(paste0("Estimated number of missing studies on the ", x$side, " side: ")))
      cat(mstyle$result(paste0(x$k0, " (SE = ", ifelse(is.na(x$se.k0), NA, .fcf(x$se.k0, digits[["se"]])), ")")))
      cat("\n")
      if (x$k0.est == "R0") {
         cat(mstyle$text(paste0("Test of H0: no missing studies on the ", x$side, " side:     ")))
         cat(paste0(rep(" ", nchar(x$k0)), collapse=""))
         cat(mstyle$result(paste0("p-val ", .pval(x$p.k0, digits[["pval"]], showeq=TRUE, sep=" "))))
         cat("\n")
      }
      .space(FALSE)
   }

   .space()

   if (x$model == "rma.ls") {

      cat(mstyle$section("Location-Scale Model"))
      cat(mstyle$section(paste0(" (k = ", x$k, "; ")))
      if (x$tau2.fix) {
         cat(mstyle$section("user-specified tau^2 value)"))
      } else {
         cat(mstyle$section(paste0("tau^2 estimator: ", x$method, ")")))
      }

   } else {

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
         if (x$tau2.fix) {
            cat(mstyle$section("user-specified tau^2 value)"))
         } else {
            cat(mstyle$section(paste0("tau^2 estimator: ", x$method, ")")))
         }
      }

   }

   cat("\n")

   if (showfit) {
      if (x$method == "REML") {
         fs <- .fcf(x$fit.stats$REML, digits[["fit"]])
      } else {
         fs <- .fcf(x$fit.stats$ML, digits[["fit"]])
      }
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      #tmp[1] <- paste0(tmp[1], "\u200b")
      .print.table(tmp, mstyle)
   }

   cat("\n")

   if (x$model == "rma.uni" || x$model == "rma.uni.selmodel") {

      if (!is.element(x$method, c("FE","EE","CE"))) {
         if (x$int.only) {
            cat(mstyle$text("tau^2 (estimated amount of total heterogeneity): "))
            cat(mstyle$result(paste0(.fcf(x$tau2, ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits[["var"]])), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , .fcf(x$se.tau2, digits[["sevar"]]), ")")))))
            cat("\n")
            cat(mstyle$text("tau (square root of estimated tau^2 value):      "))
            cat(mstyle$result(paste0(ifelse(x$tau2>=0, .fcf(sqrt(x$tau2), ifelse(x$tau2 <= .Machine$double.eps*10,0,digits[["var"]])), NA))))
            cat("\n")
         } else {
            cat(mstyle$text("tau^2 (estimated amount of residual heterogeneity):     "))
            cat(mstyle$result(paste0(.fcf(x$tau2, ifelse(abs(x$tau2) <= .Machine$double.eps*10,0,digits[["var"]])), ifelse(is.na(x$se.tau2), "", paste0(" (SE = " , .fcf(x$se.tau2, digits[["sevar"]]), ")")))))
            cat("\n")
            cat(mstyle$text("tau (square root of estimated tau^2 value):             "))
            cat(mstyle$result(paste0(ifelse(x$tau2>=0, .fcf(sqrt(x$tau2), ifelse(x$tau2 <= .Machine$double.eps*10,0,digits[["var"]])), NA))))
            cat("\n")
         }
      }

      if (x$int.only) {
         if (!is.na(x$I2)) {
            cat(mstyle$text("I^2 (total heterogeneity / total variability):   "))
            cat(mstyle$result(paste0(ifelse(is.na(x$I2), NA, .fcf(x$I2, 2)), "%")))
            cat("\n")
         }
         if (!is.na(x$H2)) {
            cat(mstyle$text("H^2 (total variability / sampling variability):  "))
            cat(mstyle$result(paste0(ifelse(is.na(x$H2), NA, .fcf(x$H2, 2)))))
            cat("\n")
         }
      } else {
         if (!is.na(x$I2)) {
            cat(mstyle$text("I^2 (residual heterogeneity / unaccounted variability): "))
            cat(mstyle$result(paste0(ifelse(is.na(x$I2), NA, .fcf(x$I2, 2)), "%")))
            cat("\n")
         }
         if (!is.na(x$H2)) {
            cat(mstyle$text("H^2 (unaccounted variability / sampling variability):   "))
            cat(mstyle$result(paste0(ifelse(is.na(x$H2), NA, .fcf(x$H2, 2)))))
            cat("\n")
         }
      }

      if (!x$int.only && !is.null(x$R2)) {
         cat(mstyle$text("R^2 (amount of heterogeneity accounted for):            "))
         cat(mstyle$result(paste0(ifelse(is.na(x$R2), NA, .fcf(x$R2, 2)), "%")))
         cat("\n")
      }

      if (!is.element(x$method, c("FE","EE","CE")) || !is.na(x$I2) || !is.na(x$H2) || !is.null(x$R2))
         cat("\n")

   }

   if (inherits(x, "rma.gen")) {
      cat(mstyle$section("Parameter Estimates:"))
      cat("\n\n")
      res.table <- data.frame(as.list(.fcf(x$pars, digits[["var"]])))
      colnames(res.table) <- names(x$pars)
      res.table <- res.table[1,,drop=FALSE]
      tmp <- capture.output(.print.vector(res.table))
      .print.table(tmp, mstyle)
      cat("\n")
   }

   if (!is.na(x$QE)) {
      if (x$int.only) {
         cat(mstyle$section("Test for Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("Q(df = ", x$k-x$p, ") = ", .fcf(x$QE, digits[["test"]]), ", p-val ", .pval(x$QEp, digits[["pval"]], showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$section("Test for Residual Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("QE(df = ", x$k-x$p, ") = ", .fcf(x$QE, digits[["test"]]), ", p-val ", .pval(x$QEp, digits[["pval"]], showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   if (x$model == "rma.uni.selmodel" && !is.na(x$LRT.tau2)) {
      if (x$int.only) {
         cat(mstyle$section("Test for Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("LRT(df = 1) = ", .fcf(x$LRT.tau2, digits[["test"]]), ", p-val ", .pval(x$LRTp.tau2, digits[["pval"]], showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$section("Test for Residual Heterogeneity:"))
         cat("\n")
         cat(mstyle$result(paste0("LRT(df = 1) = ", .fcf(x$LRT.tau2, digits[["test"]]), ", p-val ", .pval(x$LRTp.tau2, digits[["pval"]], showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   if (inherits(x, "robust.rma")) {

      cat(mstyle$text("Number of estimates:   "))
      cat(mstyle$result(x$k))
      cat("\n")
      cat(mstyle$text("Number of clusters:    "))
      cat(mstyle$result(x$n))
      cat("\n")

      cat(mstyle$text("Estimates per cluster: "))
      if (all(x$tcl[1] == x$tcl)) {
         cat(mstyle$result(x$tcl[1]))
      } else {
         cat(mstyle$result(paste0(min(x$tcl), "-", max(x$tcl), " (mean: ", .fcf(mean(x$tcl), digits=2), ", median: ", round(median(x$tcl), digits=2), ")")))
      }
      cat("\n\n")

   }

   if (x$p > 1L && !is.na(x$QM)) {
      if (x$model == "rma.ls") {
         cat(mstyle$section(paste0("Test of Location Coefficients (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):")))
      } else {
         cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):", ifelse(inherits(x, "robust.rma"), footsym[1], ""))))
      }
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$QMdf[1], ", df2 = ", round(x$QMdf[2], 2), ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits[["pval"]], showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QM(df = ", x$QMdf[1], ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits[["pval"]], showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   if (is.element(x$test, c("knha","adhoc","t"))) {
      res.table <- data.frame(estimate=.fcf(c(x$beta), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), tval=.fcf(x$zval, digits[["test"]]), df=round(x$ddf,2), pval=.pval(x$pval, digits[["pval"]]), ci.lb=.fcf(x$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
      if (inherits(x, "robust.rma"))
         colnames(res.table)[c(2:7)] <- paste0(colnames(res.table)[c(2:7)], footsym[1])
   } else {
      res.table <- data.frame(estimate=.fcf(c(x$beta), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), pval=.pval(x$pval, digits[["pval"]]), ci.lb=.fcf(x$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
   }
   rownames(res.table) <- rownames(x$beta)
   signif <- symnum(x$pval, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
   if (signif.stars) {
      res.table <- cbind(res.table, signif)
      colnames(res.table)[ncol(res.table)] <- ""
   }

   if (.isTRUE(ddd$num))
      rownames(res.table) <- paste0(seq_len(nrow(res.table)), ") ", rownames(res.table))

   if (x$int.only)
      res.table <- res.table[1,]

   if (x$model == "rma.uni" || x$model == "rma.uni.selmodel") {
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
   #tmp[1] <- paste0(tmp[1], "\u200b")
   .print.table(tmp, mstyle)

   if (x$model == "rma.ls") {

      if (x$q > 1L && !is.na(x$QS)) {
         cat("\n")
         cat(mstyle$section(paste0("Test of Scale Coefficients (coefficient", ifelse(x$m.alpha == 1, " ", "s "), .format.btt(x$att),"):")))
         cat("\n")
         if (is.element(x$test, c("knha","adhoc","t"))) {
            cat(mstyle$result(paste0("F(df1 = ", x$QSdf[1], ", df2 = ", round(x$QSdf[2], 2), ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits[["pval"]], showeq=TRUE, sep=" "))))
         } else {
            cat(mstyle$result(paste0("QS(df = ", x$QSdf[1], ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits[["pval"]], showeq=TRUE, sep=" "))))
         }
         cat("\n")
      }

      if (is.element(x$test, c("knha","adhoc","t"))) {
         res.table <- data.frame(estimate=.fcf(c(x$alpha), digits[["est"]]), se=.fcf(x$se.alpha, digits[["se"]]), tval=.fcf(x$zval.alpha, digits[["test"]]), df=round(x$ddf.alpha, 2), pval=.pval(x$pval.alpha, digits[["pval"]]), ci.lb=.fcf(x$ci.lb.alpha, digits[["ci"]]), ci.ub=.fcf(x$ci.ub.alpha, digits[["ci"]]), stringsAsFactors=FALSE)
      } else {
         res.table <- data.frame(estimate=.fcf(c(x$alpha), digits[["est"]]), se=.fcf(x$se.alpha, digits[["se"]]), zval=.fcf(x$zval.alpha, digits[["test"]]), pval=.pval(x$pval.alpha, digits[["pval"]]), ci.lb=.fcf(x$ci.lb.alpha, digits[["ci"]]), ci.ub=.fcf(x$ci.ub.alpha, digits[["ci"]]), stringsAsFactors=FALSE)
      }
      rownames(res.table) <- rownames(x$alpha)
      signif <- symnum(x$pval.alpha, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
      if (signif.stars) {
         res.table <- cbind(res.table, signif)
         colnames(res.table)[ncol(res.table)] <- ""
      }

      for (j in seq_len(nrow(res.table))) {
         res.table[j, is.na(res.table[j,])]  <- ifelse(x$alpha.fix[j], "---", "NA")
         res.table[j, res.table[j,] == "NA"] <- ifelse(x$alpha.fix[j], "---", "NA")
      }

      if (.isTRUE(ddd$num))
         rownames(res.table) <- paste0(seq_len(nrow(res.table)), ") ", rownames(res.table))

      if (length(x$alpha) == 1L)
         res.table <- res.table[1,]

      cat("\n")
      cat(mstyle$section("Model Results (Scale):"))
      cat("\n\n")

      if (length(x$alpha) == 1L) {
         tmp <- capture.output(.print.vector(res.table))
      } else {
         tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
      }
      #tmp[1] <- paste0(tmp[1], "\u200b")
      .print.table(tmp, mstyle)

   }

   if (x$model == "rma.uni.selmodel") {

      if (!is.na(x$LRT)) {
         cat("\n")
         cat(mstyle$section("Test for Selection Model Parameters:"))
         cat("\n")
         cat(mstyle$result(paste0("LRT(df = ", x$LRTdf, ") = ", .fcf(x$LRT, digits[["test"]]), ", p-val ", .pval(x$LRTp, digits[["pval"]], showeq=TRUE, sep=" "))))
         cat("\n")
      }

      res.table <- data.frame(estimate=.fcf(c(x$delta), digits[["est"]]), se=.fcf(x$se.delta, digits[["se"]]), zval=.fcf(x$zval.delta, digits[["test"]]), pval=.pval(x$pval.delta, digits[["pval"]]), ci.lb=.fcf(x$ci.lb.delta, digits[["ci"]]), ci.ub=.fcf(x$ci.ub.delta, digits[["ci"]]), stringsAsFactors=FALSE)

      if (x$type == "stepfun") {
         rownames(res.table) <- rownames(x$ptable)
         res.table <- cbind(k=x$ptable$k, res.table)
      } else {
         rownames(res.table) <- paste0("delta.", seq_along(x$delta))
      }
      #if (x$test == "t")
      #   colnames(res.table)[3] <- "tval"
      signif <- symnum(x$pval.delta, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
      if (signif.stars) {
         res.table <- cbind(res.table, signif)
         colnames(res.table)[ncol(res.table)] <- ""
      }

      for (j in seq_len(nrow(res.table))) {
         res.table[j, is.na(res.table[j,])]  <- ifelse(x$delta.fix[j], "---", "NA")
         res.table[j, res.table[j,] == "NA"] <- ifelse(x$delta.fix[j], "---", "NA")
      }

      if (length(x$delta) == 1L)
         res.table <- res.table[1,]

      cat("\n")
      cat(mstyle$section("Selection Model Results:"))
      cat("\n\n")

      if (length(x$delta) == 1L) {
         tmp <- capture.output(.print.vector(res.table))
      } else {
         tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
      }
      #tmp[1] <- paste0(tmp[1], "\u200b")
      .print.table(tmp, mstyle)

   }

   if (signif.legend || legend) {
      cat("\n")
      cat(mstyle$legend("---"))
   }

   if (signif.legend) {
      cat("\n")
      cat(mstyle$legend("Signif. codes: "), mstyle$legend(attr(signif, "legend")))
      cat("\n")
   }

   if (inherits(x, "robust.rma") && legend) {
      cat("\n")
      cat(mstyle$legend(paste0(footsym[2], " results based on cluster-robust inference (var-cov estimator: ", x$vbest)))
      if (x$robumethod == "default") {
         cat(mstyle$legend(","))
         cat("\n")
         cat(mstyle$legend(paste0("   approx. ", ifelse(x$int.only, "t-test and confidence interval", "t/F-tests and confidence intervals"), ", dfs = residual method)")))
      } else {
         if (x$coef_test == "Satterthwaite" && x$conf_test == "Satterthwaite" && x$wald_test == "HTZ") {
            cat(mstyle$legend(","))
            cat("\n")
            cat(mstyle$legend(paste0("   approx. ", ifelse(x$int.only, "t-test and confidence interval", "t/F-tests and confidence intervals"), ", dfs = Satterthwaite method)")))
         } else {
            cat(mstyle$legend(")"))
         }
      }
      cat("\n")
   }

   .space()

   invisible()

}
