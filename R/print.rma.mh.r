print.rma.mh <- function(x, digits, showfit=FALSE, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "rma.mh"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"rma.mh\"."))

   if (missing(digits))
      digits <- x$digits

   cat("\n")

   cat(mstyle$section("Fixed-Effects Model"))
   cat(mstyle$section(paste0(" (k = ", x$k, ")")))

   if (showfit) {
      cat("\n")
      if (anyNA(x$fit.stats$ML)) {
         fs <- x$fit.stats$ML
      } else {
         fs <- c(formatC(x$fit.stats$ML, digits=digits, format="f"))
      }
      names(fs) <- c("logLik", "deviance", "AIC", "BIC", "AICc")
      cat("\n")
      tmp <- capture.output(print(fs, quote=FALSE, print.gap=2))
      .print.table(tmp, mstyle)
      cat("\n")
   } else {
      cat("\n\n")
   }

   if (!is.na(x$QE)) {
      cat(mstyle$section("Test for Heterogeneity:"), "\n")
      cat(mstyle$result(paste0("Q(df = ", ifelse(x$k.yi-1 >= 0, x$k.yi-1, 0), ") = ", formatC(x$QE, digits=digits, format="f"), ", p-val ", .pval(x$QEp, digits=digits, showeq=TRUE, sep=" "))))
   }

   if (is.element(x$measure, c("OR","RR","IRR"))) {

      res.table     <- c(x$beta, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)
      res.table.exp <- c(exp(x$beta), exp(x$ci.lb), exp(x$ci.ub))

      if (!is.na(x$beta)) {
         res.table    <- formatC(res.table, digits=digits, format="f")
         res.table[4] <- .pval(x$pval, digits=digits)
      }

      if (!is.na(x$beta))
         res.table.exp <- formatC(res.table.exp, digits=digits, format="f")

      names(res.table)     <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
      names(res.table.exp) <- c("estimate", "ci.lb", "ci.ub")

      cat("\n\n")
      cat(mstyle$section("Model Results (log scale):"))
      cat("\n\n")
      tmp <- capture.output(.print.vector(res.table))
      .print.table(tmp, mstyle)

      cat("\n")
      cat(mstyle$section(paste0("Model Results (", x$measure, " scale):")))
      cat("\n\n")
      tmp <- capture.output(.print.vector(res.table.exp))
      .print.table(tmp, mstyle)
      cat("\n")

      if (x$measure == "OR") {
         MH <- ifelse(is.na(x$MH), NA, formatC(x$MH, digits=digits, format="f"))
         TA <- ifelse(is.na(x$TA), NA, formatC(x$TA, digits=digits, format="f"))
         if (is.na(MH) && is.na(TA)) {
            width <- 1
         } else {
            width <- max(nchar(MH), nchar(TA), na.rm=TRUE)
         }
         if (is.na(MH)) {
            cat(mstyle$text("Cochran-Mantel-Haenszel Test:    "))
            cat(mstyle$result("test value not computable for these data"))
            cat("\n")
         } else {
            cat(mstyle$text("Cochran-Mantel-Haenszel Test:    "))
            cat(mstyle$result(paste0("CMH = ", formatC(MH, width=width), ", df = 1,", paste(rep(" ", nchar(x$k.pos)-1, collapse="")), " p-val ", .pval(x$MHp, digits=digits, showeq=TRUE, sep=" ", add0=TRUE))))
            cat("\n")
         }
         if (is.na(TA)) {
            cat(mstyle$text("Tarone's Test for Heterogeneity: "))
            cat(mstyle$result("test value not computable for these data"))
         } else {
            cat(mstyle$text("Tarone's Test for Heterogeneity: "))
            cat(mstyle$result(paste0("X^2 = ", formatC(TA, width=width), ", df = ", x$k.pos-1, ", p-val ", .pval(x$TAp, digits=digits, showeq=TRUE, sep=" ", add0=TRUE))))
         }
         cat("\n\n")
      }

      if (x$measure == "IRR") {
         if (is.na(x$MH)) {
            cat(mstyle$text("Mantel-Haenszel Test: "))
            cat(mstyle$result("test value not computable for these data"))
         } else {
            cat(mstyle$text("Mantel-Haenszel Test: "))
            cat(mstyle$result(paste0("MH = ", formatC(x$MH, digits, format="f"), ", df = 1, p-val ", .pval(x$MHp, digits=digits, showeq=TRUE, sep=" "))))
         }
         cat("\n\n")
      }

   } else {

      res.table <- c(x$beta, x$se, x$zval, x$pval, x$ci.lb, x$ci.ub)

      if (!is.na(x$beta)) {
         res.table    <- formatC(res.table, digits=digits, format="f")
         res.table[4] <- .pval(x$pval, digits=digits)
      }

      names(res.table) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

      cat("\n\n")
      cat(mstyle$section("Model Results:"))
      cat("\n\n")
      tmp <- capture.output(.print.vector(res.table))
      .print.table(tmp, mstyle)
      cat("\n")

   }

   invisible()

}
