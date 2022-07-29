print.permutest.rma.uni <- function(x, digits=x$digits, signif.stars=getOption("show.signif.stars"), signif.legend=signif.stars, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="permutest.rma.uni")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   .space()

   ddd <- list(...)

   .chkdots(ddd, c("num", "legend"))

   if (is.null(ddd$legend)) {
      legend <- TRUE
   } else {
      legend <- .isTRUE(ddd$legend)
   }

   footsym <- .get.footsym()

   if (!x$int.only) {
      if (inherits(x, "permutest.rma.ls")) {
         cat(mstyle$section(paste0("Test of Location Coefficients (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):", ifelse(x$skip.beta, "", footsym[1]))))
      } else {
         cat(mstyle$section(paste0("Test of Moderators (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):", ifelse(x$skip.beta, "", footsym[1]))))
      }
      cat("\n")
      if (is.element(x$test, c("knha","adhoc","t"))) {
         cat(mstyle$result(paste0("F(df1 = ", x$QMdf[1], ", df2 = ", round(x$QMdf[2], 2), ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
      } else {
         cat(mstyle$result(paste0("QM(df = ", x$QMdf[1], ") = ", .fcf(x$QM, digits[["test"]]), ", p-val ", .pval(x$QMp, digits[["pval"]], showeq=TRUE, sep=" "))))
      }
      cat("\n\n")
   }

   if (is.element(x$test, c("knha","adhoc","t"))) {
      res.table <- data.frame(estimate=.fcf(c(x$beta), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), tval=.fcf(x$zval, digits[["test"]]), df=round(x$ddf,2), "pval"=.pval(x$pval, digits[["pval"]]), ci.lb=.fcf(x$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
      if (!x$skip.beta)
         res.table <- .addfootsym(res.table, 5, footsym[1])
      if (x$permci)
         res.table <- .addfootsym(res.table, 6:7, footsym[1])
   } else {
      res.table <- data.frame(estimate=.fcf(c(x$beta), digits[["est"]]), se=.fcf(x$se, digits[["se"]]), zval=.fcf(x$zval, digits[["test"]]), "pval"=.pval(x$pval, digits[["pval"]]), ci.lb=.fcf(x$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$ci.ub, digits[["ci"]]), stringsAsFactors=FALSE)
      if (!x$skip.beta)
         res.table <- .addfootsym(res.table, 4, footsym[1])
      if (x$permci)
         res.table <- .addfootsym(res.table, 5:6, footsym[1])
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

   if (inherits(x, "permutest.rma.ls")) {
      cat(mstyle$section("Model Results (Location):"))
   } else {
      cat(mstyle$section("Model Results:"))
   }
   cat("\n\n")
   if (x$int.only) {
      tmp <- capture.output(.print.vector(res.table))
   } else {
      tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
   }
   #tmp[1] <- paste0(tmp[1], "\u200b")
   .print.table(tmp, mstyle)

   if (inherits(x, "permutest.rma.ls")) {

      cat("\n")

      if (!x$Z.int.only) {
         cat(mstyle$section(paste0("Test of Scale Coefficients (coefficient", ifelse(x$m.alpha == 1, " ", "s "), .format.btt(x$att),"):", ifelse(x$skip.alpha, "", footsym[1]))))
         cat("\n")
         if (is.element(x$test, c("knha","adhoc","t"))) {
            cat(mstyle$result(paste0("F(df1 = ", x$QSdf[1], ", df2 = ", round(x$QSdf[2], 2), ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits=digits[["pval"]], showeq=TRUE, sep=" "))))
         } else {
            cat(mstyle$result(paste0("QS(df = ", x$QSdf[1], ") = ", .fcf(x$QS, digits[["test"]]), ", p-val ", .pval(x$QSp, digits[["pval"]], showeq=TRUE, sep=" "))))
         }
         cat("\n\n")
      }

      if (is.element(x$test, c("knha","adhoc","t"))) {
         res.table <- data.frame(estimate=.fcf(c(x$alpha), digits[["est"]]), se=.fcf(x$se.alpha, digits[["se"]]), tval=.fcf(x$zval.alpha, digits[["test"]]), df=round(x$ddf.alpha,2), "pval"=.pval(x$pval.alpha, digits[["pval"]]), ci.lb=.fcf(x$ci.lb.alpha, digits[["ci"]]), ci.ub=.fcf(x$ci.ub.alpha, digits[["ci"]]), stringsAsFactors=FALSE)
         if (!x$skip.alpha)
            res.table <- .addfootsym(res.table, 5, footsym[1])
      } else {
         res.table <- data.frame(estimate=.fcf(c(x$alpha), digits[["est"]]), se=.fcf(x$se.alpha, digits[["se"]]), zval=.fcf(x$zval.alpha, digits[["test"]]), "pval"=.pval(x$pval.alpha, digits[["pval"]]), ci.lb=.fcf(x$ci.lb.alpha, digits[["ci"]]), ci.ub=.fcf(x$ci.ub.alpha, digits[["ci"]]), stringsAsFactors=FALSE)
         if (!x$skip.alpha)
            res.table <- .addfootsym(res.table, 4, footsym[1])
      }
      rownames(res.table) <- rownames(x$alpha)
      signif <- symnum(x$pval.alpha, corr=FALSE, na=FALSE, cutpoints=c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
      if (signif.stars) {
         res.table <- cbind(res.table, signif)
         colnames(res.table)[ncol(res.table)] <- ""
      }

      if (.isTRUE(ddd$num))
         rownames(res.table) <- paste0(seq_len(nrow(res.table)), ") ", rownames(res.table))

      if (x$Z.int.only)
         res.table <- res.table[1,]

      cat(mstyle$section("Model Results (Scale):"))
      cat("\n\n")
      if (x$Z.int.only) {
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

   if (legend) {
      cat("\n")
      if (inherits(x, "permutest.rma.ls")) {
         cat(mstyle$legend(paste0(footsym[2], " p-values based on permutation testing")))
      } else {
         cat(mstyle$legend(paste0(footsym[2], " p-value", ifelse(x$int.only, "", "s"), ifelse(x$permci, " and CI bounds", ""), " based on permutation testing")))
      }
      cat("\n")
   }

   .space()

   invisible()

}
