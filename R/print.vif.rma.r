print.vif.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(x, "vif.rma"))
      stop(mstyle$stop("Argument 'x' must be an object of class \"vif.rma\"."))

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   if (!exists(".rmspace"))
      cat("\n")

   if (is.null(x$gvif)) {

      if (x$table) {

         res.table <- cbind(estimate=.fcf(x$vif$estimate, digits[["est"]]), se=.fcf(x$vif$se, digits[["se"]]), zval=.fcf(x$vif$zval, digits[["test"]]), "pval"=.pval(x$vif$pval, digits[["pval"]]), ci.lb=.fcf(x$vif$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$vif$ci.ub, digits[["ci"]]), vif=.fcf(x$vif$vif, digits[["est"]]), sif=.fcf(x$vif$sif, digits[["est"]]))
         rownames(res.table) <- rownames(x$vif)
         if (is.element(x$test, c("knha","adhoc","t")))
            colnames(res.table)[3] <- "tval"

         ddd <- list(...)

         if (.isTRUE(ddd$num))
            rownames(res.table) <- paste0(1:nrow(res.table), ") ", rownames(res.table))

         tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
         .print.table(tmp, mstyle)

      } else {

         print(.fcf(x$vif, digits[["est"]]), quote=FALSE, right=TRUE)

      }

   } else {

      cat(mstyle$section(paste0("Collinearity Diagnostics (coefficient", ifelse(x$m == 1, " ", "s "), .format.btt(x$btt),"):\n")))
      cat(mstyle$result(paste0("GVIF = ", .fcf(x$gvif, digits[["est"]]), ", GSIF = ", .fcf(x$gsif, digits[["est"]]), "\n")))

   }

   if (!exists(".rmspace"))
      cat("\n")

   invisible()

}
