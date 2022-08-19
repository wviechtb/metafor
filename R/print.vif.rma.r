print.vif.rma <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="vif.rma")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   ddd <- list(...)

   .chkdots(ddd, c("num"))

   .space()

   if (!is.null(x$alpha)) {

      cat(mstyle$section(paste0("Location Coefficients:\n")))

      print(x[[1]], digits=digits, ...)

      .space(FALSE)

      cat(mstyle$section(paste0("Scale Coefficients:\n")))

      print(x[[2]], digits=digits, ...)

   } else {

      if (isTRUE(x$bttspec) || isTRUE(x$attspec)) {

         if (length(x$vif) == 1L) {

            if (x$vif[[1]]$m == 1) {
               cat(mstyle$section(paste0("Collinearity Diagnostics (coefficient ", x$vif[[1]]$coefs,"):\n")))
               cat(mstyle$result(paste0("VIF = ", .fcf(x$vif[[1]]$vif, digits[["est"]]), ", SIF = ", .fcf(x$vif[[1]]$sif, digits[["est"]]))))
            } else {
               cat(mstyle$section(paste0("Collinearity Diagnostics (coefficients ", x$vif[[1]]$coefs,"):\n")))
               cat(mstyle$result(paste0("GVIF = ", .fcf(x$vif[[1]]$vif, digits[["est"]]), ", GSIF = ", .fcf(x$vif[[1]]$sif, digits[["est"]]))))
            }
            if (!is.null(x$sim))
               cat(mstyle$result(paste0(", prop = ", .fcf(x$prop, 2))))
            cat("\n")

         } else {

            res.table <- do.call(rbind, x$vif)
            res.table$vif <- .fcf(res.table$vif, digits[["est"]])
            res.table$sif <- .fcf(res.table$sif, digits[["est"]])

            res.table$coefname <- NULL

            if (!is.null(x$sim))
               res.table$prop <- .fcf(x$prop, 2)

            # if all btt/att specifications are numeric, remove the 'spec' column
            if (all(substr(res.table$spec, 1, 1) %in% as.character(1:9)))
                res.table$spec <- NULL

            # just use numbers for row names
            rownames(res.table) <- NULL

            tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=1))
            .print.table(tmp, mstyle)

         }

      } else {

         vifs <- sapply(x$vif, function(x) x$vif)
         sifs <- sapply(x$vif, function(x) x$sif)

         if (is.null(x$table)) {

            if (is.null(x$sim)) {
               tmp <- .fcf(vifs, digits[["est"]])
               tmp <- capture.output(.print.vector(tmp))
               .print.table(tmp, mstyle)
            } else {
               res.table <- data.frame(vif=vifs)
               res.table$prop <- .fcf(x$prop, 2)
               res.table$vif  <- .fcf(res.table$vif, digits[["est"]])
               tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=2))
               .print.table(tmp, mstyle)
            }

         } else {

            if (length(vifs) != length(x$table$estimate)) {
               vifs <- c(NA, vifs)
               sifs <- c(NA, sifs)
               x$prop <- c(NA, x$prop)
            }

            if (is.element(x$test, c("knha","adhoc","t"))) {
               res.table <- data.frame(estimate=.fcf(x$table$estimate, digits[["est"]]), se=.fcf(x$table$se, digits[["se"]]), tval=.fcf(x$table$tval, digits[["test"]]), df=round(x$table$df,2), "pval"=.pval(x$table$pval, digits[["pval"]]), ci.lb=.fcf(x$table$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$table$ci.ub, digits[["ci"]]), vif=.fcf(vifs, digits[["est"]]), sif=.fcf(sifs, digits[["est"]]), stringsAsFactors=FALSE)
            } else {
               res.table <- data.frame(estimate=.fcf(x$table$estimate, digits[["est"]]), se=.fcf(x$table$se, digits[["se"]]), zval=.fcf(x$table$zval, digits[["test"]]), "pval"=.pval(x$table$pval, digits[["pval"]]), ci.lb=.fcf(x$table$ci.lb, digits[["ci"]]), ci.ub=.fcf(x$table$ci.ub, digits[["ci"]]), vif=.fcf(vifs, digits[["est"]]), sif=.fcf(sifs, digits[["est"]]), stringsAsFactors=FALSE)
            }
            rownames(res.table) <- rownames(x$table)

            if (!is.null(x$sim))
               res.table$prop <- .fcf(x$prop, 2)

            if (.isTRUE(ddd$num))
               rownames(res.table) <- paste0(seq_len(nrow(res.table)), ") ", rownames(res.table))

            tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE, print.gap=1))
            .print.table(tmp, mstyle)

         }

      }

      .space()

   }

   invisible()

}
