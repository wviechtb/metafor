print.hc.rma.uni <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="hc.rma.uni")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   res.table <- data.frame(method   = c(x$method.rma, x$method),
                           tau2     = .fcf(c(x$tau2.rma, x$tau2), digits[["var"]]),
                           estimate = .fcf(c(x$beta.rma, x$beta), digits[["est"]]),
                           se       = .fcf(c(x$se.rma, x$se), digits[["se"]]),
                           ci.lb    = .fcf(c(x$ci.lb.rma, x$ci.lb), digits[["ci"]]),
                           ci.ub    = .fcf(c(x$ci.ub.rma, x$ci.ub), digits[["ci"]]), stringsAsFactors=FALSE)

   if (is.na(res.table$se[1]))
      res.table$se <- NULL

   rownames(res.table) <- c("rma", "hc")

   if (!exists(".rmspace"))
      cat("\n")

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   if (!exists(".rmspace"))
      cat("\n")

   invisible(res.table)

}
