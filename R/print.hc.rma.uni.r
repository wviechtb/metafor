print.hc.rma.uni <- function(x, digits=x$digits, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(x), must="hc.rma.uni")

   digits <- .get.digits(digits=digits, xdigits=x$digits, dmiss=FALSE)

   res.table <- data.frame(method   = c(x$method.rma, x$method),
                           tau2     = fmtx(c(x$tau2.rma, x$tau2), digits[["var"]]),
                           estimate = fmtx(c(x$beta.rma, x$beta), digits[["est"]]),
                           se       = fmtx(c(x$se.rma, x$se), digits[["se"]]),
                           ci.lb    = fmtx(c(x$ci.lb.rma, x$ci.lb), digits[["ci"]]),
                           ci.ub    = fmtx(c(x$ci.ub.rma, x$ci.ub), digits[["ci"]]), stringsAsFactors=FALSE)

   if (is.na(x$se[1]))
      res.table$se <- NULL

   rownames(res.table) <- c("rma", "hc")

   .space()

   tmp <- capture.output(print(res.table, quote=FALSE, right=TRUE))
   .print.table(tmp, mstyle)

   .space()

   invisible(res.table)

}
