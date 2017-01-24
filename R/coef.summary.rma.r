coef.summary.rma <- function(object, ...) {

   if (!inherits(object, "summary.rma"))
      stop("Argument 'object' must be an object of class \"summary.rma\".")

   x <- object

   res.table <- data.frame(estimate=x$beta, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)

   if (is.element(x$test, c("knha","adhoc","t")))
      colnames(res.table)[3] <- "tval"

   if (inherits(x, "rma.ls")) {
      res.table <- list(beta=res.table)
      res.table$alpha <- data.frame(estimate=x$alpha, se=x$se.alpha, zval=x$zval.alpha, pval=x$pval.alpha, ci.lb=x$ci.lb.alpha, ci.ub=x$ci.ub.alpha)
      if (is.element(x$test, c("t")))
         colnames(res.table$alpha)[3] <- "tval"
   }

   return(res.table)

}
