# Note: Works with "robust.rma" objects.

coef.summary.rma <- function(object, ...) {

   if (!is.element("summary.rma", class(object)))
      stop("Argument 'object' must be an object of class \"summary.rma\".")

   x <- object

   if (is.element("robust.rma", class(x))) ### so that code below works with x$zval
      x$zval <- x$tval

   res.table <- data.frame(estimate=x$b, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)

   if ((is.logical(x$knha) && x$knha) || is.character(x$knha))
      colnames(res.table)[3] <- "tval"

   return(res.table)

}
