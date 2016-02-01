coef.permutest.rma.uni <- function(object, ...) {

   if (!inherits(object, "permutest.rma.uni"))
      stop("Argument 'object' must be an object of class \"permutest.rma.uni\".")

   x <- object

   res.table <- data.frame(estimate=x$b, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)

   if ((is.logical(x$knha) && x$knha) || is.character(x$knha))
      colnames(res.table)[3] <- "tval"

   return(res.table)

}
