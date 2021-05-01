coef.permutest.rma.uni <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="permutest.rma.uni")

   x <- object

   if (is.element(x$test, c("knha","adhoc","t"))) {
      res.table <- data.frame(estimate=x$beta, se=x$se, tval=x$zval, df=x$ddf, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   } else {
      res.table <- data.frame(estimate=x$beta, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   }

   return(res.table)

}
