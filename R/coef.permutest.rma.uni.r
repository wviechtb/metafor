coef.permutest.rma.uni <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   .chkclass(class(object), must="permutest.rma.uni")

   x <- object

   if (is.element(x$test, c("knha","adhoc","t"))) {
      res.table <- data.frame(estimate=x$beta, se=x$se, tval=x$zval, df=x$ddf, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   } else {
      res.table <- data.frame(estimate=x$beta, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   }

   if (inherits(x, "permutest.rma.ls")) {

      if (is.element(x$test, c("knha","adhoc","t"))) {
         res.table.alpha <- data.frame(estimate=x$alpha, se=x$se.alpha, tval=x$zval.alpha, df=x$ddf.alpha, pval=x$pval.alpha, ci.lb=x$ci.lb.alpha, ci.ub=x$ci.ub.alpha)
      } else {
         res.table.alpha <- data.frame(estimate=x$alpha, se=x$se.alpha, zval=x$zval.alpha, pval=x$pval.alpha, ci.lb=x$ci.lb.alpha, ci.ub=x$ci.ub.alpha)
      }

      res.table <- list(beta=res.table, alpha=res.table.alpha)

   }

   return(res.table)

}
