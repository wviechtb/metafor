coef.summary.rma <- function(object, ...) {

   mstyle <- .get.mstyle()

   .chkclass(class(object), must="summary.rma")

   ddd <- list(...)

   x <- object

   if (is.element(x$test, c("knha","adhoc","t"))) {
      res.table <- data.frame(estimate=x$beta, se=x$se, tval=x$zval, df=x$ddf, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   } else {
      res.table <- data.frame(estimate=x$beta, se=x$se, zval=x$zval, pval=x$pval, ci.lb=x$ci.lb, ci.ub=x$ci.ub)
   }

   if (isTRUE(ddd$type=="beta"))
      return(res.table)

   if (inherits(x, "rma.ls")) {
      res.table <- list(beta=res.table)
      if (is.element(x$test, c("knha","adhoc","t"))) {
         res.table$alpha <- data.frame(estimate=x$alpha, se=x$se.alpha, tval=x$zval.alpha, df=x$ddf.alpha, pval=x$pval.alpha, ci.lb=x$ci.lb.alpha, ci.ub=x$ci.ub.alpha)
      } else {
         res.table$alpha <- data.frame(estimate=x$alpha, se=x$se.alpha, zval=x$zval.alpha, pval=x$pval.alpha, ci.lb=x$ci.lb.alpha, ci.ub=x$ci.ub.alpha)
      }
      if (isTRUE(ddd$type=="alpha"))
         return(res.table$alpha)
   }

   if (inherits(x, "rma.uni.selmodel")) {
      res.table <- list(beta=res.table)
      res.table$delta <- data.frame(estimate=x$delta, se=x$se.delta, zval=x$zval.delta, pval=x$pval.delta, ci.lb=x$ci.lb.delta, ci.ub=x$ci.ub.delta)
      if (length(x$delta) == 1L) {
         rownames(res.table$delta) <- "delta"
      } else {
         rownames(res.table$delta) <- paste0("delta.", seq_along(x$delta))
      }
      if (isTRUE(ddd$type=="delta"))
         return(res.table$delta)
   }

   return(res.table)

}
