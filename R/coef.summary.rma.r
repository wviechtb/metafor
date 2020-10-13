coef.summary.rma <- function(object, ...) {

   mstyle <- .get.mstyle("crayon" %in% .packages())

   if (!inherits(object, "summary.rma"))
      stop(mstyle$stop("Argument 'object' must be an object of class \"summary.rma\"."))

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

   if (inherits(x, "rma.uni.selmodel")) {
      res.table <- list(beta=res.table)
      res.table$delta <- data.frame(estimate=x$delta, se=x$se.delta, zval=x$zval.delta, pval=x$pval.delta, ci.lb=x$ci.lb.delta, ci.ub=x$ci.ub.delta)
      if (length(x$delta) == 1L) {
         rownames(res.table$delta) <- "delta"
      } else {
         rownames(res.table$delta) <- paste0("delta.", 1:length(x$delta))
      }
   }

   return(res.table)

}
